{-# LANGUAGE FlexibleContexts #-}
module App.STHarmonics.WeightSharing
  ( ConvolutionTypeST(..)
  , timeReversal
  , shareWeightST
  , makeFilter
  , module Src.Utils.DFT
  ) where

import           Control.Monad.Parallel  as MP
import           Data.Array.Repa         as R
import           Data.Complex
import           Data.List               as L
import           Data.Vector.Storable    as VS
import           Data.Vector.Unboxed     as VU
import           Src.Array.Transform
import           Src.Array.TwoHalfDArray
import           Src.Utils.DFT

data ConvolutionTypeST
  = ConvolutionST
  | CrossCorrelationST
  deriving (Read, Show)
  
data TimeReversal
  = SourceST
  | SinkST
  deriving (Read, Show)
  

{-# INLINE checkN #-}

checkN :: Int -> Int -> Int
checkN maxN n
  | n < 0 = checkN maxN (n + maxN)
  | n >= maxN = checkN maxN (n - maxN)
  | otherwise = n 

{-# INLINE rotateST #-}

rotateST :: (Double, Double)
         -> Array U DIM3 Double
         -> Int
         -> Array U DIM3 Double
rotateST _ arr 0 = arr
rotateST valRange arr n' =
  let (Z :. nf :. _ :. _) = extent arr
      n = checkN nf n'
      deg = 360 / (fromIntegral nf) * fromIntegral n
  in rotate25D deg valRange $
     R.backpermute
       (extent arr)
       (\(Z :. k :. i :. j) ->
          let x = nf - n
          in if k >= x
               then (Z :. k - x :. i :. j)
               else (Z :. k + n :. i :. j))
       arr
       

{-# INLINE timeReversal #-}

timeReversal :: Array U DIM3 Double -> Array U DIM3 Double
timeReversal arr =
  let (Z :. nf :. _ :. _) = extent arr
      n = div nf 2
  in computeS $
     R.backpermute
       (extent arr)
       (\(Z :. k :. i :. j) ->
          let x = nf - n
          in if k >= x
               then (Z :. k - x :. i :. j)
               else (Z :. k + n :. i :. j))
       arr

{-# INLINE convolve #-}

convolve
  :: DFTPlan
  -> ConvolutionTypeST
  -> Array U DIM3 Double
  -> Array U DIM2 Double
  -> IO (VU.Vector Double)
convolve dftPlan convolutionType arr3D arr2D = do
  let (Z :. orientations :. rows :. cols) = extent arr3D
      planID2D = DFTPlanID DFT1DG [rows, cols] [0, 1]
      planID3D = DFTPlanID DFT1DG [orientations, rows, cols] [1, 2]
      inversePlanID = DFTPlanID IDFT1DG [orientations, rows, cols] [1, 2]
  vec2DF <-
    dftExecute dftPlan planID2D . VU.convert . VU.map (:+ 0) . toUnboxed $ arr2D
  vec3DF <-
    dftExecute dftPlan planID3D . VU.convert . VU.map (:+ 0) . toUnboxed $ arr3D
  let vec2DF' =
        case convolutionType of
          ConvolutionST -> vec2DF
          CrossCorrelationST -> VS.map conjugate vec2DF
      arr2DF = fromUnboxed (Z :. rows :. cols) . VS.convert $ vec2DF'
      arr3DF =
        fromUnboxed (Z :. orientations :. rows :. cols) . VS.convert $ vec3DF
  vec <-
    dftExecute
      dftPlan
      inversePlanID
      (VU.convert . toUnboxed . computeS $ zipWith25D (*) arr3DF arr2DF)
  return . VU.map realPart . VS.convert $ vec

{-# INLINE shareWeightST #-}

shareWeightST
  :: DFTPlan
  -> ConvolutionTypeST
  -> (Double, Double)
  -> Array U DIM3 Double
  -> Array U DIM3 Double
  -> IO (Array U DIM3 Double)
shareWeightST dftPlan convolutionType valRange arr arrG = do
  let (Z :. nf :. _ :. _) = extent arrG
  xs <-
    MP.mapM
      (\i ->
         convolve
           dftPlan
           convolutionType
           (computeS . makeFilter . rotateST valRange arrG $ i)
           (computeS . R.slice arr $ (Z :. i :. All :. All)))
      [0 .. nf - 1]
  return . fromUnboxed (extent arrG) . L.foldl1' (VU.zipWith (+)) $ xs

{-# INLINE makeFilter #-}

makeFilter :: (R.Source s e) => R.Array s DIM3 e -> R.Array D DIM3 e 
makeFilter arr =
  let (Z :. _ :. rows :. cols) = extent arr
  in R.backpermute
       (extent arr)
       (\(Z :. k :. i :. j) ->
          let halfRows = div rows 2
              halfCols = div cols 2
              x =
                if i < halfRows
                  then i + halfRows
                  else i - halfRows
              y =
                if j < halfCols
                  then j + halfCols
                  else j - halfCols
          in (Z :. k :. x :. y))
       arr
