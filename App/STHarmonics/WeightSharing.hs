{-# LANGUAGE FlexibleContexts #-}
module App.STHarmonics.WeightSharing
  ( ConvolutionTypeST(..)
  , shareWeightST
  , module Src.Utils.DFT
  ) where

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

{-# INLINE rotateST #-}

rotateST
  :: Int -> (Double, Double) -> Array U DIM3 Double -> Array U DIM3 Double
rotateST 0 _ arr = arr
rotateST n valRange arr =
  let (Z :. nf :. rows :. cols) = extent arr
      deg = 360 / (fromIntegral nf) * fromIntegral n
  in rotate25D deg valRange $
     R.backpermute
       (extent arr)
       (\(Z :. k :. i :. j) ->
          if k <= n - 1
            then (Z :. (nf - n + k) :. i :. j)
            else (Z :. (k - n) :. i :. j))
       arr

{-# INLINE convolve #-}

convolve
  :: DFTPlan
  -> ConvolutionTypeST
  -> Array U DIM3 Double
  -> Array U DIM3 Double
  -> IO (VU.Vector Double)
convolve dftPlan convolutionType arr arrG = do
  let (Z :. orientations :. rows :. cols) = extent arrG
      planID = DFTPlanID DFT1DG [orientations, rows, cols] [1, 2]
      inversePlanID = DFTPlanID IDFT1DG [orientations, rows, cols] [1, 2]
  arrF <-
    dftExecute dftPlan planID . VU.convert . VU.map (:+ 0) . toUnboxed $ arr
  arrGF <-
    dftExecute dftPlan planID . VU.convert . VU.map (:+ 0) . toUnboxed $ arrG
  let arrGF' =
        case convolutionType of
          ConvolutionST -> arrGF
          CrossCorrelationST -> VS.map (conjugate) arrGF
  vec <- dftExecute dftPlan inversePlanID (VS.zipWith (*) arrF arrGF')
  let outputArr =
        fromUnboxed (extent arrG) . VU.map magnitude . VS.convert $ vec
  return . L.foldl1' (VU.zipWith (+)) . get25DArray2DVector $ outputArr

{-# INLINE _shareWeightST #-}

_shareWeightST
  :: DFTPlan
  -> ConvolutionTypeST
  -> Int
  -> (Double, Double)
  -> Array U DIM3 Double
  -> Array U DIM3 Double
  -> IO [VU.Vector Double]
_shareWeightST dftPlan convolutionType n valRange arr arrG
  | nf == n = return []
  | otherwise = do
    let rotatedG = rotateST n valRange arrG
    x <- convolve dftPlan convolutionType arr rotatedG
    xs <- _shareWeightST dftPlan convolutionType (n + 1) valRange arr rotatedG
    return (x : xs)
  where
    (Z :. nf :. _ :. _) = extent arr

{-# INLINE shareWeightST #-}

shareWeightST
  :: DFTPlan
  -> ConvolutionTypeST
  -> (Double, Double)
  -> Array U DIM3 Double
  -> Array U DIM3 Double
  -> IO [VU.Vector Double]
shareWeightST dftPlan convolutionType valRange arr arrG =
  _shareWeightST dftPlan convolutionType 0 valRange arr arrG
