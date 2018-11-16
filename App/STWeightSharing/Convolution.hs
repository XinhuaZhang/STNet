module App.STWeightSharing.Convolution where

import           Control.Monad               as M
import           Control.Parallel.Strategies
import           Data.Array                  as Arr
import           Data.Array.Repa             as R
import           Data.Complex
import           Data.List                   as L
import           Data.Vector.Storable        as VS
import           Data.Vector.Unboxed         as VU
import           Src.Array.TwoHalfDArray
import           Src.Utils.Coordinates
import           Src.Utils.DFT
import           System.Random

{-# INLINE generateHarmonicCoefficient' #-}
generateHarmonicCoefficient' :: Double -> Int -> Complex Double
generateHarmonicCoefficient' theta freq = exp $ 0 :+ (-1) * (fromIntegral freq) * theta

{-# INLINE generateHarmonicCoefficients #-}
generateHarmonicCoefficients :: Double -> [Int] -> [Complex Double]
generateHarmonicCoefficients theta = L.map (generateHarmonicCoefficient' theta)

{-# INLINE initialDist #-}
initialDist :: Int
            -> [Int]
            -> [(Double, Int, Int)]
            -> [R.Array U DIM2 (Complex Double)]
initialDist size freqs xs =
  L.map
    (\freq ->
        fromListUnboxed (Z :. size :. size) .
        Arr.elems .
        accumArray (+) (0 :: Complex Double) ((0, 0), (size - 1, size - 1)) .
        L.map
          (\(theta, x, y) -> ((x, y), generateHarmonicCoefficient' theta freq)) $
        xs)
    freqs

{-# INLINE makeFilter #-}
makeFilter :: (R.Source s e) => R.Array s DIM3 e -> R.Array D DIM3 e
makeFilter arr =
  let (Z :. rows :. cols :. _) = extent arr
   in R.backpermute
        (extent arr)
        (\(Z :. i :. j :. k) ->
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
            in (Z :. x :. y :. k))
        arr

{-# INLINE crosscorrelation #-}
crosscorrelation ::
     DFTPlan
  -> R.Array U DIM2 (Complex Double)
  -> R.Array U DIM3 (Complex Double)
  -> IO (R.Array U DIM3 (Complex Double))
crosscorrelation plan img harmonics = do
  let (Z :. rows :. cols :. freqs) = extent harmonics
      planID2D = DFTPlanID DFT1DG ([rows, cols]) [0, 1]
      planID = DFTPlanID DFT1DG ([rows, cols, freqs]) [0, 1]
      inversePlanID = DFTPlanID IDFT1DG ([rows, cols, freqs]) [0, 1]
      imgVec = VU.convert . toUnboxed $ img
      harmonicsVec = VU.convert . toUnboxed . computeS . makeFilter $ harmonics
  imgVecF <- dftExecute plan planID2D imgVec
  harmonicsVecF <- dftExecute plan planID harmonicsVec
  let convolvedVecF =
        computeS .
        R.traverse2
          (fromUnboxed (extent harmonics) . VS.convert . VS.map conjugate $ harmonicsVecF)
          (fromUnboxed (extent img) . VS.convert $ imgVecF)
          const $
        (\f3d f2d idx@(Z :. i :. j :. k) -> f2d (Z :. i :. j) * f3d idx)
  convolvedVec <- dftExecute plan inversePlanID . VU.convert . toUnboxed $ convolvedVecF
  return . fromUnboxed (extent harmonics) . VS.convert $ convolvedVec


{-# INLINE generateDFTPlan #-}
generateDFTPlan :: DFTPlan -> R.Array U DIM3 (Complex Double) -> IO DFTPlan
generateDFTPlan plan arr = do
  let (Z :. rows :. cols :. orientations) = extent arr
  lock <- getFFTWLock
  vecTemp1 <-
    VS.fromList <$> M.replicateM (orientations * rows * cols) randomIO :: IO (VS.Vector Double)
  vecTemp2 <-
    VS.fromList <$> M.replicateM (orientations * rows * cols) randomIO :: IO (VS.Vector Double)
  (plan1, vecTemp3) <-
    dft1dGPlan
      lock
      plan
      ([rows, cols, orientations])
      [0, 1]
      (VS.zipWith mkPolar vecTemp1 vecTemp2)
  (plan2, _) <-
    idft1dGPlan lock plan1 ([rows, cols, orientations]) [0, 1] vecTemp3
  vecTemp4 <-
    VS.fromList <$> M.replicateM (rows * cols) randomIO :: IO (VS.Vector Double)
  vecTemp5 <-
    VS.fromList <$> M.replicateM (rows * cols) randomIO :: IO (VS.Vector Double)
  (plan3, _) <-
    dft1dGPlan
      lock
      plan2
      [rows, cols]
      [0, 1]
      (VS.zipWith mkPolar vecTemp4 vecTemp5)
  return plan3

{-# INLINE timeReversal #-}

timeReversal
  :: (Unbox e)
  => R.Array D DIM3 e -> R.Array U DIM3 e
timeReversal arr =
  let (Z :. _ :. _ :. nf) = extent arr
      n = div nf 2
  in computeS $
     R.backpermute
       (extent arr)
       (\(Z :. i :. j :. k) ->
           let x = nf - n
           in if k >= x
                then (Z :. i :. j :. k - x)
                else (Z :. i :. j :. k + n))
       arr
