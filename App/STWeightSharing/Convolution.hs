{-# LANGUAGE FlexibleContexts #-}
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
  vecTemp6 <-
    VS.fromList <$> M.replicateM (orientations * rows * cols) randomIO :: IO (VS.Vector Double)
  vecTemp7 <-
    VS.fromList <$> M.replicateM (orientations * rows * cols) randomIO :: IO (VS.Vector Double)
  (plan4, vecTemp8) <-
    dft1dGPlan
      lock
      plan3
      ([rows, cols, orientations])
      [2]
      (VS.zipWith mkPolar vecTemp6 vecTemp7)
  (plan5, _) <- idft1dGPlan lock plan4 ([rows, cols, orientations]) [2] vecTemp8
  return plan5

{-# INLINE timeReversal #-}

timeReversal :: (R.Source s e, Unbox e) => R.Array s DIM3 e -> R.Array U DIM3 e
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

{-# INLINE freqDomain2S1 #-}
freqDomain2S1 ::
     (R.Source s (Complex Double))
  => Int
  -> R.Array s DIM3 (Complex Double)
  -> (R.Array U DIM3 (Complex Double))
freqDomain2S1 oris arr =
  let (Z :. rows :. cols :. n) = extent arr
      freq = div (n - 1) 2
      xs = L.map (\i -> R.slice arr (Z :. All :. All :. i)) [0 .. n - 1]
      ys =
        L.map
          (\m ->
             fromFunction (Z :. rows :. cols :. oris) $ \(Z :. _ :. _ :. i) ->
               let theta = fromIntegral i * 2 * pi / (fromIntegral oris)
                in exp $ 0 :+ (-1) * fromIntegral m * theta)
          [-freq .. freq]
   in computeS . L.foldl1' (R.zipWith (+)) $ 
      L.zipWith
        (\x y ->
           traverse2
             y
             x
             const
             (\f3d f2d idx@(Z :. i :. j :. k) -> f3d idx * f2d (Z :. i :. j)))
        xs
        ys


{-# INLINE makeFilter1D #-}
makeFilter1D :: (R.Source s e) => R.Array s DIM3 e -> R.Array D DIM3 e
makeFilter1D arr =
  let (Z :. _ :. _ :. n) = extent arr
   in R.backpermute
        (extent arr)
        (\(Z :. i :. j :. k) ->
           let half = div n 2
               z =
                 if k < half
                   then k + half
                   else k - half
            in (Z :. i :. j :. z))
        arr

{-# INLINE timeReversal1D #-}
timeReversal1D ::
     (R.Source s (Complex Double))
  => R.Array s DIM3 (Complex Double)
  -> R.Array U DIM3 (Complex Double)
timeReversal1D arr =
  let (Z :. _ :. _ :. n) = extent arr
      freq = div (n - 1) 2
   in computeS . R.traverse arr id $ \f idx@(Z :. _ :. _ :. k) ->
        let m = fromIntegral $ k - freq
         in f idx * (exp (0 :+ m * pi))

completion ::
     (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  => DFTPlan
  -> R.Array s1 DIM3 (Complex Double)
  -> R.Array s2 DIM3 (Complex Double)
  -> IO (R.Array U DIM3 (Complex Double))
completion plan arr1 arr2 = do
  let (Z :. rows :. cols :. oris) = extent arr1
      planID = DFTPlanID DFT1DG [rows, cols, oris] [2]
      inversePlanID = DFTPlanID IDFT1DG [rows, cols, oris] [2]
      vec1 = VU.convert . toUnboxed . computeS . delay $ arr1
      vec2 = VU.convert . toUnboxed . computeS . makeFilter1D $ arr2
  vec1F <- dftExecute plan planID vec1
  vec2F <- dftExecute plan planID vec2
  convolvedVec <- dftExecute plan inversePlanID $ VS.zipWith (*) vec1F vec2F
  return . fromUnboxed (extent arr1) . VS.convert $ convolvedVec
