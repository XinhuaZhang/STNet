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
import           Src.Array.UnboxedArray      as UA
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

{-# INLINE freqDomainR2S1 #-}
freqDomainR2S1 ::
     (R.Source s (Complex Double))
  => Int
  -> R.Array s DIM3 (Complex Double)
  -> (R.Array U DIM3 (Complex Double))
freqDomainR2S1 oris arr =
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


{-# INLINE makeFilterR2S1 #-}
makeFilterR2S1 :: (R.Source s e) => R.Array s DIM3 e -> R.Array D DIM3 e
makeFilterR2S1 arr =
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

{-# INLINE timeReversalR2S1 #-}
timeReversalR2S1 ::
     (R.Source s (Complex Double))
  => R.Array s DIM3 (Complex Double)
  -> R.Array D DIM3 (Complex Double)
timeReversalR2S1 arr =
  let (Z :. _ :. _ :. n) = extent arr
      freq = div (n - 1) 2
   in R.traverse arr id $ \f idx@(Z :. _ :. _ :. k) ->
        let m = fromIntegral $ k - freq
         in f idx * (exp (0 :+ m * pi))

{-# INLINE completionR2S1 #-}
completionR2S1 ::
     (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  => DFTPlan
  -> R.Array s1 DIM3 (Complex Double)
  -> R.Array s2 DIM3 (Complex Double)
  -> IO (R.Array U DIM3 (Complex Double))
completionR2S1 plan arr1 arr2 = do
  let (Z :. rows :. cols :. oris) = extent arr1
      planID = DFTPlanID DFT1DG [rows, cols, oris] [2]
      inversePlanID = DFTPlanID IDFT1DG [rows, cols, oris] [2]
      vec1 = VU.convert . toUnboxed . computeS . delay $ arr1
      vec2 = VU.convert . toUnboxed . computeS . makeFilterR2S1 $ arr2
  vec1F <- dftExecute plan planID vec1
  vec2F <- dftExecute plan planID vec2
  convolvedVec <- dftExecute plan inversePlanID $ VS.zipWith (*) vec1F vec2F
  return . fromUnboxed (extent arr1) . VS.convert $ convolvedVec


{-# INLINE generateHarmonicCoefficientR2S1RP' #-}
generateHarmonicCoefficientR2S1RP' ::
     Double -> Int -> Double -> Int -> Double -> Complex Double
generateHarmonicCoefficientR2S1RP' theta angularFreq scale radialFreq maxScale =
  exp $
  0 :+
  ((-1) *
   ((fromIntegral angularFreq) * theta +
    (fromIntegral radialFreq) * 2 * pi * scale / maxScale))

{-# INLINE generateHarmonicCoefficientsR2S1RP #-}
generateHarmonicCoefficientsR2S1RP ::
     Double -> [Int] -> Double -> [Int] -> Double -> [Complex Double]
generateHarmonicCoefficientsR2S1RP theta angularFreqs scale radialFreqs maxScale =
  L.map
    (\(t, s) -> generateHarmonicCoefficientR2S1RP' theta t scale s maxScale)
    [(t, s) | t <- angularFreqs, s <- radialFreqs]
    
{-# INLINE initialDistR2S1RP #-}
initialDistR2S1RP ::
     Int
  -> [Int]
  -> [Int]
  -> Double
  -> [(Int, Int, Double, Double)]
  -> [R.Array U DIM2 (Complex Double)]
initialDistR2S1RP size angularFreqs radialFreqs maxScale xs =
  L.map
    (\(af, rf) ->
       fromUnboxed (Z :. size :. size) .
       toUnboxedVector .
       UA.accum (+) (0 :: Complex Double) ((0, 0), (size - 1, size - 1)) .
       L.map
         (\(x, y, theta, scale) ->
            ( (x, y)
            , generateHarmonicCoefficientR2S1RP' theta af scale rf maxScale)) $
       xs)
    [(af, rf) | af <- angularFreqs, rf <- radialFreqs]


{-# INLINE freqDomainR2S1RP #-}
freqDomainR2S1RP ::
     (R.Source s (Complex Double))
  => Int
  -> Int
  -> Double
  -> R.Array s DIM4 (Complex Double)
  -> (R.Array U DIM4 (Complex Double))
freqDomainR2S1RP oris scales maxScale arr =
  let (Z :. rows :. cols :. m :. n) = extent arr
      angularFreq = div (m - 1) 2
      radialFreq = div (n - 1) 2
      xs =
        L.map
          (\(i, j) -> R.slice arr (Z :. All :. All :. i :. j))
          [(i, j) | i <- [0 .. m - 1], j <- [0 .. n - 1]]
      ys =
        L.map (fromUnboxed (extent arr)) $
        parMap
          rdeepseq
          (\(i, j) ->
             VU.concat .
             L.replicate (rows * cols) .
             toUnboxed . computeUnboxedS . fromFunction (Z :. oris :. scales) $ \(Z :. a :. b) ->
               let theta = fromIntegral a * 2 * pi / (fromIntegral oris)
                   s = fromIntegral b * maxScale / fromIntegral scales
                in exp $
                   0 :+
                   ((-1) * fromIntegral i * theta +
                    (-1) * fromIntegral j * s * 2 * pi / maxScale))
          [ (i, j)
          | i <- [-angularFreq .. angularFreq]
          , j <- [-radialFreq .. radialFreq]
          ]
   in fromUnboxed (extent arr) .
      L.foldl1' (VU.zipWith (+)) .
      parMap
        rdeepseq
        (\(x, y) ->
           toUnboxed . computeS $
           traverse2
             y
             x
             const
             (\f4d f2d idx@(Z :. i :. j :. k :. l) ->
                f4d idx * f2d (Z :. i :. j))) $
      L.zip xs ys


{-# INLINE makeFilterR2S1RP #-}
makeFilterR2S1RP :: (R.Source s e) => R.Array s DIM4 e -> R.Array D DIM4 e
makeFilterR2S1RP arr =
  let (Z :. _ :. _ :. m :. n) = extent arr
   in R.backpermute
        (extent arr)
        (\(Z :. i :. j :. k :. l) ->
           let mHalf = div m 2
               nHalf = div n 2
               a =
                 if k < mHalf
                   then k + mHalf
                   else k - mHalf
               b =
                 if l < nHalf
                   then l + nHalf
                   else l - nHalf
            in (Z :. i :. j :. a :. b))
        arr

{-# INLINE timeReversalR2S1RP #-}
timeReversalR2S1RP ::
     (R.Source s (Complex Double))
  => R.Array s DIM4 (Complex Double)
  -> R.Array D DIM4 (Complex Double)
timeReversalR2S1RP arr =
  let (Z :. _ :. _ :. n :. _) = extent arr
      freq = div (n - 1) 2
   in R.traverse arr id $ \f idx@(Z :. _ :. _ :. k :. _) ->
        let m = fromIntegral $ k - freq
         in f idx * (exp (0 :+ m * pi))


{-# INLINE generateDFTPlanR2S1RP #-}
generateDFTPlanR2S1RP :: DFTPlan -> R.Array U DIM4 (Complex Double) -> IO DFTPlan
generateDFTPlanR2S1RP plan arr = do
  let (Z :. rows :. cols :. orientations :. scales) = extent arr
  lock <- getFFTWLock
  vecTemp1 <-
    VS.fromList <$> M.replicateM (orientations * rows * cols * scales) randomIO :: IO (VS.Vector Double)
  vecTemp2 <-
    VS.fromList <$> M.replicateM (orientations * rows * cols * scales) randomIO :: IO (VS.Vector Double)
  (plan1, vecTemp3) <-
    dft1dGPlan
      lock
      plan
      ([rows, cols, orientations, scales])
      [0, 1]
      (VS.zipWith mkPolar vecTemp1 vecTemp2)
  (plan2, _) <-
    idft1dGPlan lock plan1 ([rows, cols, orientations, scales]) [0, 1] vecTemp3
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
    VS.fromList <$> M.replicateM (orientations * rows * cols * scales) randomIO :: IO (VS.Vector Double)
  vecTemp7 <-
    VS.fromList <$> M.replicateM (orientations * rows * cols * scales) randomIO :: IO (VS.Vector Double)
  (plan4, vecTemp8) <-
    dft1dGPlan
      lock
      plan3
      ([rows, cols, orientations, scales])
      [2, 3]
      (VS.zipWith mkPolar vecTemp6 vecTemp7)
  (plan5, _) <-
    idft1dGPlan lock plan4 ([rows, cols, orientations, scales]) [2, 3] vecTemp8
  return plan5


{-# INLINE crosscorrelationR2S1RP #-}
crosscorrelationR2S1RP ::
     DFTPlan
  -> R.Array U DIM2 (Complex Double)
  -> R.Array U DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
crosscorrelationR2S1RP plan img harmonics = do
  let (Z :. rows :. cols :. m :. n) = extent harmonics
      planID2D = DFTPlanID DFT1DG ([rows, cols]) [0, 1]
      planID = DFTPlanID DFT1DG ([rows, cols, m, n]) [0, 1]
      inversePlanID = DFTPlanID IDFT1DG ([rows, cols, m, n]) [0, 1]
      imgVec = VU.convert . toUnboxed $ img
      harmonicsVec =
        VU.convert . toUnboxed . computeS . makeFilterR2S1RP $ harmonics
  imgVecF <- dftExecute plan planID2D imgVec
  harmonicsVecF <- dftExecute plan planID harmonicsVec
  let convolvedVecF =
        computeS .
        R.traverse2
          (fromUnboxed (extent harmonics) . VS.convert . VS.map conjugate $
           harmonicsVecF)
          (fromUnboxed (extent img) . VS.convert $ imgVecF)
          const $
        (\f4d f2d idx@(Z :. i :. j :. _ :. _) -> f2d (Z :. i :. j) * f4d idx)
  convolvedVec <-
    dftExecute plan inversePlanID . VU.convert . toUnboxed $ convolvedVecF
  return . fromUnboxed (extent harmonics) . VS.convert $ convolvedVec


{-# INLINE completionR2S1RP #-}
completionR2S1RP ::
     (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  => DFTPlan
  -> R.Array s1 DIM4 (Complex Double)
  -> R.Array s2 DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
completionR2S1RP plan arr1 arr2 = do
  let (Z :. rows :. cols :. oris :. scales) = extent arr1
      planID = DFTPlanID DFT1DG [rows, cols, oris, scales] [2, 3]
      inversePlanID = DFTPlanID IDFT1DG [rows, cols, oris, scales] [2, 3]
      vec1 = VU.convert . toUnboxed . computeS . delay $ arr1
      vec2 = VU.convert . toUnboxed . computeS . makeFilterR2S1RP $ arr2
  vec1F <- dftExecute plan planID vec1
  vec2F <- dftExecute plan planID vec2
  convolvedVec <- dftExecute plan inversePlanID $ VS.zipWith (*) vec1F vec2F
  return . fromUnboxed (extent arr1) . VS.convert $ convolvedVec
