module App.STHarmonics.GeneralizedConvolution where

import           App.STHarmonics.R2S1RPHarmonics
import           Control.Monad                   as M
import           Data.Array                      as Arr
import           Data.Array.Repa                 as R
import           Data.Complex
import           Data.List                       as L
import           Data.Vector.Storable            as VS
import           Data.Vector.Unboxed             as VU
import           Src.Array.CoordinatesTransform
import           Src.Array.TwoHalfDArray
import           Src.Utils.Coordinates
import           Src.Utils.DFT
import           System.Random
import Control.Parallel.Strategies

data ConvolutionTypeST
  = ConvolutionST
  | CrossCorrelationST
  deriving (Read, Show)
  
-- For filters

{-# INLINE project #-}

project
  :: (Shape sh, R.Source r e, Num e)
  => R.Array r sh e -> [R.Array r sh e] -> [e]
project arr = parMap rseq (R.sumAllS . R.zipWith (*) arr)

{-# INLINE recover #-}

recover
  :: (Shape sh, R.Source r e, Num e, Unbox e)
  => [e] -> [R.Array r sh e] -> R.Array U sh e
recover xs =
  computeS .
  L.foldl1' (R.zipWith (+)) . L.zipWith (\x arr -> R.map (* x) arr) xs
  
{-# INLINE recoverP #-}

recoverP
  :: (Shape sh, R.Source r e, Num e, Unbox e)
  => [e] -> [R.Array r sh e] -> IO (R.Array U sh e)
recoverP xs =
  computeP .
  L.foldl1' (R.zipWith (+)) . L.zipWith (\x arr -> R.map (* x) arr) xs
  
{-# INLINE projectFilterR2S1 #-}

projectFilterR2S1 :: R.Array U DIM3 (Complex Double)
                  -> [R2S1RPHarmonic]
                  -> [Complex Double]
projectFilterR2S1 arr =
  project arr -- . L.map (computeS . makeFilterR2S1)
  . L.map (\(R2S1Harmonic (HarmoicArray harmonic)) -> harmonic)


{-# INLINE recoverFilterR2S1 #-}

recoverFilterR2S1 :: [Complex Double]
                  -> [R2S1RPHarmonic]
                  -> R.Array U DIM3 (Complex Double)
recoverFilterR2S1 xs =
  recover xs . L.map (\(R2S1Harmonic (ConjugateHarmoicArray arr)) -> arr)
  
{-# INLINE recoverFilterR2S1P #-}

recoverFilterR2S1P :: [Complex Double]
                   -> [R2S1RPHarmonic]
                   -> IO (R.Array U DIM3 (Complex Double))
recoverFilterR2S1P xs =
  recoverP xs . L.map (\(R2S1Harmonic (ConjugateHarmoicArray arr)) -> arr)

-- For inputs
{-# INLINE projectR2S1 #-}

projectR2S1
  :: DFTPlan
  -> R.Array U DIM3 (Complex Double)
  -> [R2S1RPHarmonic]
  -> IO [R.Array U DIM2 (Complex Double)]
projectR2S1 plan input harmonics = do
  let (Z :. rows :. cols :. oris) = extent input
      planID = DFTPlanID DFT1DG ([rows, cols, oris]) [0, 1]
      inversePlanID = DFTPlanID IDFT1DG ([rows, cols, oris]) [0, 1]
      inputVec = VU.convert . toUnboxed $ input
      harmonicVecList =
        L.map
          (\(R2S1Harmonic (HarmoicArray harmonic)) ->
             VU.convert . toUnboxed . computeS . makeFilterR2S1 $ harmonic) $
        harmonics
  inputVecF <- dftExecute plan planID inputVec
  harmonicVecListF <- dftExecuteBatchP plan planID harmonicVecList
  convolvedVecList <-
    dftExecuteBatchP plan inversePlanID .
    parMap rdeepseq (VS.zipWith (*) inputVecF) $
    harmonicVecListF
  return . L.map (R.sumS . fromUnboxed (extent input) . VS.convert) $
    convolvedVecList

{-# INLINE recoverR2S1 #-}

recoverR2S1
  :: DFTPlan
  -> [R.Array U DIM2 (Complex Double)]
  -> [R2S1RPHarmonic]
  -> IO (R.Array U DIM3 (Complex Double))
recoverR2S1 plan arrs harmonics = do
  let (Z :. rows :. cols :. oris) =
        (\(R2S1Harmonic (ConjugateHarmoicArray harmonic)) -> extent harmonic) .
        L.head $
        harmonics
      planID2D = DFTPlanID DFT1DG [rows, cols] [0, 1]
      planID = DFTPlanID DFT1DG [rows, cols, oris] [0, 1]
      inversePlanID = DFTPlanID IDFT1DG [rows, cols, oris] [0, 1]
      inputVecList = L.map (VU.convert . toUnboxed) arrs
      harmonicVecList =
        L.map
          (\(R2S1Harmonic (ConjugateHarmoicArray harmonic)) ->
             VU.convert . toUnboxed . computeS . makeFilterR2S1 $ harmonic)
          harmonics
  inputVecListF <- dftExecuteBatchP plan planID2D inputVecList
  harmonicVecListF <- dftExecuteBatchP plan planID harmonicVecList
  let convolvedVecListF =
        withStrategy (parList rdeepseq) .
        L.zipWith
          (\inputVecF harmonicVecF ->
             let inputArrF =
                   fromUnboxed (Z :. rows :. cols) . VS.convert $ inputVecF
                 harmonicsArrF =
                   fromUnboxed (Z :. rows :. cols :. oris) .
                   VS.convert . VS.map conjugate $
                   harmonicVecF
             in VU.convert .
                toUnboxed .
                computeS . R.traverse2 inputArrF harmonicsArrF (\_ sh -> sh) $
                (\f2d f3d idx@(Z :. i :. j :. k) -> f2d (Z :. i :. j) * f3d idx))
          inputVecListF $
        harmonicVecListF
  convolvedVecList <- dftExecuteBatchP plan inversePlanID convolvedVecListF
  return .
    fromUnboxed (Z :. rows :. cols :. oris) .
    VS.convert . L.foldl1' (VS.zipWith (+)) $
    convolvedVecList
    

{-# INLINE makeFilterR2S1 #-}

makeFilterR2S1 :: (R.Source s e) => R.Array s DIM3 e -> R.Array D DIM3 e 
makeFilterR2S1 arr =
  let (Z :. rows :. cols :. _ ) = extent arr
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
       
{-# INLINE makeFilterR2S1' #-}

makeFilterR2S1' :: (R.Source s e) => R.Array s DIM2 e -> R.Array D DIM2 e 
makeFilterR2S1' arr =
  let (Z :. rows :. cols) = extent arr
  in R.backpermute
       (extent arr)
       (\(Z :. i :. j) ->
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
          in (Z :. x :. y))
       arr
       
{-# INLINE convolveR2S1 #-}

convolveR2S1
  :: DFTPlan
  -> R.Array U DIM3 (Complex Double)
  -> R.Array U DIM3 (Complex Double)
  -> [R2S1RPHarmonic]
  -> [R2S1RPHarmonic]
  -> [R2S1RPHarmonic]
  -> IO (R.Array U DIM3 (Complex Double))
convolveR2S1 plan filter input harmonics harmonics1 conjugateHarmonics = do
  let coefficients = projectFilterR2S1 filter harmonics1
  coefficientsR2 <- projectR2S1 plan input harmonics
  let convolvedR2 =
        L.zipWith
          (\c arr -> computeS . R.map (* c) $ arr)
          coefficients
          coefficientsR2
  recoverR2S1 plan convolvedR2 conjugateHarmonics

{-# INLINE generateR2S1DFTPlan #-}

generateR2S1DFTPlan :: DFTPlan -> R.Array U DIM3 (Complex Double) -> IO DFTPlan
generateR2S1DFTPlan plan arr = do
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
      ([rows,cols,orientations])
      [0, 1]
      (VS.zipWith mkPolar vecTemp1 vecTemp2)
  (plan2, _) <-
    idft1dGPlan
      lock
      plan1
      ([rows,cols,orientations])
      [0, 1]
      vecTemp3
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

{-# INLINE timeReversalR2S1 #-}

timeReversalR2S1 :: R.Array U DIM3 Double -> R.Array U DIM3 Double
timeReversalR2S1 arr =
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
       

{-# INLINE projectFilterR2S1RP #-}

projectFilterR2S1RP :: R.Array U DIM4 (Complex Double)
                    -> [R2S1RPHarmonic]
                    -> [Complex Double]
projectFilterR2S1RP arr =
  project arr . L.map (\(R2S1RPHarmonic (HarmoicArray harmonic)) -> harmonic)

{-# INLINE recoverFilterR2S1RP #-}

recoverFilterR2S1RP :: [Complex Double]
                    -> [R2S1RPHarmonic]
                    -> R.Array U DIM4 (Complex Double)
recoverFilterR2S1RP xs =
  recover xs . L.map (\(R2S1RPHarmonic (ConjugateHarmoicArray arr)) -> arr)

-- For inputs
{-# INLINE projectR2S1RP #-}

projectR2S1RP
  :: DFTPlan
  -> R.Array U DIM4 (Complex Double)
  -> [R2S1RPHarmonic]
  -> IO [R.Array U DIM2 (Complex Double)]
projectR2S1RP plan arr harmonics = do
  let planID = DFTPlanID DFT1DG (L.reverse . listOfShape . extent $ arr) [0, 1]
      inversePlanID =
        DFTPlanID IDFT1DG (L.reverse . listOfShape . extent $ arr) [0, 1]
      inputVec = VU.convert . toUnboxed $ arr
      harmonicVecList =
        L.map
          (\(R2S1RPHarmonic (HarmoicArray harmonic)) ->
             VU.convert . toUnboxed $ harmonic) $
        harmonics
  inputVecF <- dftExecute plan planID inputVec
  harmonicVecListF <- dftExecuteBatchP plan planID harmonicVecList
  convolvedVecList <-
    dftExecuteBatchP plan inversePlanID .
    parMap rdeepseq (VS.zipWith (*) inputVecF) $
    harmonicVecListF
  return . L.map (R.sumS . R.sumS . fromUnboxed (extent arr) . VS.convert) $
    convolvedVecList

{-# INLINE recoverR2S1RP #-}

recoverR2S1RP
  :: DFTPlan
  -> [R.Array U DIM2 (Complex Double)]
  -> [R2S1RPHarmonic]
  -> IO (R.Array U DIM4 (Complex Double))
recoverR2S1RP plan arrs harmonics = do
  let (Z :. rows :. cols :. oris :. scales) =
        (\(R2S1RPHarmonic (ConjugateHarmoicArray harmonic)) -> extent harmonic) .
        L.head $
        harmonics
      planID2D = DFTPlanID DFT1DG [rows, cols] [0, 1]
      planID = DFTPlanID DFT1DG [rows, cols, oris, scales] [0, 1]
      inversePlanID = DFTPlanID IDFT1DG [rows, cols, oris, scales] [0, 1]
      inputVecList = L.map (VU.convert . toUnboxed) arrs
      harmonicVecList =
        L.map
          (\(R2S1RPHarmonic (ConjugateHarmoicArray harmonic)) ->
             VU.convert . toUnboxed $ harmonic)
          harmonics
  inputVecListF <- dftExecuteBatchP plan planID2D inputVecList
  harmonicVecListF <- dftExecuteBatchP plan planID harmonicVecList
  let convolvedVecListF =
        withStrategy (parList rdeepseq) .
        L.zipWith
          (\inputVecF harmonicVecF ->
             let inputArrF =
                   fromUnboxed (Z :. rows :. cols) . VS.convert $ inputVecF
                 harmonicsArrF =
                   fromUnboxed (Z :. rows :. cols :. oris :. scales) .
                   VS.convert . VS.map conjugate $
                   harmonicVecF
             in VU.convert .
                toUnboxed .
                computeS . R.traverse2 inputArrF harmonicsArrF (\_ sh -> sh) $
                (\f2d f4d idx@(Z :. i :. j :. k :. l) ->
                   f2d (Z :. i :. j) * f4d idx))
          inputVecListF $
        harmonicVecListF
  convolvedVecList <- dftExecuteBatchP plan inversePlanID convolvedVecListF
  return .
    fromUnboxed (Z :. rows :. cols :. oris :. scales) .
    VS.convert . L.foldl1' (VS.zipWith (+)) $
    convolvedVecList


{-# INLINE makeFilterR2S1RP #-}

makeFilterR2S1RP :: (R.Source s e) => R.Array s DIM4 e -> R.Array D DIM4 e 
makeFilterR2S1RP arr =
  let (Z :. rows :. cols :. _ :. _ ) = extent arr
  in R.backpermute
       (extent arr)
       (\(Z :. i :. j :. k :. l) ->
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
          in (Z :. x :. y :. k :. l))
       arr

{-# INLINE convolveR2S1RP #-}

convolveR2S1RP
  :: DFTPlan
  -> R.Array U DIM4 (Complex Double)
  -> R.Array U DIM4 (Complex Double)
  -> [R2S1RPHarmonic]
  -> [R2S1RPHarmonic]
  -> IO (R.Array U DIM4 (Complex Double))
convolveR2S1RP plan filter input harmonics conjugateHarmonics = do
  let coefficients = projectFilterR2S1RP filter harmonics
  coefficientsR2 <- projectR2S1RP plan input harmonics
  let convolvedR2 =
        L.zipWith
          (\c arr -> computeS . R.map (* c) $ arr)
          coefficients
          coefficientsR2
  recoverR2S1RP plan convolvedR2 conjugateHarmonics

{-# INLINE generateR2S1RPDFTPlan #-}

generateR2S1RPDFTPlan :: DFTPlan -> R.Array U DIM4 (Complex Double) -> IO DFTPlan
generateR2S1RPDFTPlan plan arr = do
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
      (L.reverse . listOfShape . extent $ arr)
      [0, 1]
      (VS.zipWith mkPolar vecTemp1 vecTemp2)
  (plan2, _) <-
    idft1dGPlan
      lock
      plan1
      (L.reverse . listOfShape . extent $ arr)
      [0, 1]
      vecTemp3
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

{-# INLINE timeReversalR2S1RP #-}

timeReversalR2S1RP :: R.Array U DIM4 Double -> R.Array U DIM4 Double
timeReversalR2S1RP arr =
  let (Z :. _ :. _ :. nf :. _) = extent arr
      n = div nf 2
  in computeS $
     R.backpermute
       (extent arr)
       (\(Z :. i :. j :. k :. l) ->
          let x = nf - n
          in if k >= x
               then (Z :. i :. j :. k - x :. l)
               else (Z :. i :. j :. k + n :. l))
       arr
