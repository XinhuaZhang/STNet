module App.STHarmonics.GeneralizedConvolution where

import           App.STHarmonics.Generator
import           Src.Utils.Coordinates
import           Control.Monad                  as M
import           Data.Array.Repa                as R
import Data.Array as Arr
import           Data.Complex
import           Data.List                      as L
import           Data.Vector.Storable           as VS
import           Data.Vector.Unboxed            as VU
import           Src.Array.CoordinatesTransform
import           Src.Array.TwoHalfDArray
import           Src.Utils.DFT
import           System.Random

data ConvolutionTypeST
  = ConvolutionST
  | CrossCorrelationST
  deriving (Read, Show)


{-# INLINE convolve25D #-}

convolve25D :: DFTPlan -> R2S1Array -> R2S1Array -> IO R2S1Array
convolve25D plan arr1 arr2 = do
  let planID = DFTPlanID DFT1DG (L.reverse . listOfShape . extent $ arr1) [1, 2]
      inversePlanID =
        DFTPlanID IDFT1DG (L.reverse . listOfShape . extent $ arr1) [1, 2]
      arrList1 = L.map VU.convert . get25DArrayFeatureVector $ arr1
      arrList2 = L.map VU.convert . get25DArrayFeatureVector $ arr2
  arrList1' <- dftExecuteBatch plan planID arrList1
  arrList2' <- dftExecuteBatch plan planID arrList2
  let arrList = L.zipWith (VS.zipWith (*)) arrList1' arrList2'
  arrList' <- dftExecuteBatch plan inversePlanID arrList
  return . fromUnboxed (extent arr1) . VS.convert . VS.concat $ arrList'
  
{-# INLINE projectR2S1 #-}

projectR2S1 :: R2S1Array -> [R2S1Array] -> [R.Array U DIM1 (Complex Double)]
projectR2S1 arr = L.map (sumS . sumS . R.zipWith (*) arr)

{-# INLINE projectR2S1' #-}

projectR2S1' :: DFTPlan -> R2S1Array -> [R2S1Array] -> IO [R2S1Array]
projectR2S1' plan arr = M.mapM (convolve25D plan arr)

{-# INLINE recoverR2S1 #-}

recoverR2S1 :: DFTPlan -> [R2S1Array] -> [R2S1Array] -> IO R2S1Array
recoverR2S1 plan arrList1 arrList2 =
  (computeS . L.foldl1' (R.zipWith (+)) . L.map delay) <$>
  M.zipWithM (convolve25D plan) arrList1 arrList2


{-# INLINE generateDFTPlan #-}

generateDFTPlan :: DFTPlan -> R2S1Array -> IO DFTPlan
generateDFTPlan plan arr = do
  let (Z :. orientations :. rows :. cols) = extent arr
  lock <- getFFTWLock
  vecTemp1 <-
    VS.fromList <$> M.replicateM (orientations * rows * cols) randomIO :: IO (VS.Vector Double)
  vecTemp2 <-
    VS.fromList <$> M.replicateM (orientations * rows * cols) randomIO :: IO (VS.Vector Double)
  (plan1, vecTemp3) <-
    dft1dGPlan
      lock
      plan
      (L.reverse . listOfShape . extent $ arr)
      [1, 2]
      (VS.zipWith mkPolar vecTemp1 vecTemp2)
  (plan2, _) <-
    idft1dGPlan
      lock
      plan1
      (L.reverse . listOfShape . extent $ arr)
      [1, 2]
      (VS.zipWith mkPolar vecTemp1 vecTemp2)
  return plan2



{-# INLINE makeFilter #-}

makeFilter :: (R.Source s e) => R.Array s DIM3 e -> R.Array D DIM3 e 
makeFilter arr =
  let (Z :. rows :. cols :. orientations) = extent arr
  in R.backpermute
       (extent arr)
       (\(Z :. i :. j :. k) ->
          let halfRows = div rows 2
              halfCols = div cols 2
              halfOris = div orientations 2
              x =
                if i < halfRows
                  then i + halfRows
                  else i - halfRows
              y =
                if j < halfCols
                  then j + halfCols
                  else j - halfCols
              z =
                if k < halfOris
                  then k + halfOris
                  else k - halfOris
          in (Z :. x :. y :. z))
       arr

{-# INLINE convolve #-}

convolve
  :: DFTPlan
  -> ConvolutionTypeST
  -> R.Array U DIM3 Double
  -> R.Array U DIM3 Double
  -> IO (R.Array U DIM3 Double)
convolve dftPlan convolutionType arr1 arr2 = do
  let (Z :. rows :. cols :. orientations) = extent arr1
      planID = DFTPlanID DFT1DG [rows, cols, orientations] [0, 1, 2]
      inversePlanID = DFTPlanID IDFT1DG [rows, cols, orientations] [0, 1, 2]
  vec1 <-
    dftExecute dftPlan planID . VU.convert . VU.map (:+ 0) . toUnboxed $ arr1
  vec2 <-
    dftExecute dftPlan planID . VU.convert . VU.map (:+ 0) . toUnboxed $ arr2
  let vec1' =
        case convolutionType of
          ConvolutionST -> vec1
          CrossCorrelationST -> VS.map conjugate vec1
  vec <- dftExecute dftPlan inversePlanID (VS.zipWith (*) vec1' vec2)
  return .
    fromUnboxed (Z :. rows :. cols :. orientations) .
    VU.map realPart . VS.convert $
    vec


{-# INLINE timeReversal #-}

timeReversal :: R.Array U DIM3 Double -> R.Array U DIM3 Double
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


convert' :: R.Array U DIM3 Double -> R.Array U DIM3 Double
convert' arr =
  let (Z :. rows :. cols :. oris) = extent arr
      centerR = div rows 2
      centerC = div cols 2
      deltaTheta2 = (2 * pi) / (fromIntegral oris) :: Double
  in fromListUnboxed (extent arr) .
     Arr.elems .
     accumArray (+) 0 ((0, 0, 0), (rows - 1, cols - 1, oris - 1)) .
     L.filter
       (\((x, y, theta), _) ->
          if (x < rows) &&
             (y < cols) &&
             (theta < oris) && (x >= 0) && (y >= 0) && (theta >= 0)
            then True
            else False) .
     R.toList . R.traverse arr id $
     (\f idx@(Z :. i :. j :. k) ->
        let x' = fromIntegral $ i - centerR
            y' = fromIntegral $ j - centerC
            r = sqrt $ x' ** 2 + y' ** 2 :: Double
            theta1 = angleFunctionRad x' y' :: Double
            theta2 = fromIntegral k * deltaTheta2
            theta1' = (theta1 + theta2) / 1
            theta2' = thetaCheck $ (theta1 - theta2) / 1
            x = (round $ r * cos theta1') + centerR
            y = (round $ r * sin theta1') + centerC
        in ((x, y, floor $ theta2' / deltaTheta2 * fromIntegral oris), f idx))

{-# INLINE thetaCheck #-}

thetaCheck :: Double -> Double
thetaCheck theta =
  if theta < 0
    then theta + 2 * pi
    else if theta >= 2 * pi
           then theta - 2 * pi
           else theta
