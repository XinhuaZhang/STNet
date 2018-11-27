import           App.STHarmonics.WeightSharing
import           Control.Monad                 as M
import           Data.Array                    as Arr
import           Data.Array.Repa               as R
import           Data.Complex
import           Data.List                     as L
import           Data.Vector.Storable          as VS
import           Data.Vector.Unboxed           as VU
import           Src.Array.Transform           (rotate3D)
import           Src.Array.TwoHalfDArray
import           Src.FokkerPlanck.MonteCarlo
import           Src.Image.ImageIO
import           System.Environment
import           System.Random

data STPoint
  = Source (Double, Int, Int)
  | Sink (Double, Int, Int)
  deriving (Read,Show)

instance Eq STPoint where
  (==) (Source _) (Source _) = True
  (==) (Sink _) (Sink _)     = True
  (==) (Source _) (Sink _)   = False
  (==) (Sink _) (Source _)   = False

instance Ord STPoint where
  compare (Source x) (Source y) = compare x y
  compare (Sink x) (Sink y)     = compare x y
  compare (Source _) (Sink _)   = LT
  compare (Sink _) (Source _)   = GT


{-# INLINE reduceContrast #-}

reduceContrast
  :: Int -> VU.Vector Double -> VU.Vector Double
reduceContrast idx vec =
  let x = L.head . L.drop idx . L.reverse . L.sort . VU.toList $ vec
  in VU.map
       (\y ->
          if y >= x
            then x
            else y)
       vec

{-# INLINE checkDeg #-}

checkDeg :: Double -> Double
checkDeg deg =
  if deg < 0
    then checkDeg (deg + 360)
    else if deg >= 360
           then checkDeg (deg - 360)
           else deg

{-# INLINE combineArgs #-}

combineArgs :: [String] -> [String]
combineArgs (x:y:[]) = [x L.++ " " L.++ y]
combineArgs (x:y:xs) = (x L.++ " " L.++ y) : combineArgs xs

main = do
  args <- getArgs
  let (orientations:size:len:trails:threads:contrastN:_) =
        L.map (\x -> read x :: Int) . L.take 6 $ args
      (sigma:_) = L.map (\x -> read x :: Double) . L.drop 6 $ args
      -- Generate initial distribution: arrRepa
      (xs:ys:[]) =
        L.group .
        L.sort . L.map (\x -> read x :: STPoint) . combineArgs . L.drop 7 $
        args
      arrSource =
        accumArray
          (+)
          (0 :: Double)
          ((0, 0, 0), (orientations - 1, size - 1, size - 1)) .
        L.map (\idx -> (idx, 1 / (fromIntegral . L.length $ xs))) .
        L.map
          (\(Source (t, x, y)) ->
             (floor ((checkDeg t) / 360 * fromIntegral orientations), x, y)) $
        xs
      arrSourceRepa =
        fromListUnboxed (Z :. orientations :. size :. size) . Arr.elems $
        arrSource
      arrSink =
        accumArray
          (+)
          (0 :: Double)
          ((0, 0, 0), (orientations - 1, size - 1, size - 1)) .
        L.map (\idx -> (idx, 1 / (fromIntegral . L.length $ xs))) .
        L.map
          (\(Sink (t, x, y)) ->
             (floor ((checkDeg t) / 360 * fromIntegral orientations), x, y)) $
        ys
      arrSinkRepa =
        fromListUnboxed (Z :. orientations :. size :. size) . Arr.elems $
        arrSink
  print xs
  print ys
  -- Generate DFT plan
  lock <- getFFTWLock
  vecTemp1 <-
    VS.fromList <$> M.replicateM (orientations * size * size) randomIO :: IO (VS.Vector Double)
  vecTemp2 <-
    VS.fromList <$> M.replicateM (orientations * size * size) randomIO :: IO (VS.Vector Double)
  vecTemp2D1 <-
    VS.fromList <$> M.replicateM (size * size) randomIO :: IO (VS.Vector Double)
  vecTemp2D2 <-
    VS.fromList <$> M.replicateM (size * size) randomIO :: IO (VS.Vector Double)
  (dftPlan1, vecTemp4) <-
    dft1dGPlan
      lock
      getEmptyPlan
      [orientations, size, size]
      [1, 2]
      (VS.zipWith mkPolar vecTemp1 vecTemp2)
  (dftPlan2, _) <-
    idft1dGPlan lock dftPlan1 [orientations, size, size] [1, 2] vecTemp4
  (dftPlan, _) <-
    dft1dGPlan
      lock
      dftPlan2
      [size, size]
      [0, 1]
      (VS.zipWith mkPolar vecTemp2D1 vecTemp2D2)
  -- Generate Green's function
  arrG <-
    solveMonteCarloR2S1
      threads
      trails
      size
      orientations
      sigma
      len
      (0, 0, 0 / 360 * pi, 0)
  let arrG'
        -- computeS .
        -- makeFilter .
       =
        fromUnboxed (Z :. (1 :: Int) :. size :. size) .
        reduceContrast contrastN . toUnboxed . R.sumS $
        arrG
  plotImageRepa "green.png" . Image 8 $ arrG'
  -- Source field
  outputSource <-
    shareWeightST dftPlan CrossCorrelationST (0, 1) arrSourceRepa .
    computeS . rotate3D $
    arrG
  plotImageRepa "source.png" .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) .
    reduceContrast contrastN . L.foldl1' (VU.zipWith (+)) . get25DArray2DVector $
    outputSource
  -- Sink field
  outputSink
    -- timeReversal <$>
     <-
    (shareWeightST dftPlan CrossCorrelationST (0, 1) arrSinkRepa .
     computeS . rotate3D $
     arrG)
  -- let outputSink = rotateST' (0,100) outputSink' (div orientations 2)
  plotImageRepa "sink.png" .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) .
    reduceContrast contrastN . L.foldl1' (VU.zipWith (+)) . get25DArray2DVector $
    outputSink
  -- Completion field
  plotImageRepa "completionField.png" .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) .
    -- reduceContrast contrastN .
    L.foldl1' (VU.zipWith (+)) . get25DArray2DVector $
    R.zipWith (*) outputSource (timeReversal outputSink)
  plotImageRepa "completionField1.png" .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) $
    VU.zipWith
      (*)
      (L.foldl1' (VU.zipWith (+)) . get25DArray2DVector $ outputSource)
      (L.foldl1' (VU.zipWith (+)) . get25DArray2DVector $ outputSink)
