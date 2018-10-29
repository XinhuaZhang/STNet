import           App.STHarmonics.GeneralizedConvolution
import           App.STHarmonics.Utils
import           Control.Monad                          as M
import           Data.Array                             as Arr
import           Data.Array.Repa                        as R
import           Data.Complex
import           Data.List                              as L
import           Data.Vector.Storable                   as VS
import           Data.Vector.Unboxed                    as VU
import           Src.Array.TwoHalfDArray
import           Src.FokkerPlanck.MonteCarlo
import           Src.Image.ImageIO
import           Src.Utils.Coordinates
import           Src.Utils.DFT
import           System.Directory
import           System.Environment
import           System.FilePath
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
          ((0, 0, 0), (size - 1, size - 1, orientations - 1)) .
        L.map (\idx -> (idx, 1 / (fromIntegral . L.length $ xs))) .
        L.map
          (\(Source (t, x, y)) ->
             let theta1 =
                   angleFunctionDeg
                     (fromIntegral (x - div size 2))
                     (fromIntegral (y - div size 2))
                 theta2 = checkDeg t
                 r =
                   sqrt $
                   (fromIntegral (x - div size 2)) ** 2 +
                   (fromIntegral (y - div size 2)) ** 2
                 theta3 = (theta1 - theta2) / 2
                 x' = round $ r * (cos . deg2Rad $ theta3)
                 y' = round $ r * (sin . deg2Rad $ theta3)
             in ( x' + div size 2
                , y' + div size 2
                , floor $
                  (checkDeg ((theta1 + theta2) / 2)) / 360 *
                  fromIntegral orientations)) $
        xs
      arrSourceRepa =
        fromListUnboxed (Z :. size :. size :. orientations) . Arr.elems $
        arrSource
      arrSink =
        accumArray
          (+)
          (0 :: Double)
          ((0, 0, 0), (size - 1, size - 1, orientations - 1)) .
        L.map (\idx -> (idx, 1 / (fromIntegral . L.length $ xs))) .
        L.map
          (\(Sink (t, x, y)) ->
             let theta1 =
                   angleFunctionDeg
                     (fromIntegral (x - div size 2))
                     (fromIntegral (y - div size 2))
                 theta2 = checkDeg t
                 r =
                   sqrt $
                   (fromIntegral (x - div size 2)) ** 2 +
                   (fromIntegral (y - div size 2)) ** 2
                 theta3 = (theta1 - theta2) / 2
                 x' = round $ r * (cos . deg2Rad $ theta3)
                 y' = round $ r * (sin . deg2Rad $ theta3)
             in ( x' + div size 2
                , y' + div size 2
                , floor $
                  (checkDeg ((theta1 + theta2) / 2)) / 360 *
                  fromIntegral orientations)) $
        ys
      arrSinkRepa =
        fromListUnboxed (Z :. size :. size :. orientations) . Arr.elems $
        arrSink
      folderName = "GreensData"
  print xs
  print ys
  arrG <-
    readRepaArray
      (folderName </>
       ("GeneralizedGreens_" L.++ show size L.++ "_" L.++ show orientations L.++
        "_" L.++
        show sigma L.++
        ".dat"))
  -- Generate DFT plan
  lock <- getFFTWLock
  vecTemp1 <-
    VS.fromList <$> M.replicateM (orientations * size * size) randomIO :: IO (VS.Vector Double)
  vecTemp2 <-
    VS.fromList <$> M.replicateM (orientations * size * size) randomIO :: IO (VS.Vector Double)
  (dftPlan1, vecTemp4) <-
    dft1dGPlan
      lock
      getEmptyPlan
      [size, size, orientations]
      [0, 1, 2]
      (VS.zipWith mkPolar vecTemp1 vecTemp2)
  (dftPlan, _) <-
    idft1dGPlan lock dftPlan1 [size, size, orientations] [0, 1, 2] vecTemp4
  createDirectoryIfMissing True "GeneralizedConvolutionTest"
  plotImageRepa "GeneralizedConvolutionTest/green.png" .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) .
    R.toUnboxed . R.sumS . convert' $
    arrG
  -- Source Field
  outputSource <-
    convolve
      dftPlan
      CrossCorrelationST
      (computeS . makeFilter $ arrG)
      arrSourceRepa
  plotImageRepa "GeneralizedConvolutionTest/source.png" .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) .
    R.toUnboxed . R.sumS . convert' $
    outputSource
  -- Sink Field
  outputSink <-
    convolve
      dftPlan
      CrossCorrelationST
      (computeS . makeFilter $ arrG)
      arrSinkRepa
  plotImageRepa "GeneralizedConvolutionTest/sink.png" .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) .
    R.toUnboxed . R.sumS . convert' $
    outputSink
  -- Completion field
  plotImageRepa "GeneralizedConvolutionTest/completionField.png" .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . R.toUnboxed . R.sumS $
    R.zipWith (*) (convert' outputSource) (timeReversal . convert' $ outputSink)
