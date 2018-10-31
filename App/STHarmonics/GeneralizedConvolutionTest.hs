import           App.STHarmonics.GeneralizedConvolution
import           App.STHarmonics.R2S1RPHarmonics
import           App.STHarmonics.Utils
import           Control.Monad                          as M
import           Control.Parallel.Strategies
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
  print args
  let (orientations:size:len:trails:threads:contrastN:_) =
        L.map (\x -> read x :: Int) . L.take 6 $ args
      (sigma:angularFreq:radialFreq:alpha:_) =
        L.map (\x -> read x :: Double) . L.drop 6 $ args
      -- Generate initial distribution: arrRepa
      (xs:ys:[]) =
        L.group .
        L.sort . L.map (\x -> read x :: STPoint) . combineArgs . L.drop 10 $
        args
      arrSource =
        accumArray
          (+)
          (0 :: Double)
          ((0, 0, 0), (size - 1, size - 1, orientations - 1)) .
        L.map (\idx -> (idx, 1 / (fromIntegral . L.length $ xs))) .
        L.map
          (\(Source (t, x, y)) ->
             (x, y, floor ((checkDeg t) / 360 * fromIntegral orientations))) $
        xs
      arrSourceRepa =
        fromListUnboxed (Z :. size :. size :. orientations) .
        L.map (:+ 0) . Arr.elems $
        arrSource
      arrSink =
        accumArray
          (+)
          (0 :: Double)
          ((0, 0, 0), (size - 1, size - 1, orientations - 1)) .
        L.map (\idx -> (idx, 1 / (fromIntegral . L.length $ xs))) .
        L.map
          (\(Sink (t, x, y)) ->
             (x, y, floor ((checkDeg t) / 360 * fromIntegral orientations))) $
        ys
      arrSinkRepa =
        fromListUnboxed (Z :. size :. size :. orientations) .
        L.map (:+ 0) . Arr.elems $
        arrSink
      r2s1Harmonics =
        parMap
          rseq
          (\(af1, af2, rf) ->
             generateR2S1HarmonicArray
               orientations
               size
               size
               af1
               af2
               rf
               alpha
               Harmonic)
          [ (af1, af2, rf)
          | af1 <- [-angularFreq .. angularFreq]
          , af2 <- [-angularFreq .. angularFreq]
          , rf <- [0 .. radialFreq]
          ]
      r2s1HarmonicsC =
        parMap
          rseq
          (\(af1, af2, rf) ->
             generateR2S1HarmonicArray
               orientations
               size
               size
               af1
               af2
               rf
               alpha
               ConjugateHarmonic)
          [ (af1, af2, rf)
          | af1 <- [-angularFreq .. angularFreq]
          , af2 <- [-angularFreq .. angularFreq]
          , rf <- [0 .. radialFreq]
          ]
      folderName = "GreensData"
  print xs
  print ys
  arrG <-
    (computeS . R.map (:+ 0)) <$>
    readRepaArray
      (folderName </>
       ("Greens_" L.++ show size L.++ "_" L.++ show orientations L.++ "_" L.++
        show sigma L.++
        ".dat")) :: IO (R.Array U DIM3 (Complex Double))
  -- Generate DFT plan
  dftPlan <- generateR2S1DFTPlan getEmptyPlan arrSinkRepa
  createDirectoryIfMissing True "GeneralizedConvolutionTest"
  plotImageRepa "GeneralizedConvolutionTest/green.png" .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) .
    R.toUnboxed . R.sumS . R.map magnitude $
    arrG
  -- Source Field
  outputSource <-
    (computeS . R.map magnitude) <$>
    convolveR2S1
      dftPlan
      (computeS . makeFilterR2S1 $ arrG)
      arrSourceRepa
      r2s1Harmonics
      r2s1HarmonicsC :: IO (R.Array U DIM3 Double)
  plotImageRepa "GeneralizedConvolutionTest/source.png" .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . R.toUnboxed . R.sumS $
    outputSource
  -- Sink Field
  outputSink <-
    (timeReversalR2S1 . computeS . R.map magnitude) <$>
    convolveR2S1
      dftPlan
      (computeS . makeFilterR2S1 $ arrG)
      arrSinkRepa
      r2s1Harmonics
      r2s1HarmonicsC
  plotImageRepa "GeneralizedConvolutionTest/sink.png" .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . R.toUnboxed . R.sumS $
    outputSink
  -- Completion field
  plotImageRepa "GeneralizedConvolutionTest/completionField.png" .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . R.toUnboxed . R.sumS $
    R.zipWith (*) outputSource outputSink
