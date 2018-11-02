import           App.PinwheelBasis.Convolution
import           App.PinwheelBasis.Pinwheel
import           Control.Monad                 as M
import           Control.Parallel.Strategies
import           Data.Array                    as Arr
import           Data.Array.Repa               as R
import           Data.Complex
import           Data.List                     as L
import           Data.Vector.Storable          as VS
import           Data.Vector.Unboxed           as VU
import           Src.FokkerPlanck.MonteCarlo
import           Src.Image.ImageIO
import           Src.Utils.Coordinates
import           Src.Utils.DFT
import           System.Directory
import           System.Environment
import           System.FilePath
import           System.Random

{-# INLINE combineArgs #-}

combineArgs :: [String] -> [String]
combineArgs (x:y:[]) = [x L.++ y]
combineArgs (x:y:xs) = (x L.++ y) : (combineArgs xs)

{-# INLINE shift #-}
shift :: (Int,Int) ->  (Int,Int) -> (Int,Int)
shift (a,b) (c,d) = (a + c , b + d)

main = do
  args <- getArgs
  print args
  let (orientations:size:len:trails:threads:contrastN:_) =
        L.map (\x -> read x :: Int) . L.take 6 $ args
      (sigma:angularFreq:radialFreq:alpha:_) =
        L.map (\x -> read x :: Double) . L.drop 6 $ args
      (s:_) = L.map (\x -> read x :: (Int,Int)) . L.drop 10 $args
      xs = L.map (\x -> read x :: (Int, Int)) . L.drop 11 $ args
      img =
        fromListUnboxed (Z :. size :. size) .
        L.map (:+ 0) .
        Arr.elems .
        accumArray (+) (0 :: Double) ((0, 0), (size - 1, size - 1)) .
        L.map (\idx -> (idx, 255)) $
        xs
      img90 =
        fromListUnboxed (Z :. size :. size) .
        L.map (:+ 0) .
        Arr.elems .
        accumArray (+) (0 :: Double) ((0, 0), (size - 1, size - 1)) .
        L.map (\idx -> (shift s idx, 255)) $
        [(64, 64), (64, 65), (64, 66)]
      img45 =
        fromListUnboxed (Z :. size :. size) .
        L.map (:+ 0) .
        Arr.elems .
        accumArray (+) (0 :: Double) ((0, 0), (size - 1, size - 1)) .
        L.map (\idx -> (shift s idx, 255))  $
        [(64, 64), (63, 65), (62, 66)]
      harmonics@(HarmoicArray harmonicsArr) =
        generateR2S1HarmonicArray
          size
          size
          [-angularFreq .. angularFreq]
          [-radialFreq .. radialFreq]
          alpha
          Harmonic
      harmonicsC =
        generateR2S1HarmonicArray
          size
          size
          [-angularFreq .. angularFreq]
          [-radialFreq .. radialFreq]
          alpha
          ConjugateHarmonic
      folderName = "PinwheelShiftTwistWeightSharing"
  createDirectoryIfMissing False folderName
  arrG <-
    (computeS . R.map (:+ 0) . sumS) <$>
    solveMonteCarloR2S1
      threads
      trails
      size
      orientations
      sigma
      len
      (0, 0, 0 / 360 * 2 * pi, 1)
  plan <- generateDFTPlan getEmptyPlan harmonicsArr
  output <- convolve plan arrG img harmonics harmonicsC
  plotImageRepa (folderName </> "green.png") .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . VU.map magnitude . toUnboxed $
    arrG
  plotImageRepa (folderName </> "image.png") .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . VU.map magnitude . toUnboxed $
    img
  plotImageRepa (folderName </> "completionField.png") .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . VU.map magnitude . toUnboxed $
    output
  let freqFilter = projectFilter arrG harmonics
      freqImage45 = projectFilter img45 harmonics
      freqImage90 = projectFilter img90 harmonics
      convolved45 =
        recoverFilter
          (computeS $ R.zipWith (\a b -> conjugate a * b) freqFilter freqImage45)
          harmonicsC
      convolved90 =
        recoverFilter
          (computeS $ R.zipWith (\a b -> conjugate a * b) freqFilter freqImage90)
          harmonicsC
  plotImageRepa (folderName </> "convolution45.png") .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . VU.map magnitude . toUnboxed $
    convolved45
  plotImageRepa (folderName </> "convolution90.png") .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . VU.map magnitude . toUnboxed $
    convolved90
  output45 <- convolve plan arrG img45  harmonics harmonicsC
  output90 <- convolve plan arrG img90   harmonics harmonicsC
  plotImageRepa (folderName </> "completionField45.png") .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . VU.map magnitude . toUnboxed $
    output45
  plotImageRepa (folderName </> "completionField90.png") .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . VU.map magnitude . toUnboxed $
    output90
