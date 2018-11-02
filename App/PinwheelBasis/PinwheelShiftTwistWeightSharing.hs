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

main = do
  args <- getArgs
  print args
  let (orientations:size:len:trails:threads:contrastN:_) =
        L.map (\x -> read x :: Int) . L.take 6 $ args
      (sigma:angularFreq:radialFreq:alpha:_) =
        L.map (\x -> read x :: Double) . L.drop 6 $ args
      xs = L.map (\x -> read x :: (Int, Int)) . combineArgs . L.drop 10 $ args
      img =
        fromListUnboxed (Z :. size :. size) .
        L.map (:+ 0) .
        Arr.elems .
        accumArray (+) (0 :: Double) ((0, 0), (size - 1, size - 1)) .
        L.map (\idx -> (idx, 255)) $
        xs
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
