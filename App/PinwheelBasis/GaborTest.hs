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

main = do
  args <- getArgs
  -- print args
  let (orientations:size:len:trails:threads:_) =
        L.map (\x -> read x :: Int) . L.take 5 $ args
      (freq:theta:sigmaGabor:sigmaGreen:angularFreq:radialFreq:alpha:_) =
        L.map (\x -> read x :: Double) . L.drop 5 $ args
      xs = L.map (\x -> read x :: (Int, Int)) . L.drop 12 $ args
      img =
        fromListUnboxed (Z :. size :. size) .
        L.map (:+ 0) .
        Arr.elems .
        accumArray (+) (0 :: Double) ((0, 0), (size - 1, size - 1)) .
        L.map (\idx -> (idx, 255)) $
        xs
      harmonics@(HarmonicArray harmonicsArray) =
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
      gaborImg = gabor size size freq (theta / 360 * 2 * pi) sigmaGabor
      gaborImgF = projectOntoBasis gaborImg harmonics harmonicsC
      folderName = "GaborTest"
  createDirectoryIfMissing True "GaborTest"
  arrG <-
    (computeS . R.map (:+ 0) . sumS) <$>
    solveMonteCarloR2S1
      threads
      trails
      size
      orientations
      sigmaGreen
      len
      (0, 0, 0 / 360 * 2 * pi, 1)
  plotImageRepa (folderName </> "green.png") .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . VU.map magnitude . toUnboxed $
    arrG
  plotImageRepa (folderName </> "gabor.png") .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . VU.map realPart . toUnboxed $
    gaborImg
  plotImageRepa (folderName </> "image.png") .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . VU.map magnitude . toUnboxed $
    img
  plan <- generateDFTPlan getEmptyPlan gaborImgF
  -- convolvedGabor <- convolveGabor plan img gaborImgF
  -- projectedGabor <- projectFilterP (computeS . R.map (\x -> realPart x :+ 0) $ gaborImg) (HarmonicArray gaborImgF)
  projectedGabor <- projectFilterP img (HarmonicArray gaborImgF)
  projectedFilter <- projectFilterP arrG harmonics
  convolved <-
    recoverFilterP
      (computeS $
       R.zipWith (\a b -> a * conjugate b) projectedGabor projectedFilter)
      harmonicsC
  -- convolved <- recoverFilterP projectedGabor harmonicsC
  plotImageRepa (folderName </> ("GaborRecon_" L.++ show theta L.++ ".png")) .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . VU.map magnitude . toUnboxed $
    convolved
