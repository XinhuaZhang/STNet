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
  let (size:_) = L.map (\x -> read x :: Int) . L.take 1 $ args
      (r:theta:sigma:angularFreq:radialFreq:alpha:_) =
        L.map (\x -> read x :: Double) . L.drop 1 $ args
      harmonics =
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
      coefficients =
        impulsiveCoefficients
          (log r)
          (theta / 360 * 2 * pi)
          (fromIntegral $ div size 2)
          [-angularFreq .. angularFreq]
          [-radialFreq .. radialFreq]
      folderName = "LocalizationTest"
      localizedImg =
        recoverFilter'
          (fromListUnboxed (Z :. (L.length coefficients)) coefficients)
          harmonicsC
  when (r <= 0) (error "r <= 0")
  plotImageRepa
    (folderName </>
     ("Impulse_" L.++ show r L.++ "_" L.++ show theta L.++ ".png")) .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) . VU.map magnitude . toUnboxed $
    localizedImg
