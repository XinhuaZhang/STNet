import           App.STWeightSharing.Convolution
import           Control.Monad                   as M
import           Control.Parallel.Strategies
import           Data.Array                      as Arr
import           Data.Array.Repa                 as R
import           Data.Complex
import           Data.List                       as L
import           Data.Vector.Storable            as VS
import           Data.Vector.Unboxed             as VU
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
  print args
  let (orientations:size:len:trails:threads:contrastN:_) =
        L.map (\x -> read x :: Int) . L.take 6 $ args
      (sigma:angularFreq':radialFreq:alpha:theta:_) =
        L.map (\x -> read x :: Double) . L.drop 6 $ args
      init = (0, 0, 0, 1)
      angularFreq = round angularFreq'
      freqs =
        generateHarmonicCoefficients
          (theta / 360 * 2 * pi)
          [-angularFreq .. angularFreq]
  arrs <-
    M.mapM
      (\f ->
         solveMonteCarloR2S1' threads trails size orientations f sigma len init)
      [-angularFreq .. angularFreq]
  let arr =
        R.sumS . L.foldl1' (R.zipWith (+)) $
        L.zipWith (\x y -> R.map (* x) y) freqs arrs
  plotImageRepa "test.png" .
    Image 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map magnitude $
    arr
