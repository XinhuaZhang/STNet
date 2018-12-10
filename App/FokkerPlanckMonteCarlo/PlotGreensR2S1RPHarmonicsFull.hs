import           Control.Monad               as M
import           Control.Monad.Parallel      as MP
import           Data.Array                  as Arr
import           Data.Array.Repa             as R
import           Data.Complex
import qualified Data.Image                  as IM
import           Data.List                   as L
import           Src.FokkerPlanck.MonteCarlo
import           Src.Image.ImageIO
import           System.Directory
import           System.Environment
import           System.FilePath

{-# INLINE norm #-}
norm ::  [Complex Double] -> [Complex Double]
norm vec =
  let maxV = L.maximum . L.map magnitude $ vec
      minV = L.minimum . L.map magnitude $ vec
   in L.map
        (uncurry mkPolar .
         (\(m, p) -> (1, p)) . polar)
        vec

main = do
  args@(numPointStr:numOrientationStr:angularFreqStr:angularFreq1Str:radialFreqStr:radialFreq1Str:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:lenStr:initStr:numTrailStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      angularFreq = read angularFreqStr :: Int
      angularFreq1 = read angularFreq1Str :: Int
      radialFreq = read radialFreqStr :: Int
      radialFreq1 = read radialFreq1Str :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
      tao = read taoStr :: Double
      len = read lenStr :: Int
      init = read initStr :: (Double, Double, Double, Double)
      numTrail = read numTrailStr :: Int
      numThread = read numThreadStr :: Int
      angularFreqs = [-angularFreq .. angularFreq]
      radialFreqs = [-radialFreq .. radialFreq]
  arr <-
    solveMonteCarloR2S1RP''
      numThread
      numTrail
      numPoint
      angularFreqs
      [-angularFreq1 .. angularFreq1]
      radialFreqs
      [-radialFreq1 .. radialFreq1]
      thetaSigma
      scaleSigma
      maxScale
      tao
      len
      init
  let arr' = R.sumS . R.sumS $ arr :: R.Array U DIM2 (Complex Double)
      folderName = "GreensR2S1RPHarmonicsFull"
  removePathForcibly folderName
  createDirectoryIfMissing True folderName
  MP.mapM_
    (\(m, n) ->
       IM.writeImage
         (folderName </> thetaSigmaStr L.++ "_" L.++
          (show $ m + L.last angularFreqs) L.++
          "_" L.++
          show m L.++
          "_" L.++
          scaleSigmaStr L.++
          "_" L.++
          (show $ n + L.last radialFreqs) L.++
          "_" L.++
          show n L.++
          ".ppm")
         (IM.arrayToImage .
          listArray ((0, 0), (numPoint - 1, numPoint - 1)) .
          R.toList . R.slice arr $
          (Z :. All :. All :. (m + L.last angularFreqs) :.
           (n + L.last radialFreqs)) :: IM.ComplexImage))
    [(m, n) | m <- angularFreqs, n <- radialFreqs]
