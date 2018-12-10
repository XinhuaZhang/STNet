import           Control.Monad               as M
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
  (numPointStr:numOrientationStr:freqStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:lenStr:initStr:numTrailStr:numThreadStr:_) <-
    getArgs
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      freq = read freqStr :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
      tao = read taoStr :: Double
      len = read lenStr :: Int
      init = read initStr :: (Double, Double, Double, Double)
      numTrail = read numTrailStr :: Int
      numThread = read numThreadStr :: Int
  arr <-
    solveMonteCarloR2S1'
      numThread
      numTrail
      numPoint
      numOrientation
      freq
      thetaSigma
      tao
      len
      init
  let arr' = R.sumS arr :: R.Array U DIM2 (Complex Double)
      arr'' = R.map (\x -> x * exp (0 :+ fromIntegral freq * pi / 4)) arr
      folderName = "GreensR2S1Harmonics"
  removePathForcibly folderName
  createDirectoryIfMissing True folderName
  IM.writeImage
    (folderName </> "GreensR2S1Harmonics_" L.++ show (numOrientation + 1) L.++
     ".ppm") $
    (IM.arrayToImage .
     listArray ((0, 0), (numPoint - 1, numPoint - 1)) . R.toList $
     arr' :: IM.ComplexImage)
  IM.writeImage
    (folderName </> "GreensR2S1Harmonics_0_" L.++ show (numOrientation + 1) L.++
     ".ppm") $
    (IM.arrayToImage .
     listArray ((0, 0), (numPoint - 1, numPoint - 1)) . R.toList . R.sumS $
     arr'' :: IM.ComplexImage)
  M.mapM_
    (\i -> do
       IM.writeImage
         (folderName </> folderName L.++ "_" L.++ ((show $ i + 1) L.++ ".ppm"))
         (IM.arrayToImage .
          listArray ((0, 0), (numPoint - 1, numPoint - 1)) .
          R.toList . R.slice arr $
          (Z :. All :. All :. i) :: IM.ComplexImage)
       IM.writeImage
         (folderName </> folderName L.++ "_0_" L.++ ((show $ i + 1) L.++ ".ppm"))
         (IM.arrayToImage .
          listArray ((0, 0), (numPoint - 1, numPoint - 1)) .
          R.toList . R.slice arr'' $
          (Z :. All :. All :. i) :: IM.ComplexImage))
    [0 .. numOrientation - 1]
