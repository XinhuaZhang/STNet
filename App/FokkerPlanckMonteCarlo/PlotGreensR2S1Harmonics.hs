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
  (numPointStr:numOrientationStr:freqStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:lenStr:initStr:numTrailStr:numThreadStr:_) <-
    getArgs
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      freq = read freqStr :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
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
      len
      init
  let arr'
        -- computeS .
        -- reduceContrast 10 .
        -- R.extend (Z :. (1 :: Int) :. All :. All) .
       = R.sumS arr :: R.Array U DIM2 (Complex Double)
      folderName = "GreensR2S1Harmonics"
  createDirectoryIfMissing True folderName
  -- plotImageRepa (folderName </> "GreensR2S1RP.png") . Image 8 $ arr'
  IM.writeImage (folderName </> "GreensR2S1RP.ppm") $
    (IM.arrayToImage .
     listArray ((0, 0), (numPoint - 1, numPoint - 1)) . R.toList $
     arr' :: IM.ComplexImage)
  M.mapM_
    (\i ->
       IM.writeImage
         (folderName </> ((show $ i + 1) L.++ ".ppm"))
         (IM.arrayToImage .
          listArray ((0, 0), (numPoint - 1, numPoint - 1)) .
          R.toList . R.slice arr $
          (Z :. All :. All :. i) :: IM.ComplexImage))
    [0 .. numOrientation - 1]
