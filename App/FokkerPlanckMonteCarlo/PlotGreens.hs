import           Src.FokkerPlanck.MonteCarlo
import           Src.Image.ImageIO
import           System.Environment
import Data.Array.Repa as R

main = do
  (numPointStr:numOrientationStr:sigmaStr:lenStr:initStr:numTrailStr:numThreadStr:_) <-
    getArgs
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      sigma = read sigmaStr :: Double
      len = read lenStr :: Int
      init = read initStr :: (Double, Double, Double)
      numTrail = read numTrailStr :: Int
      numThread = read numThreadStr :: Int
  arr <-
    solveMonteCarlo numThread numTrail numPoint numOrientation sigma len init
  let arr' = computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS $ arr
  plotImageRepa "test.png" . Image 8 $ arr'
