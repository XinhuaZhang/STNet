{-# LANGUAGE FlexibleContexts #-}
import           Control.Monad               as M
import           Data.Array.Repa             as R
import           Data.List                   as L
import           Src.FokkerPlanck.MonteCarlo
import           Src.Image.ImageIO
import           System.Directory
import           System.Environment
import           System.FilePath

{-# INLINE reduceContrast #-}

reduceContrast
  :: (R.Source s Double)
  => Int -> Array s DIM3 Double -> Array D DIM3 Double
reduceContrast idx arr =
  let x = L.head . L.drop idx . L.reverse . L.sort . R.toList $ arr
  in R.map
       (\y ->
          if y >= x
            then x
            else y)
       arr

main = do
  (numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:lenStr:initStr:numTrailStr:numThreadStr:_) <-
    getArgs
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
      len = read lenStr :: Int
      init = read initStr :: (Double, Double, Double, Double)
      numTrail = read numTrailStr :: Int
      numThread = read numThreadStr :: Int
  arr <-
    solveMonteCarloR2S1RP
      numThread
      numTrail
      numPoint
      numOrientation
      numScale
      thetaSigma
      scaleSigma
      maxScale
      len
      init
  let arr' =
        computeS .
        reduceContrast 10 .
        R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . R.sumS $
        arr
  createDirectoryIfMissing True "GreensR2S1RP"
  plotImageRepa "GreensR2S1RP/GreensR2S1RP.png" . Image 8 $ arr'
  M.mapM_
    (\i ->
       plotImageRepa ("GreensR2S1RP/" L.++ (show $ i + 1) L.++ ".png") .
       Image 8 .
       computeS .
       reduceContrast 10 .
       R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . R.slice arr $
       (Z :. All :. All :. All :. i))
    [0 .. numScale - 1]
