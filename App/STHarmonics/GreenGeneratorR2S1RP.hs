import           App.STHarmonics.Utils
import           Control.Monad               as M
import           Data.Array                  as Arr
import           Data.Array.Repa             as R
import           Data.Complex
import           Data.List                   as L
import           Data.Vector.Storable        as VS
import           Data.Vector.Unboxed         as VU
import           Src.Array.Transform         (rotate3D)
import           Src.Array.TwoHalfDArray
import           Src.FokkerPlanck.MonteCarlo
import           Src.Image.ImageIO
import           System.Directory
import           System.Environment
import           System.FilePath
import           System.Random

main = do
  args <- getArgs
  let (orientations:scales:size:len:trails:threads:contrastN:_) =
        L.map (\x -> read x :: Int) . L.take 7 $ args
      (thetaSigma:scaleSigma:maxScale:_) =
        L.map (\x -> read x :: Double) . L.drop 7 $ args
      folderName = "GreensData"
  arrG <-
    solveMonteCarloR2S1RP
      threads
      trails
      size
      orientations
      scales
      thetaSigma
      scaleSigma
      maxScale
      len
      (0, 0, 0 / 360 * 2 * pi, 1)
  createDirectoryIfMissing True folderName
  writeRepaArray
    (folderName </>
     ("Greens_" L.++ show size L.++ "_" L.++ show orientations L.++ "_" L.++
      show scales L.++
      "_" L.++
      show thetaSigma L.++
      "_" L.++
      show scaleSigma L.++
      "_" L.++
      show maxScale L.++
      ".dat"))
    arrG
