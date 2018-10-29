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
  let (orientations:size:len:trails:threads:contrastN:_) =
        L.map (\x -> read x :: Int) . L.take 6 $ args
      (sigma:_) = L.map (\x -> read x :: Double) . L.drop 6 $ args
      folderName = "GreensData"
  arrG <-
    solveMonteCarlo'
      threads
      trails
      size
      orientations
      sigma
      len
      (0, 0, 0 / 360 * 2 * pi)
  createDirectoryIfMissing True folderName
  writeRepaArray
    (folderName </>
     ("GeneralizedGreens_" L.++ show size L.++ "_" L.++ show orientations L.++ "_" L.++
      show sigma L.++
      ".dat"))
    arrG
