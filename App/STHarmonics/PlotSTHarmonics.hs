import           App.STHarmonics.Generator
import           Data.Array.Repa           as R
import           Data.List                 as L
import           Src.Image.Utils
import           System.Directory
import           System.Environment
import           System.FilePath

main = do
  args <- getArgs
  let (rows:cols:numTheta:anglarFreq:radialFreq:_) =
        L.map (\x -> read x :: Int) . L.take 5 $ args
      (alpha:_) = L.map (\x -> read x :: Double) . L.take 1 . L.drop 5 $ args
      (inverseFlag:_) =
        L.map (\x -> read x :: STHarmonicFlag) . L.take 1 . L.drop 6 $ args
      folderPath = L.last args
      arr =
        generateSTHarmonicArray
          numTheta
          rows
          cols
          anglarFreq
          radialFreq
          alpha
          inverseFlag
  removePathForcibly folderPath
  createDirectoryIfMissing True folderPath
  writeComplexImageRepa folderPath arr
