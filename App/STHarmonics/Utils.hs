module App.STHarmonics.Utils where

import           Control.Monad        as M
import           Data.Array.Repa      as R
import           Data.Binary
import           Data.ByteString      as BS
import           Data.ByteString.Lazy as BL
import           Data.Int
import           Data.List            as L
import           GHC.Generics
import           System.IO

writeRepaArray :: FilePath -> R.Array U DIM3 Double -> IO ()
writeRepaArray filePath arr =
  withBinaryFile filePath WriteMode $ \h -> do
    let x = encode (listOfShape . extent $ arr, R.toList $ arr)
        len = fromIntegral . BL.length $ x :: Word64
    BL.hPut h . encode $ len
    BL.hPut h x

readRepaArray :: FilePath -> IO (R.Array U DIM3 Double)
readRepaArray filePath =
  withBinaryFile filePath ReadMode $ \h -> do
    lenBS <- BL.hGet h 8
    let len = fromIntegral (decode lenBS :: Int64) :: Int
    bs <- BL.hGet h len
    return . uncurry (\xs ys -> fromListUnboxed (shapeOfList xs) ys) . decode $
      bs
