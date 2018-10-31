{-# LANGUAGE FlexibleContexts #-}
module App.STHarmonics.Utils where

import           Control.Monad                as M
import           Data.Array.Repa              as R
import           Data.Array.Repa.Repr.Unboxed
import           Data.Binary
import           Data.ByteString.Lazy         as BL
import           Data.Int
import           Data.List                    as L
import           GHC.Generics
import           System.IO

writeRepaArray
  :: (R.Source r e, Shape sh, Binary e)
  => FilePath -> R.Array r sh e -> IO ()
writeRepaArray filePath arr =
  withBinaryFile filePath WriteMode $ \h -> do
    let x = encode (listOfShape . extent $ arr, R.toList $ arr)
        len = fromIntegral . BL.length $ x :: Word64
    BL.hPut h . encode $ len
    BL.hPut h x

readRepaArray
  :: (R.Source U e, Shape sh, Binary e, Unbox e)
  => FilePath -> IO (R.Array U sh e)
readRepaArray filePath =
  withBinaryFile filePath ReadMode $ \h -> do
    lenBS <- BL.hGet h 8
    let len = fromIntegral (decode lenBS :: Int64) :: Int
    bs <- BL.hGet h len
    return . uncurry (\xs ys -> fromListUnboxed (shapeOfList xs) ys) . decode $
      bs
