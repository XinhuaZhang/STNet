{-# LANGUAGE FlexibleContexts #-}
module Src.Image.Utils where

import           Control.Monad       as M
import           Data.Array
import           Data.Array.Repa     as R
import           Data.Complex        as C
import           Data.Image          as Img
import           Data.List           as L
import           Data.Vector.Unboxed as VU
import           System.FilePath

-- Because the element is complex valued, the Repa array can only be 2D

writeComplexImageRepa2D
  :: (R.Source s (C.Complex Double))
  => FilePath -> R.Array s DIM2 (C.Complex Double) -> IO ()
writeComplexImageRepa2D filePath arr = do
  let (Z :. rows :. cols) = extent arr
      img =
        arrayToImage . listArray ((0, 0), (rows - 1, cols - 1)) . R.toList $ arr :: ComplexImage
  writeImage filePath img

writeComplexImageRepa
  :: (R.Source s (C.Complex Double))
  => FilePath -> R.Array s DIM3 (C.Complex Double) -> IO ()
writeComplexImageRepa folderPath arr = do
  let (Z :. nf :. _ :. _) = extent arr
  M.mapM_
    (\i ->
       writeComplexImageRepa2D (folderPath </> ((show i) L.++ ".ppm")) .
       R.slice arr $
       (Z :. (i - 1) :. All :. All))
    [1 .. nf]
