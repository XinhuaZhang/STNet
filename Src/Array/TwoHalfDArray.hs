{-# LANGUAGE FlexibleContexts #-}
module Src.Array.TwoHalfDArray where

import           Control.Monad       as M
import           Data.Array.Repa     as R
import           Data.List           as L
import           Data.Vector.Unboxed as VU

-- 2.5D array: features x rows x cols, where column index changes fastest.
-- 2D operations on this array are replicated to the feature dimension.

newtype TwoHalfDArray a =
  TwoHalfDArray (Array U DIM3 a)

{-# INLINE get25DArrayFeatures #-}

get25DArrayFeatures
  :: (Unbox a)
  => TwoHalfDArray a -> Int
get25DArrayFeatures (TwoHalfDArray arr) =
  let (Z :. nf :. _ :. _) = extent arr
  in nf

{-# INLINE get25DArrayRows #-}

get25DArrayRows
  :: (Unbox a)
  => TwoHalfDArray a -> Int
get25DArrayRows (TwoHalfDArray arr) =
  let (Z :. _ :. rows :. _) = extent arr
  in rows

{-# INLINE get25DArrayCols #-}

get25DArrayCols
  :: (Unbox a)
  => TwoHalfDArray a -> Int
get25DArrayCols (TwoHalfDArray arr) =
  let (Z :. _ :. _ :. cols) = extent arr
  in cols

{-# INLINE get25DArrayShape #-}

get25DArrayShape
  :: (Unbox a)
  => TwoHalfDArray a -> DIM3
get25DArrayShape (TwoHalfDArray arr) =
   extent arr


{-# INLINE get25DArray #-}

get25DArray :: TwoHalfDArray a -> Array U DIM3 a
get25DArray (TwoHalfDArray arr) = arr

{-# INLINE get25DArrayFeatureVector #-}

get25DArrayFeatureVector
  :: (Unbox a)
  => TwoHalfDArray a -> [VU.Vector a]
get25DArrayFeatureVector arr =
  [ toUnboxed . computeUnboxedS . R.slice (get25DArray arr) $
  (Z :. All :. i :. j)
  | i <- [0 .. get25DArrayRows arr - 1]
  , j <- [0 .. get25DArrayCols arr - 1]
  ]

{-# INLINE get25DArrayFeatureList #-}

get25DArrayFeatureList
  :: (Unbox a)
  => TwoHalfDArray a -> [[a]]
get25DArrayFeatureList arr =
  [ R.toList . computeUnboxedS . R.slice (get25DArray arr) $
  (Z :. All :. i :. j)
  | i <- [0 .. get25DArrayRows arr - 1]
  , j <- [0 .. get25DArrayCols arr - 1]
  ]
  
{-# INLINE get25DArray2DVector #-}

get25DArray2DVector
  :: (Unbox a)
  => TwoHalfDArray a -> [VU.Vector a]
get25DArray2DVector arr =
  L.map
    (\k ->
       toUnboxed . computeS . R.slice (get25DArray arr) $ (Z :. k :. All :. All))
    [0 .. (get25DArrayFeatures arr) - 1]
    

{-# INLINE get25DArray2DList #-}

get25DArray2DList
  :: (Unbox a)
  => TwoHalfDArray a -> [[a]]
get25DArray2DList arr =
  L.map
    (\k -> R.toList . R.slice (get25DArray arr) $ (Z :. k :. All :. All))
    [0 .. (get25DArrayFeatures arr) - 1]
    

-- element-wise operations
   
{-# INLINE map25DArray #-}

map25DArray
  :: (Unbox a, Unbox b)
  => (a -> b) -> TwoHalfDArray a -> TwoHalfDArray b
map25DArray f (TwoHalfDArray arr) = TwoHalfDArray . computeS . R.map f $ arr

{-# INLINE zipWith25DArray #-}

zipWith25DArray
  :: (Unbox a, Unbox b, Unbox c)
  => (a -> b -> c) -> TwoHalfDArray a -> TwoHalfDArray b -> TwoHalfDArray c
zipWith25DArray f (TwoHalfDArray arr1) (TwoHalfDArray arr2) =
  TwoHalfDArray . computeS $ R.zipWith f arr1 arr2


-- Replicating a 2D operation on the feature dimension.

{-# INLINE mapArray #-}

mapArray
  :: (Unbox a)
  => (Array D DIM2 a -> Array U DIM2 a) -> TwoHalfDArray a -> TwoHalfDArray a
mapArray f arr =
  TwoHalfDArray .
  fromUnboxed (get25DArrayShape arr) .
  VU.concat .
  L.map
    (\i -> toUnboxed . f . R.slice (get25DArray arr) $ (Z :. i :. All :. All)) $
  [0 .. get25DArrayFeatures arr - 1]


{-# INLINE mapArrayM #-}

mapArrayM
  :: (Unbox a, Monad m)
  => (Array D DIM2 a -> m (Array U DIM2 a))
  -> TwoHalfDArray a
  -> m (TwoHalfDArray a)
mapArrayM f arr = do
  xs <-
    M.mapM
      (\i ->
         fmap toUnboxed . f . R.slice (get25DArray arr) $ (Z :. i :. All :. All)) $
    [0 .. get25DArrayFeatures arr - 1]
  return . TwoHalfDArray . fromUnboxed (get25DArrayShape arr) . VU.concat $ xs


{-# INLINE mapVector #-}

mapVector
  :: (Unbox a)
  => (VU.Vector a -> VU.Vector a) -> TwoHalfDArray a -> TwoHalfDArray a
mapVector f arr =
  TwoHalfDArray .
  fromUnboxed (get25DArrayShape arr) .
  VU.concat .
  L.map
    (\i ->
       f . toUnboxed . computeS . R.slice (get25DArray arr) $
       (Z :. i :. All :. All)) $
  [0 .. get25DArrayFeatures arr - 1]


{-# INLINE mapVectorM #-}

mapVectorM
  :: (Unbox a, Monad m)
  => (VU.Vector a -> m (VU.Vector a)) -> TwoHalfDArray a -> m (TwoHalfDArray a)
mapVectorM f arr = do
  xs <-
    M.mapM
      (\i ->
         f . toUnboxed . computeS . R.slice (get25DArray arr) $
         (Z :. i :. All :. All)) $
    [0 .. get25DArrayFeatures arr - 1]
  return . TwoHalfDArray . fromUnboxed (get25DArrayShape arr) . VU.concat $ xs
