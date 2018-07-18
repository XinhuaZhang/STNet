{-# LANGUAGE FlexibleContexts #-}
module Src.Array.TwoHalfDArray where

import           Control.Monad       as M
import           Data.Array.Repa     as R
import           Data.List           as L
import           Data.Vector.Unboxed as VU

-- 2.5D array: features x rows x cols, where column index changes fastest.
-- 2D operations on this array are replicated to the feature dimension.

{-# INLINE get25DArrayFeatures #-}

get25DArrayFeatures
  :: (Source s e)
  => Array s DIM3 e -> Int
get25DArrayFeatures arr =
  let (Z :. nf :. _ :. _) = extent arr
  in nf

{-# INLINE get25DArrayRows #-}

get25DArrayRows
  :: (Source s e)
  => Array s DIM3 e -> Int
get25DArrayRows arr =
  let (Z :. _ :. rows :. _) = extent arr
  in rows

{-# INLINE get25DArrayCols #-}

get25DArrayCols
  :: (Source s e)
  => Array s DIM3 e -> Int
get25DArrayCols arr =
  let (Z :. _ :. _ :. cols) = extent arr
  in cols

{-# INLINE get25DArrayShape #-}

get25DArrayShape
  :: (Source s e)
  => Array s DIM3 e -> DIM3
get25DArrayShape arr = extent arr

{-# INLINE get25DArrayFeatureVector #-}

get25DArrayFeatureVector
  :: (Source s e, Unbox e)
  => Array s DIM3 e -> [VU.Vector e]
get25DArrayFeatureVector arr =
  [ toUnboxed . computeS . R.slice arr $ (Z :. All :. i :. j)
  | i <- [0 .. get25DArrayRows arr - 1]
  , j <- [0 .. get25DArrayCols arr - 1]
  ]

{-# INLINE get25DArrayFeatureList #-}

get25DArrayFeatureList
  :: (Source s e)
  => Array s DIM3 e -> [[e]]
get25DArrayFeatureList arr =
  [ R.toList . R.slice arr $ (Z :. All :. i :. j)
  | i <- [0 .. get25DArrayRows arr - 1]
  , j <- [0 .. get25DArrayCols arr - 1]
  ]
  
{-# INLINE get25DArray2DVector #-}

get25DArray2DVector
  :: (Source s e, Unbox e)
  => Array s DIM3 e -> [VU.Vector e]
get25DArray2DVector arr =
  L.map
    (\k -> toUnboxed . computeS . R.slice arr $ (Z :. k :. All :. All))
    [0 .. (get25DArrayFeatures arr) - 1]
    

{-# INLINE get25DArray2DList #-}

get25DArray2DList
  :: (Source s e)
  => Array s DIM3 e -> [[e]]
get25DArray2DList arr =
  L.map
    (\k -> R.toList . R.slice arr $ (Z :. k :. All :. All))
    [0 .. (get25DArrayFeatures arr) - 1]

-- Replicating a 2D operation on the feature dimension.

{-# INLINE mapArray #-}

mapArray
  :: (Source s e, Unbox e)
  => (Array D DIM2 e -> Array U DIM2 e) -> Array s DIM3 e -> Array U DIM3 e
mapArray f arr =
  fromUnboxed (get25DArrayShape arr) .
  VU.concat . L.map (\i -> toUnboxed . f . R.slice arr $ (Z :. i :. All :. All)) $
  [0 .. get25DArrayFeatures arr - 1]


{-# INLINE mapArrayM #-}

mapArrayM
  :: (Source s e, Unbox e, Monad m)
  => (Array D DIM2 e -> m (Array U DIM2 e))
  -> Array s DIM3 e
  -> m (Array U DIM3 e)
mapArrayM f arr = do
  xs <-
    M.mapM (\i -> fmap toUnboxed . f . R.slice arr $ (Z :. i :. All :. All)) $
    [0 .. get25DArrayFeatures arr - 1]
  return . fromUnboxed (get25DArrayShape arr) . VU.concat $ xs


{-# INLINE mapVector #-}

mapVector
  :: (Source s e, Unbox e)
  => (VU.Vector e -> VU.Vector e) -> Array s DIM3 e -> Array U DIM3 e
mapVector f arr =
  fromUnboxed (get25DArrayShape arr) .
  VU.concat .
  L.map (\i -> f . toUnboxed . computeS . R.slice arr $ (Z :. i :. All :. All)) $
  [0 .. get25DArrayFeatures arr - 1]


{-# INLINE mapVectorM #-}

mapVectorM
  :: (Source s e, Unbox e, Monad m)
  => (VU.Vector e -> m (VU.Vector e)) -> Array s DIM3 e -> m (Array U DIM3 e)
mapVectorM f arr = do
  xs <-
    M.mapM
      (\i -> f . toUnboxed . computeS . R.slice arr $ (Z :. i :. All :. All)) $
    [0 .. get25DArrayFeatures arr - 1]
  return . fromUnboxed (get25DArrayShape arr) . VU.concat $ xs
