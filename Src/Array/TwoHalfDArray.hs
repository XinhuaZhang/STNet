{-# LANGUAGE FlexibleContexts #-}
module Src.Array.TwoHalfDArray where

import           Control.Monad               as M
import           Control.Monad.Parallel      as MP
import           Control.Parallel.Strategies
import           Data.Array.Repa             as R
import           Data.List                   as L
import           Data.Vector.Unboxed         as VU

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

{-# INLINE map25DArray #-}

map25DArray
  :: (Source s e, Unbox e)
  => (Array D DIM2 e -> Array U DIM2 e) -> Array s DIM3 e -> Array U DIM3 e
map25DArray f arr =
  fromUnboxed (get25DArrayShape arr) .
  VU.concat . L.map (\i -> toUnboxed . f . R.slice arr $ (Z :. i :. All :. All)) $
  [0 .. get25DArrayFeatures arr - 1]

{-# INLINE map25DArrayP #-}

map25DArrayP
  :: (Source s e, Unbox e)
  => (Array D DIM2 e -> Array U DIM2 e) -> Array s DIM3 e -> Array U DIM3 e
map25DArrayP f arr =
  fromUnboxed (get25DArrayShape arr) .
  VU.concat .
  parMap rdeepseq (\i -> toUnboxed . f . R.slice arr $ (Z :. i :. All :. All)) $
  [0 .. get25DArrayFeatures arr - 1]


{-# INLINE map25DArrayM #-}

map25DArrayM
  :: (Source s e, Unbox e, Monad m)
  => (Array D DIM2 e -> m (Array U DIM2 e))
  -> Array s DIM3 e
  -> m (Array U DIM3 e)
map25DArrayM f arr = do
  xs <-
    M.mapM (\i -> fmap toUnboxed . f . R.slice arr $ (Z :. i :. All :. All)) $
    [0 .. get25DArrayFeatures arr - 1]
  return . fromUnboxed (get25DArrayShape arr) . VU.concat $ xs
  
{-# INLINE map25DArrayMP #-}

map25DArrayMP
  :: (Source s e, Unbox e, MonadParallel m)
  => (Array D DIM2 e -> m (Array U DIM2 e))
  -> Array s DIM3 e
  -> m (Array U DIM3 e)
map25DArrayMP f arr = do
  xs <-
    MP.mapM (\i -> fmap toUnboxed . f . R.slice arr $ (Z :. i :. All :. All)) $
    [0 .. get25DArrayFeatures arr - 1]
  return . fromUnboxed (get25DArrayShape arr) . VU.concat $ xs


{-# INLINE map25DVector #-}

map25DVector
  :: (Source s e, Unbox e)
  => (VU.Vector e -> VU.Vector e) -> Array s DIM3 e -> Array U DIM3 e
map25DVector f arr =
  fromUnboxed (get25DArrayShape arr) .
  VU.concat .
  L.map (\i -> f . toUnboxed . computeS . R.slice arr $ (Z :. i :. All :. All)) $
  [0 .. get25DArrayFeatures arr - 1]
  
{-# INLINE map25DVectorP #-}

map25DVectorP
  :: (Source s e, Unbox e)
  => (VU.Vector e -> VU.Vector e) -> Array s DIM3 e -> Array U DIM3 e
map25DVectorP f arr =
  fromUnboxed (get25DArrayShape arr) .
  VU.concat .
  parMap
    rdeepseq
    (\i -> f . toUnboxed . computeS . R.slice arr $ (Z :. i :. All :. All)) $
  [0 .. get25DArrayFeatures arr - 1]


{-# INLINE map25DVectorM #-}

map25DVectorM
  :: (Source s e, Unbox e, Monad m)
  => (VU.Vector e -> m (VU.Vector e)) -> Array s DIM3 e -> m (Array U DIM3 e)
map25DVectorM f arr = do
  xs <-
    M.mapM
      (\i -> f . toUnboxed . computeS . R.slice arr $ (Z :. i :. All :. All)) $
    [0 .. get25DArrayFeatures arr - 1]
  return . fromUnboxed (get25DArrayShape arr) . VU.concat $ xs
  
{-# INLINE map25DVectorMP #-}

map25DVectorMP
  :: (Source s e, Unbox e, MonadParallel m)
  => (VU.Vector e -> m (VU.Vector e)) -> Array s DIM3 e -> m (Array U DIM3 e)
map25DVectorMP f arr = do
  xs <-
    MP.mapM
      (\i -> f . toUnboxed . computeS . R.slice arr $ (Z :. i :. All :. All)) $
    [0 .. get25DArrayFeatures arr - 1]
  return . fromUnboxed (get25DArrayShape arr) . VU.concat $ xs


{-# INLINE zipWith25D #-}

zipWith25D
  :: (Source s1 a, Source s2 b, Unbox a, Unbox b, Unbox c)
  => (a -> b -> c) -> Array s1 DIM3 a -> Array s2 DIM2 b -> Array D DIM3 c
zipWith25D f arr1 arr2 =
  R.traverse2
    arr1
    arr2
    const
    (\f3D f2D (Z :. k :. i :. j) ->
       f (f3D (Z :. k :. i :. j)) (f2D (Z :. i :. j)))
