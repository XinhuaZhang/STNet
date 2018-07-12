{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Src.Array.Transform
  ( module Src.Array.Transform
  , module Src.Array.TwoHalfDArray
  ) where

import           Data.Array.Repa              as R
import           Data.Array.Repa.Stencil      as R
import           Data.Array.Repa.Stencil.Dim2 as R
import           Data.Complex
import           Data.List                    as L
import           Data.Vector                  as V
import           Data.Vector.Unboxed          as VU
import           Foreign.Storable
import           Prelude                      as P
import           Src.Array.TwoHalfDArray
import           Src.Utils.Coordinates

-- factor = 2^n, n = 0,1,..
-- the first factor in the list corresponds to the inner-most (right-most) dimension.

{-# INLINE downsample #-}

downsample
  :: (Source s e, Shape sh)
  => [Int] -> R.Array s sh e -> R.Array D sh e
downsample factorList arr
  | L.all (== 1) factorList = delay arr
  | L.any (< 1) newDList = error "Downsample factors are too large."
  | otherwise =
    R.backpermute
      (shapeOfList newDList)
      (shapeOfList . L.zipWith (*) factorList . listOfShape)
      arr
  where
    dList = listOfShape . extent $ arr
    newDList = L.zipWith div dList factorList

{-# INLINE downsampleUnsafe #-}

downsampleUnsafe
  :: (Source s e, Shape sh)
  => [Int] -> R.Array s sh e -> R.Array D sh e
downsampleUnsafe factorList arr =
  R.backpermute newSh (shapeOfList . L.zipWith (*) factorList . listOfShape) arr
  where
    dList = listOfShape $ extent arr
    newSh = shapeOfList $ L.zipWith div dList factorList

{-# INLINE crop #-}

crop
  :: (Source s e, Shape sh)
  => [Int] -> [Int] -> R.Array s sh e -> R.Array D sh e
crop start len arr
  | L.any (< 0) start ||
      L.or (L.zipWith3 (\x y z -> x > (z - y)) start len dList) =
    error $
    "Crop out of boundary!\n" L.++ show start L.++ "\n" L.++ show len L.++ "\n" L.++
    show dList
  | L.length start /= L.length len || L.length start /= L.length dList =
    error $
    "crop: dimension error. \n start: " L.++ show (L.length start) L.++ " len:" L.++
    show (L.length len) L.++
    " arr:" L.++
    show (L.length dList)
  | otherwise =
    R.backpermute
      (shapeOfList len)
      (shapeOfList . L.zipWith (+) start . listOfShape)
      arr
  where
    dList = listOfShape $ extent arr

{-# INLINE cropUnsafe  #-}

cropUnsafe
  :: (Source s e, Shape sh)
  => [Int] -> [Int] -> R.Array s sh e -> R.Array D sh e
cropUnsafe start len =
  R.backpermute
    (shapeOfList len)
    (shapeOfList . L.zipWith (+) start . listOfShape)

{-# INLINE pad #-}

pad
  :: (Real e, Source s e, Shape sh)
  => [Int] -> e -> R.Array s sh e -> R.Array D sh e
pad newDims padVal arr
  | L.all (== 0) diff = delay arr
  | otherwise =
    backpermuteDft
      (fromFunction (shapeOfList dimList) (const padVal))
      (\sh' ->
         let idx = L.zipWith (-) (listOfShape sh') diff
         in if L.or (L.zipWith (\i j -> i < 0 || (i >= j)) idx oldDimList)
              then Nothing
              else Just $ shapeOfList idx)
      arr
  where
    oldDimList = listOfShape . extent $ arr
    dimList = L.zipWith max newDims oldDimList
    diff =
      L.zipWith
        (\a b ->
           if a - b <= 0
             then 0
             else if a - b == 1
                    then 1
                    else div (a - b) 2)
        newDims
        oldDimList


-- 2.5D Array Operations

computeDerivativeP
  :: R.Array U DIM2 Double -> IO [R.Array U DIM2 Double]
computeDerivativeP arr = do
  let xStencil =
        makeStencil2 3 3 $ \ix ->
          case ix of
            Z :. 0 :. (-1) -> Just (-1)
            Z :. 0 :. 1 -> Just 1
            _ -> Nothing
      yStencil =
        makeStencil2 3 3 $ \ix ->
          case ix of
            Z :. -1 :. 0 -> Just (-1)
            Z :. 1 :. 0 -> Just 1
            _ -> Nothing
      xyStencil =
        makeStencil2 3 3 $ \ix ->
          case ix of
            Z :. -1 :. -1 -> Just 1
            Z :. -1 :. 1 -> Just (-1)
            Z :. 1 :. -1 -> Just (-1)
            Z :. 1 :. 1 -> Just 1
            _ -> Nothing
      ds =
        L.map
          (\s -> R.map (/ 2) $ mapStencil2 BoundClamp s arr)
          [xStencil, yStencil, xyStencil]
  ds' <- P.mapM computeP ds
  return $! (arr : ds')


{-# INLINE computeDerivativeS #-}

computeDerivativeS
  :: R.Array U DIM2 Double -> [R.Array U DIM2 Double]
computeDerivativeS arr = arr : ds'
  where
    xStencil =
      makeStencil2 3 3 $ \ix ->
        case ix of
          Z :. 0 :. (-1) -> Just (-1)
          Z :. 0 :. 1 -> Just 1
          _ -> Nothing
    yStencil =
      makeStencil2 3 3 $ \ix ->
        case ix of
          Z :. -1 :. 0 -> Just (-1)
          Z :. 1 :. 0 -> Just 1
          _ -> Nothing
    xyStencil =
      makeStencil2 3 3 $ \ix ->
        case ix of
          Z :. -1 :. -1 -> Just 1
          Z :. -1 :. 1 -> Just (-1)
          Z :. 1 :. -1 -> Just (-1)
          Z :. 1 :. 1 -> Just 1
          _ -> Nothing
    ds =
      L.map
        (\s -> R.map (/ 2) $ mapStencil2 BoundClamp s arr)
        [xStencil, yStencil, xyStencil]
    ds' = L.map computeS ds

{-# INLINE bicubicInterpolation #-}

bicubicInterpolation
  :: [R.Array U DIM2 Double] -> (Double,Double) ->  (Double,Double) -> Double
bicubicInterpolation ds (minVal, maxVal) (y, x)
  | (x < 0) ||
      (x > (fromIntegral nx - 1)) || (y < 0) || (y > (fromIntegral ny - 1)) = 0
  | result < minVal = minVal
  | result > maxVal = maxVal
  | otherwise = result
  where
    (Z :. ny :. nx) = extent . P.head $ ds
    x' = x - fromIntegral (floor x :: Int)
    y' = y - fromIntegral (floor y :: Int)
    idx =
      VU.fromListN
        4
        [ (floor y, floor x)
        , (floor y, ceiling x)
        , (ceiling y, floor x)
        , (ceiling y, ceiling x)
        ] :: VU.Vector (Int, Int)
    xs =
      VU.concat .
      P.map (\arr' -> VU.map (\(i, j) -> arr' R.! (Z :. i :. j)) idx) $
      ds
    alpha = V.map (VU.sum . VU.zipWith (*) xs) matrixA
    arr =
      fromListUnboxed (Z :. 4 :. 4) . V.toList $ alpha :: R.Array U DIM2 Double
    arr1 =
      R.traverse arr id (\f idx'@(Z :. j :. i) -> f idx' * (x' ^ i) * (y' ^ j))
    result = sumAllS arr1
    matrixA =
      V.fromListN 16 . P.map (VU.fromListN 16) $
      [ [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
      , [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
      , [-3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
      , [2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
      , [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
      , [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
      , [0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0]
      , [0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0]
      , [-3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0]
      , [0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0]
      , [9, -9, -9, 9, 6, 3, -6, -3, 6, -6, 3, -3, 4, 2, 2, 1]
      , [-6, 6, 6, -6, -3, -3, 3, 3, -4, 4, -2, 2, -2, -2, -1, -1]
      , [2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0]
      , [0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0]
      , [-6, 6, 6, -6, -4, -2, 4, 2, -3, 3, -3, 3, -2, -1, -2, -1]
      , [4, -4, -4, 4, 2, 2, -2, -2, 2, -2, 2, -2, 1, 1, 1, 1]
      ]


{-# INLINE rescale2D #-}

rescale2D
  :: (Source s Double)
  => (Int, Int)
  -> (Double, Double)
  -> R.Array s DIM2 Double
  -> R.Array D DIM2 Double
rescale2D (newNy, newNx) bound arr =
  fromFunction
    (Z :. newNy :. newNx)
    (\(Z :. j :. i) ->
       bicubicInterpolation
         ds
         bound
         (fromIntegral j * ratioY, fromIntegral i * ratioX))
  where
    ds = computeDerivativeS . computeUnboxedS . delay $ arr
    (Z :. ny' :. nx') = extent arr
    ratioX = fromIntegral (nx' - 1) / fromIntegral (newNx - 1)
    ratioY = fromIntegral (ny' - 1) / fromIntegral (newNy - 1) 

{-# INLINE rescale25D #-}

rescale25D :: (Int, Int)
           -> (Double, Double)
           -> TwoHalfDArray Double
           -> TwoHalfDArray Double
rescale25D newSize bound arr =
  mapArray (computeUnboxedS . rescale2D newSize bound) arr

{-# INLINE rescale25DC #-}

rescale25DC
  :: (Int, Int)
  -> (Double, Double)
  -> TwoHalfDArray (Complex Double)
  -> TwoHalfDArray (Complex Double)
rescale25DC newSize bound arr =
  let (x, y) =
        VU.unzip . VU.map (\(a :+ b) -> (a, b)) . toUnboxed . get25DArray $ arr
      arrX = TwoHalfDArray . fromUnboxed (get25DArrayShape arr) $ x
      arrY = TwoHalfDArray . fromUnboxed (get25DArrayShape arr) $ y
  in zipWith25DArray
       (:+)
       (rescale25D newSize bound arrX)
       (rescale25D newSize bound arrY)

-- Set the maximum value of the maximum size, the ratio is intact.

{-# INLINE resize2D #-}   

resize2D
  :: (R.Source s Double)
  => Int -> (Double, Double) -> Array s DIM2 Double -> Array D DIM2 Double
resize2D n (minVal, maxVal) arr =
  fromFunction (Z :. newNy :. newNx) $ \(Z :. j :. i) ->
    bicubicInterpolation
      ds
      (minVal, maxVal)
      (fromIntegral j * ratioY, fromIntegral i * ratioX)
  where
    (Z :. ny :. nx) = extent arr
    newNy =
      if ny >= nx
        then n
        else round
               (fromIntegral n * fromIntegral ny / fromIntegral nx :: Double)
    newNx =
      if ny >= nx
        then round
               (fromIntegral n * fromIntegral nx / fromIntegral ny :: Double)
        else n
    ratioX = fromIntegral (nx - 1) / fromIntegral (newNx - 1)
    ratioY = fromIntegral (ny - 1) / fromIntegral (newNy - 1)
    ds = computeDerivativeS . computeUnboxedS . delay $ arr

{-# INLINE resize25D #-}

resize25D :: Int
          -> (Double, Double)
          -> TwoHalfDArray Double
          -> TwoHalfDArray Double
resize25D newSize bound arr = mapArray (computeS . resize2D newSize bound) arr

{-# INLINE resize25DC #-}

resize25DC
  :: Int
  -> (Double, Double)
  -> TwoHalfDArray (Complex Double)
  -> TwoHalfDArray (Complex Double)
resize25DC newSize bound arr =
  let (x, y) =
        VU.unzip . VU.map (\(a :+ b) -> (a, b)) . toUnboxed . get25DArray $ arr
      arrX = TwoHalfDArray . fromUnboxed (get25DArrayShape arr) $ x
      arrY = TwoHalfDArray . fromUnboxed (get25DArrayShape arr) $ y
  in zipWith25DArray
       (:+)
       (resize25D newSize bound arrX)
       (resize25D newSize bound arrY)


{-# INLINE rotatePixel #-}

rotatePixel :: VU.Vector Double
            -> (Double, Double)
            -> (Double, Double)
            -> (Double, Double)
rotatePixel mat (centerY, centerX) (y, x) = (y3, x3)
  where
    x1 = x - centerX
    y1 = y - centerY
    (y2, x2) = vecMatMult (y1, x1) mat
    x3 = x2 + centerX
    y3 = y2 + centerY

{-# INLINE vecMatMult #-}

vecMatMult :: (Double, Double) -> VU.Vector Double -> (Double, Double)
vecMatMult (x, y) vec = (a * x + c * y, b * x + d * y)
  where
    a = vec VU.! 0
    b = vec VU.! 1
    c = vec VU.! 2
    d = vec VU.! 3


-- First pading a 2D array to a square array then rotating it
padResizeRotate2DArray
  :: (R.Source s Double)
  => Int -> (Double,Double) -> [Double] -> Array s DIM2 Double -> [Array U DIM2 Double]
padResizeRotate2DArray n (minVal, maxVal) degs arr =
  L.map
    (\deg ->
       computeS $
       fromFunction
         (Z :. n :. n)
         (\(Z :. j :. i) ->
            let (j', i') =
                  rotatePixel
                    (VU.fromListN 4 $
                     P.map
                       (\f -> f (deg2Rad deg))
                       [cos, sin, \x -> -(sin x), cos])
                    (center, center)
                    (fromIntegral j, fromIntegral i)
            in if j' < 0 ||
                  j' > (fromIntegral n - 1) ||
                  i' < 0 || i' > (fromIntegral n - 1)
                 then 0
                 else bicubicInterpolation
                        ds
                        (minVal, maxVal)
                        (j' * ratio, i' * ratio)))
    degs
  where
    (Z :. ny :. nx) = extent arr
    m = max nx ny
    paddedImg = pad [m, m] 0 arr
    ds = computeDerivativeS (computeUnboxedS paddedImg)
    center = fromIntegral (n - 1) / 2
    ratio = fromIntegral (m - 1) / fromIntegral (n - 1)


-- a x b x c -> c x b x a
{-# INLINE rotate3D #-}

rotate3D :: (R.Source s e) => Array s DIM3 e -> Array D DIM3 e
rotate3D arr =
  let (Z :. a :. b :. c) = extent arr
  in R.backpermute
       (Z :. c :. a :. b)
       (\(Z :. k :. i :. j) -> (Z :. i :. j :. k))
       arr
