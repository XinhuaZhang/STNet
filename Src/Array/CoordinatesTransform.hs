{-# LANGUAGE FlexibleContexts #-}

module Src.Array.CoordinatesTransform
  ( module Src.Array.CoordinatesTransform
  , module Src.Array.TwoHalfDArray
  ) where

import           Data.Array.Repa         as R
import           Data.Complex
import           Data.Vector.Unboxed     as VU
import           Src.Array.Transform
import           Src.Array.TwoHalfDArray
import           Src.Utils.Coordinates

{-# INLINE cartesian2polar2D #-}

cartesian2polar2D
  :: (R.Source s Double)
  => Int
  -> Int
  -> (Double, Double)
  -> Double
  -> (Double, Double)
  -> R.Array s DIM2 Double
  -> R.Array D DIM2 Double
cartesian2polar2D ts rs (cRow, cCol) polarR valRange arr =
  fromFunction
    (Z :. ts :. rs)
    (\(Z :. t :. r) ->
       let row =
             cRow +
             (deltaR * fromIntegral r) * cos (deltaTheta * fromIntegral t)
           col =
             cCol +
             (deltaR * fromIntegral r) * sin (deltaTheta * fromIntegral t)
       in (bicubicInterpolation ds valRange (row, col)) / (fromIntegral ts) * 2 *
          pi *
          deltaR *
          fromIntegral r)
  where
    ds = computeDerivativeS . computeS . delay $ arr
    deltaTheta = 2 * pi / fromIntegral ts
    deltaR = polarR / fromIntegral rs

{-# INLINE polar2Cartesian2D #-}

polar2Cartesian2D
  :: (R.Source s Double)
  => Int
  -> Int
  -> Double
  -> (Double, Double)
  -> (Double, Double)
  -> R.Array s DIM2 Double
  -> R.Array D DIM2 Double
polar2Cartesian2D rows cols polarR valRange (centerR, centerC) arr =
  let ds = computeDerivativeS . computeS . delay $ arr
      (Z :. ts :. rs) = extent arr
      deltaTheta = 2 * pi / fromIntegral ts
      deltaR = polarR / fromIntegral rs
  in fromFunction
       (Z :. rows :. cols)
       (\(Z :. r :. c) ->
          let x = (fromIntegral r) - centerR
              y = (fromIntegral c) - centerC
              theta = angleFunctionRad x y
              radius = sqrt $ x ^ (2 :: Int) + y ^ (2 :: Int)
          in bicubicInterpolation
               ds
               valRange
               (theta / deltaTheta, radius / deltaR))

{-# INLINE cartesian2polar #-}

cartesian2polar
  :: (Source s Double)
  => Int
  -> Int
  -> (Double, Double)
  -> Double
  -> (Double, Double)
  -> Array s DIM3 Double
  -> Array U DIM3 Double
cartesian2polar ts rs center polarR valRange =
  mapArray (computeS . cartesian2polar2D ts rs center polarR valRange)


{-# INLINE polar2Cartesian #-}

polar2Cartesian
  :: (R.Source s Double)
  => Int
  -> Int
  -> Double
  -> (Double, Double)
  -> (Double, Double)
  -> Array s DIM3 Double
  -> Array U DIM3 Double 
polar2Cartesian rows cols polarR valRange center =
  mapArray (computeS . polar2Cartesian2D rows cols polarR valRange center)


{-# INLINE cartesian2polarC #-}

cartesian2polarC
  :: (R.Source s (Complex Double))
  => Int
  -> Int
  -> (Double, Double)
  -> Double
  -> (Double, Double)
  -> Array s DIM3 (Complex Double)
  -> Array U DIM3 (Complex Double)
cartesian2polarC ts rs center polarR valRange arr =
  let (x, y) =
        VU.unzip . VU.map (\(a :+ b) -> (a, b)) . toUnboxed . computeS . delay $
        arr
      arrX = fromUnboxed (get25DArrayShape arr) $ x
      arrY = fromUnboxed (get25DArrayShape arr) $ y
  in computeS $
     R.zipWith
       (:+)
       (mapArray
          (computeS . cartesian2polar2D ts rs center polarR valRange)
          arrX)
       (mapArray
          (computeS . cartesian2polar2D ts rs center polarR valRange)
          arrY)

{-# INLINE polar2CartesianC #-}

polar2CartesianC
  :: (R.Source s (Complex Double))
  => Int
  -> Int
  -> Double
  -> (Double, Double)
  -> (Double, Double)
  -> Array s DIM3 (Complex Double)
  -> Array U DIM3 (Complex Double) 
polar2CartesianC rows cols polarR valRange center arr =
  let (x, y) =
        VU.unzip . VU.map (\(a :+ b) -> (a, b)) . toUnboxed . computeS . delay $
        arr
      arrX = fromUnboxed (get25DArrayShape arr) $ x
      arrY = fromUnboxed (get25DArrayShape arr) $ y
  in computeS $
     R.zipWith
       (:+)
       (mapArray
          (computeS . polar2Cartesian2D rows cols polarR valRange center)
          arrX)
       (mapArray
          (computeS . polar2Cartesian2D rows cols polarR valRange center)
          arrY)
