module App.STHarmonics.R2S1RPHarmonics where

import           Data.Array.Repa                as R
import           Data.Complex
import           Data.List                      as L
import           Src.Array.CoordinatesTransform
import           Src.Utils.Coordinates

data HarmonicFlag
  = Harmonic
  | ConjugateHarmonic
  deriving (Read, Show)

data HarmoicArray sh
  = HarmoicArray (Array U sh (Complex Double))
  | ConjugateHarmoicArray (Array U sh (Complex Double))

instance Show (HarmoicArray sh) where
  show (HarmoicArray arr)          = "HarmoicArray"
  show (ConjugateHarmoicArray arr) = "ConjugateHarmoicArray"

data R2S1RPHarmonic
  = R2S1Harmonic (HarmoicArray DIM3)
  | R2S1RPHarmonic (HarmoicArray DIM4)

instance Show R2S1RPHarmonic where
  show (R2S1Harmonic x)   = "R2S1Harmonic " L.++ show x
  show (R2S1RPHarmonic x) = "R2S1RPHarmonic " L.++ show x

-- (x,y,\theta) <-> (r,\theta_1, \theta_2) <-> (r, u, v)
-- (x,y) <-> (r,\theta_1)
-- u = (\theta_1 + \theta_2) / 2, v = (\theta_1 - \theta_2) / 2
-- \theta_1 = u + v, \theta_2 = u - v


-- Create Harmonics in (x, y, \theta_2) space, assuming the center of
-- polar coordinates is the center of the 2D array.

{-# INLINE generateR2S1HarmonicArray #-}

generateR2S1HarmonicArray
  :: Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> HarmonicFlag
  -> R2S1RPHarmonic
generateR2S1HarmonicArray numTheta2 rows cols angularFreq1 angularFreq2 radialFreq alpha inverseFlag =
  let centerR = div rows 2
      centerC = div cols 2
      deltaTheta2 = 2 * pi / fromIntegral numTheta2
      sign =
        case inverseFlag of
          Harmonic -> (-1)
          ConjugateHarmonic -> 1
      arr =
        computeS $
        fromFunction
          (Z :. rows :. cols :. numTheta2)
          (\(Z :. i :. j :. k) ->
             let x = fromIntegral $ i - centerR
                 y = fromIntegral $ j - centerC
                 theta1 = angleFunctionRad x y
                 theta2 = deltaTheta2 * fromIntegral k
                 r = sqrt $ x ^ (2 :: Int) + y ^ (2 :: Int)
             in if r == 0
                  then 1
                  else (((r :+ 0) ** (alpha :+ (sign * radialFreq)))) *
                       exp
                         (0 :+
                          (sign *
                           (angularFreq1 * theta2 +
                            angularFreq2 * (theta2 - theta1)))))
  in deepSeqArray arr $
     case inverseFlag of
       Harmonic -> R2S1Harmonic . HarmoicArray $ arr
       ConjugateHarmonic -> R2S1Harmonic . ConjugateHarmoicArray $ arr
       
{-# INLINE generateR2S1HarmonicArray1 #-}

generateR2S1HarmonicArray1
  :: Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> HarmonicFlag
  -> R2S1RPHarmonic
generateR2S1HarmonicArray1 numTheta2 rows cols angularFreq1 angularFreq2 radialFreq alpha inverseFlag =
  let centerR = div rows 2
      centerC = div cols 2
      deltaTheta2 = 2 * pi / fromIntegral numTheta2
      sign =
        case inverseFlag of
          Harmonic -> (-1)
          ConjugateHarmonic -> 1
      arr =
        computeS $
        fromFunction
          (Z :. rows :. cols :. numTheta2)
          (\(Z :. i :. j :. k) ->
             let x = fromIntegral $ i - centerR
                 y = fromIntegral $ j - centerC
                 theta1 = angleFunctionRad x y
                 theta2 = deltaTheta2 * fromIntegral k
                 r = sqrt $ x ^ (2 :: Int) + y ^ (2 :: Int)
             in if r == 0
                  then 1
                  else (((r :+ 0) ** (alpha :+ (sign * radialFreq)))) *
                       exp
                         (0 :+
                          (sign *
                           (-- angularFreq1 * theta1 / 2 +
                            angularFreq2 * (theta2 - theta1)))))
  in deepSeqArray arr $
     case inverseFlag of
       Harmonic -> R2S1Harmonic . HarmoicArray $ arr
       ConjugateHarmonic -> R2S1Harmonic . ConjugateHarmoicArray $ arr


{-# INLINE generateR2S1RPHarmonicArray #-}

generateR2S1RPHarmonicArray
  :: Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> HarmonicFlag
  -> R2S1RPHarmonic
generateR2S1RPHarmonicArray numTheta2 numScales rows cols angularFreq radialFreq maxScale alpha inverseFlag =
  let centerR = div rows 2
      centerC = div cols 2
      deltaTheta2 = 2 * pi / fromIntegral numTheta2
      deltaScale = maxScale / fromIntegral numScales
      sign =
        case inverseFlag of
          Harmonic -> (-1)
          ConjugateHarmonic -> 1
      arr =
        computeS $
        fromFunction
          (Z :. rows :. cols :. numTheta2 :. numScales)
          (\(Z :. i :. j :. k :. l) ->
             let x = fromIntegral $ i - centerR
                 y = fromIntegral $ j - centerC
                 theta1 = angleFunctionRad x y
                 theta2 = deltaTheta2 * fromIntegral k
                 r = sqrt $ x ^ (2 :: Int) + y ^ (2 :: Int)
             in if r <= (2 / pi * (max angularFreq radialFreq))
                  then 0
                  else ((r ** alpha) :+ 0) *
                       -- ((gaussian 4 r) :+ 0) *
                       exp
                         (0 :+
                          (sign *
                           ((angularFreq * (theta1 + theta2) / 2) +
                            (radialFreq * (theta1 - theta2) / 2)
                            -- (radialFreq *
                            --  (log r + log (fromIntegral (l + 1) * deltaScale)) /
                            --  2)
                           ))))
  in case inverseFlag of
       Harmonic -> R2S1RPHarmonic . HarmoicArray $ arr
       ConjugateHarmonic -> R2S1RPHarmonic . ConjugateHarmoicArray $ arr

{-# INLINE gaussian #-}

gaussian :: Double -> Double -> Double
gaussian sigma r =
  1 / ( (2 * pi) * sigma * sigma) * exp ((-1) * r * r / 2 / sigma / sigma)
