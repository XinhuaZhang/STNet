module App.STHarmonics.Generator where

import           Data.Array.Repa                as R
import           Data.Complex
import           Data.List                      as L
import           Src.Array.CoordinatesTransform
import           Src.Array.Transform
import           Src.Utils.Coordinates

type R2S1Array = Array U DIM3 (Complex Double)

-- (x,y,\theta) <-> (r,\theta_1, \theta_2) <-> (r, u, v)
-- (x,y) <-> (r,\theta_1)
-- u = (\theta_1 + \theta_2) / 2, v = (\theta_1 - \theta_2) / 2
-- \theta_1 = u + v, \theta_2 = u - v

data STHarmonicFlag
  = STHarmonic
  | STInverseHarmonic
  deriving (Read, Show)

-- Create Harmonics in (\theta_2, x, y) space,
-- assuming the center of polar coordinates is the center of the 2D array.
{-# INLINE generateSTHarmonicArray #-}

generateSTHarmonicArray
  :: Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> STHarmonicFlag
  -> R2S1Array
generateSTHarmonicArray numTheta2 rows cols angularFreq radialFreq alpha inverseFlag =
  let centerR = div rows 2
      centerC = div cols 2
      deltaTheta2 = 2 * pi / fromIntegral numTheta2
      sign =
        case inverseFlag of
          STHarmonic -> (-1)
          STInverseHarmonic -> 1
  in computeS $
     fromFunction
       (Z :. numTheta2 :. rows :. cols)
       (\(Z :. k :. i :. j) ->
          let x = fromIntegral $ i - centerR
              y = fromIntegral $ j - centerC
              theta1 = angleFunctionRad x y
              theta2 = deltaTheta2 * fromIntegral k
              r = sqrt $ x ^ (2 :: Int) + y ^ (2 :: Int)
          in if r <= (1 / pi * (fromIntegral . abs $ radialFreq))
               then 0
               else ((r ** alpha) :+ 0) *
                    exp
                      (0 :+
                       (sign *
                        ((fromIntegral angularFreq) * (theta1 + theta2) / 2 +
                         (fromIntegral radialFreq) * (log r)))))
