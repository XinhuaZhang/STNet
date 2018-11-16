module App.STWeightSharing.Convolution where

import           Data.Complex
import           Data.List    as L

{-# INLINE generateHarmonicCoefficient' #-}
generateHarmonicCoefficient' :: Double -> Int -> Complex Double
generateHarmonicCoefficient' theta freq = exp $ 0 :+ (-1) * (fromIntegral freq) * theta

{-# INLINE generateHarmonicCoefficients #-}
generateHarmonicCoefficients :: Double -> [Int] -> [Complex Double]
generateHarmonicCoefficients theta = L.map (generateHarmonicCoefficient' theta)
