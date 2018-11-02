module App.PinwheelBasis.Pinwheel where

import           Control.Parallel.Strategies
import           Data.Array.Repa                as R
import           Data.Complex
import           Data.List                      as L
import           Src.Utils.Coordinates


data HarmonicFlag
  = Harmonic
  | ConjugateHarmonic
  deriving (Read, Show)

data HarmoicArray 
  = HarmoicArray (Array U DIM3 (Complex Double))
  | ConjugateHarmoicArray (Array U DIM3 (Complex Double))

instance Show HarmoicArray where
  show (HarmoicArray arr)          = "HarmoicArray " L.++ (showShape . extent $ arr)
  show (ConjugateHarmoicArray arr) = "ConjugateHarmoicArray" L.++ (showShape . extent $ arr)


{-# INLINE generateR2S1HarmonicArray #-}

generateR2S1HarmonicArray ::
     Int
  -> Int
  -> [Double]
  -> [Double]
  -> Double
  -> HarmonicFlag
  -> HarmoicArray
generateR2S1HarmonicArray rows cols angularFreq radialFreq alpha inverseFlag =
  let centerR = div rows 2
      centerC = div cols 2
      numFreq = L.length angularFreq * L.length radialFreq
      sign =
        case inverseFlag of
          Harmonic -> (-1)
          ConjugateHarmonic -> 1
      freqArr =
        fromListUnboxed
          (Z :. numFreq)
          [(af, rf) | af <- angularFreq, rf <- radialFreq]
      initArr = fromFunction (Z :. rows :. cols :. numFreq) (const 0)
      arr =
        computeS $
        traverse2
          initArr
          freqArr
          const
          (\_ f1d (Z :. i :. j :. k) ->
             let x = fromIntegral $ i - centerR
                 y = fromIntegral $ j - centerC
                 theta = angleFunctionRad x y
                 r = sqrt $ x ** 2 + y ** 2
                 (af, rf) = f1d (Z :. k)
              in if r == 0
                   then 0
                   else ((r :+ 0) ** (alpha :+ (sign * af))) *
                        exp (0 :+ (sign * rf * theta)))
   in case inverseFlag of
        Harmonic -> HarmoicArray $ arr
        ConjugateHarmonic -> ConjugateHarmoicArray $ arr
