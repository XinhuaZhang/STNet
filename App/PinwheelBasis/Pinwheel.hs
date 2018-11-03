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

data HarmonicArray 
  = HarmonicArray (Array U DIM3 (Complex Double))
  | ConjugateHarmonicArray (Array U DIM3 (Complex Double))

instance Show HarmonicArray where
  show (HarmonicArray arr)          = "HarmonicArray " L.++ (showShape . extent $ arr)
  show (ConjugateHarmonicArray arr) = "ConjugateHarmonicArray" L.++ (showShape . extent $ arr)


{-# INLINE generateR2S1HarmonicArray #-}

generateR2S1HarmonicArray ::
     Int
  -> Int
  -> [Double]
  -> [Double]
  -> Double
  -> HarmonicFlag
  -> HarmonicArray
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
                   else ((r :+ 0) **
                         (alpha :+
                          (sign * rf * 2 * pi /
                           (log . fromIntegral $ div rows 2)))) *
                        exp (0 :+ (sign * af * theta)))
   in case inverseFlag of
        Harmonic -> HarmonicArray $ arr
        ConjugateHarmonic -> ConjugateHarmonicArray $ arr


{-# INLINE impulsiveCoefficients  #-}
impulsiveCoefficients ::
     Double -> Double -> Double -> [Double] -> [Double] -> [Complex Double]
impulsiveCoefficients r theta maxR angularFreqs radialFreqs =
  L.map
    (\(af, rf) -> exp $ 0 :+ ((-1) * ((2 * pi / (log maxR)) * r * rf + theta * af)))
    [(af, rf) | af <- angularFreqs, rf <- radialFreqs]

{-# INLINE gabor #-}
gabor ::
     Int -> Int -> Double -> Double -> Double -> R.Array U DIM2 (Complex Double)
gabor rows cols freq theta sigma =
  let centerR = div rows 2
      centerC = div cols 2
   in computeS . fromFunction (Z :. rows :. cols) $ \(Z :. i' :. j') ->
        let i = i' - centerR
            j = j' - centerC
         in ((exp ((-1) * (fromIntegral $ (i) ^ 2 + j ^ 2) / 2 / sigma / sigma)) :+
             0) *
            (exp $
             0 :+
             (2 * pi * freq *
              (fromIntegral i * cos theta + fromIntegral j * sin theta)))
