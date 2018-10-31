import           App.STHarmonics.GeneralizedConvolution
import           App.STHarmonics.R2S1RPHarmonics
import           App.STHarmonics.Utils
import           Control.Monad                          as M
import           Data.Array                             as Arr
import           Data.Array.Repa                        as R
import           Data.Complex
import           Data.List                              as L
import           Data.Vector.Storable                   as VS
import           Data.Vector.Unboxed                    as VU
import           Src.Array.TwoHalfDArray
import           Src.FokkerPlanck.MonteCarlo
import           Src.Image.ImageIO
import           Src.Utils.Coordinates
import           Src.Utils.DFT
import           System.Directory
import           System.Environment
import           System.FilePath
import           System.Random

data STPoint
  = Source (Double, Double, Int, Int)
  | Sink (Double, Double, Int, Int)
  deriving (Read, Show)

instance Eq STPoint where
  (==) (Source _) (Source _) = True
  (==) (Sink _) (Sink _)     = True
  (==) (Source _) (Sink _)   = False
  (==) (Sink _) (Source _)   = False

instance Ord STPoint where
  compare (Source x) (Source y) = compare x y
  compare (Sink x) (Sink y)     = compare x y
  compare (Source _) (Sink _)   = LT
  compare (Sink _) (Source _)   = GT


{-# INLINE reduceContrast #-}

reduceContrast
  :: Int -> VU.Vector Double -> VU.Vector Double
reduceContrast idx vec =
  let x = L.head . L.drop idx . L.reverse . L.sort . VU.toList $ vec
  in VU.map
       (\y ->
          if y >= x
            then x
            else y)
       vec

{-# INLINE checkDeg #-}

checkDeg :: Double -> Double
checkDeg deg =
  if deg < 0
    then checkDeg (deg + 360)
    else if deg >= 360
           then checkDeg (deg - 360)
           else deg

{-# INLINE combineArgs #-}

combineArgs :: [String] -> [String]
combineArgs (x:y:[]) = [x L.++ " " L.++ y]
combineArgs (x:y:xs) = (x L.++ " " L.++ y) : combineArgs xs

main = do
  args <- getArgs
  print args
  let (orientations:scales:size:len:trails:threads:contrastN:_) =
        L.map (\x -> read x :: Int) . L.take 7 $ args
      (thetaSigma:scaleSigma:maxScale:angularFreq:radialFreq:alpha:_) =
        L.map (\x -> read x :: Double) . L.drop 7 $ args
      -- Generate initial distribution: arrRepa
      (xs:ys:[]) =
        L.group .
        L.sort . L.map (\x -> read x :: STPoint) . combineArgs . L.drop 13 $
        args
      arrSource =
        accumArray
          (+)
          (0 :: Double)
          ((0, 0, 0, 0), (size - 1, size - 1, orientations - 1, scales - 1)) .
        L.map (\idx -> (idx, 1 / (fromIntegral . L.length $ xs))) .
        L.map
          (\(Source (t, s, x, y)) ->
             ( x
             , y
             , floor ((checkDeg t) / 360 * fromIntegral orientations)
             , floor $ s / maxScale * fromIntegral scales)) $
        xs
      arrSourceRepa =
        fromListUnboxed (Z :. size :. size :. orientations :. scales) .
        L.map (:+ 0) . Arr.elems $
        arrSource
      r2s1rpHarmonics =
        L.map
          (\(af, rf) ->
             generateR2S1RPHarmonicArray
               orientations
               scales
               size
               size
               af
               rf
               maxScale
               alpha
               Harmonic)
          [(af, rf) | af <- [0 .. angularFreq], rf <- [0 .. radialFreq]]
      r2s1rpHarmonicsC =
        L.map
          (\(af, rf) ->
             generateR2S1RPHarmonicArray
               orientations
               scales
               size
               size
               af
               rf
               maxScale
               alpha
               ConjugateHarmonic)
          [(af, rf) | af <- [0 .. angularFreq], rf <- [0 .. radialFreq]]
      folderName = "GreensData"
  print xs
  arrG <-
    (computeS . R.map (:+ 0)) <$>
    readRepaArray
      (folderName </>
       ("Greens_" L.++ show size L.++ "_" L.++ show orientations L.++ "_" L.++
        show scales L.++
        "_" L.++
        show thetaSigma L.++
        "_" L.++
        show scaleSigma L.++
        "_" L.++
        show maxScale L.++
        ".dat")) :: IO (R.Array U DIM4 (Complex Double))
  let a = projectFilterR2S1RP arrG r2s1rpHarmonics
      b = projectFilterR2S1RP arrSourceRepa r2s1rpHarmonics
      c = L.zipWith (*) a b
      d = recoverFilterR2S1RP a r2s1rpHarmonicsC
  plotImageRepa "GeneralizedConvolutionTestR2S1RP/recon1.png" .
    Image 8 .
    fromUnboxed (Z :. (1 :: Int) :. size :. size) .
    R.toUnboxed . R.sumS . R.sumS . R.map magnitude $ d
