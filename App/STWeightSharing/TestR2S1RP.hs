{-# LANGUAGE BangPatterns #-}
import           App.STHarmonics.Utils
import           App.STWeightSharing.Convolution
import           Control.Monad                   as M
import           Control.Monad.Parallel          as MP
import           Control.Parallel.Strategies
import           Data.Array                      as Arr
import           Data.Array.Repa                 as R
import           Data.Complex
import           Data.List                       as L
import           Data.Vector.Storable            as VS
import           Data.Vector.Unboxed             as VU
import           Src.FokkerPlanck.MonteCarlo
import           Src.Image.ImageIO
import           Src.Utils.Coordinates
import           Src.Utils.DFT
import           Src.Utils.Time
import           System.Directory
import           System.Environment
import           System.FilePath
import           System.Random

data STPoint
  = Source (Int, Int, Double, Double)
  | Sink (Int, Int, Double, Double)
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

{-# INLINE combineArgs #-}

combineArgs :: [String] -> [String]
combineArgs (x:y:[]) = [x L.++ " " L.++ y]
combineArgs (x:y:xs) = (x L.++ " " L.++ y) : combineArgs xs



main = do
  args <- getArgs
  print args
  let (orientations:scales:size:len:trails:threads:contrastN:angularFreq:radialFreq:_) =
        L.map (\x -> read x :: Int) . L.take 9 $ args
      (thetaSigma:scaleSigma:theta:scale:maxScale:tao:_) =
        L.map (\x -> read x :: Double) . L.drop 9 $ args
      -- Generate initial distribution: arrRepa
      (xs:ys:[]) =
        L.group .
        L.sort . L.map (\x -> read x :: STPoint) . combineArgs . L.drop 15 $
        args
      arrSourceRepa =
        initialDistR2S1RP
          size
          [-angularFreq .. angularFreq]
          [-radialFreq .. radialFreq]
          maxScale .
        L.map (\(Source (x, y, t, s)) -> (x, y, t / 360 * 2 * pi, s)) $
        xs
      arrSinkRepa =
        initialDistR2S1RP
          size
          [-angularFreq .. angularFreq]
          [-radialFreq .. radialFreq]
          maxScale .
        L.map (\(Sink (x, y, t, s)) -> (x, y, t / 360 * 2 * pi, s)) $
        ys
      init = (0, 0, 0, 1)
      freqs =
        generateHarmonicCoefficientsR2S1RP
          (theta / 360 * 2 * pi)
          [-angularFreq .. angularFreq]
          scale
          [-radialFreq .. radialFreq]
          maxScale
      freqs1 = [-5 .. 5]
  -- printCurrentTime
  -- !arrs <-
  --   solveMonteCarloR2S1RP''
  --     threads
  --     trails
  --     size
  --     freqs1
  --     [-angularFreq .. angularFreq]
  --     freqs1
  --     [-radialFreq .. radialFreq]
  --     thetaSigma
  --     scaleSigma
  --     maxScale
  --     tao
  --     len
  --     init  
  -- writeRepaArrays "greenR2S1RP.dat" arrs
  -- printCurrentTime
  printCurrentTime
  arrs <- readRepaArrays "greenR2S1RP.dat"
  printCurrentTime
  print . extent . L.head $ arrs
  let arr =
        freqDomainR2S1RP orientations scales maxScale .
        fromUnboxed (extent . L.head $ arrs) . L.foldl1' (VU.zipWith (+)) $
        L.zipWith (\x y -> VU.map (* x) . toUnboxed $ y) freqs arrs
  plotImageRepa "rotation.png" .
    Image 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) . R.map magnitude . R.sumS . R.sumS $
    arr
  printCurrentTime
  plan <- generateDFTPlanR2S1RP getEmptyPlan . L.head $ arrs
  -- Source
  sourceArrs <-
    L.foldl1' (R.zipWith (+)) . L.map delay <$>
    (MP.sequence $
     L.zipWith (\x y -> crosscorrelationR2S1RP plan x y) arrSourceRepa arrs)
  plotImageRepa "Source.png" .
    Image 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.map magnitude .
    R.sumS . R.sumS . freqDomainR2S1RP orientations scales maxScale $
    sourceArrs
  -- Sink
  sinkArrs <-
    L.foldl1' (R.zipWith (+)) . L.map delay <$>
    (MP.sequence $
     L.zipWith (\x y -> crosscorrelationR2S1RP plan x y) arrSinkRepa arrs)
  plotImageRepa "Sink.png" .
    Image 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.map magnitude .
    R.sumS . R.sumS . freqDomainR2S1RP orientations scales maxScale $
    sinkArrs
  -- Completion Field
  completionField <-
    completionR2S1RP plan sourceArrs (timeReversalR2S1RP sinkArrs)
  plotImageRepa "Completion.png" .
    Image 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.map magnitude .
    R.sumS . R.sumS . freqDomainR2S1RP orientations scales maxScale $
    completionField
