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
import           System.Directory
import           System.Environment
import           System.FilePath
import           System.Random

data STPoint
  = Source (Double, Int, Int)
  | Sink (Double, Int, Int)
  deriving (Read,Show)

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
  let (orientations:size:len:trails:threads:contrastN:_) =
        L.map (\x -> read x :: Int) . L.take 6 $ args
      (sigma:angularFreq':radialFreq:alpha:theta:_) =
        L.map (\x -> read x :: Double) . L.drop 6 $ args
      -- Generate initial distribution: arrRepa
      (xs:ys:[]) =
        L.group .
        L.sort . L.map (\x -> read x :: STPoint) . combineArgs . L.drop 11 $
        args
      arrSourceRepa =
        initialDist size [-angularFreq .. angularFreq] .
        L.map (\(Source (t, x, y)) -> (t / 360 * 2 * pi, x, y)) $
        xs
      arrSinkRepa =
        initialDist size [-angularFreq .. angularFreq] .
        L.map (\(Sink (t, x, y)) -> (t / 360 * 2 * pi, x, y)) $
        ys
      init = (0, 0, 0, 1)
      angularFreq = round angularFreq'
      freqs =
        generateHarmonicCoefficients
          (theta / 360 * 2 * pi)
          [-angularFreq .. angularFreq]
  arrs <-
    M.mapM
      (\f -> solveMonteCarloR2S1'' threads trails size 5 f sigma len init)
      [-angularFreq .. angularFreq]
  -- let arr =
  --       freqDomainR2S1 orientations . L.foldl1' (R.zipWith (+)) $
  --       L.zipWith (\x y -> R.map (* x) y) freqs arrs
  -- plotImageRepa "rotation.png" .
  --   Image 8 .
  --   computeS .
  --   R.extend (Z :. (1 :: Int) :. All :. All) . R.map magnitude . R.sumS $
  --   arr
  plan <- generateDFTPlan getEmptyPlan . L.head $ arrs
  -- Source
  sourceArrs <-
    L.foldl1' (R.zipWith (+)) . L.map delay <$>
    (MP.sequence $
     L.zipWith (\x y -> crosscorrelation plan x y) arrSourceRepa arrs)
  -- plotImageRepa "Source.png" .
  --   Image 8 .
  --   computeS .
  --   R.extend (Z :. (1 :: Int) :. All :. All) .
  --   R.map magnitude . R.sumS . freqDomainR2S1 orientations $
  --   sourceArrs
  -- Sink
  sinkArrs <-
    L.foldl1' (R.zipWith (+)) . L.map delay <$>
    (MP.sequence $
     L.zipWith (\x y -> crosscorrelation plan x y) arrSinkRepa arrs)
  -- plotImageRepa "Sink.png" .
  --   Image 8 .
  --   computeS .
  --   R.extend (Z :. (1 :: Int) :. All :. All) .
  --   R.map magnitude . R.sumS . freqDomainR2S1 orientations $
  --   sinkArrs
  -- Completion Field
  completionField <- completionR2S1 plan sourceArrs (timeReversalR2S1 sinkArrs)
  plotImageRepa "Completion.png" .
    Image 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.map magnitude . R.sumS . freqDomainR2S1 orientations $
    completionField
