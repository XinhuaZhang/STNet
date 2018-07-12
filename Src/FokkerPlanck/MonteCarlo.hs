module Src.FokkerPlanck.MonteCarlo
  ( module Src.FokkerPlanck.Types
  , solveMonteCarlo
  ) where

import           Control.Arrow
import           Control.Monad          as M
import           Control.Monad.Parallel as MP
import           Data.Array             as Arr
import           Data.Array.Repa        as R
import           Data.List              as L
import           Data.Random.Normal
import           Data.Vector            as V
import           Data.Vector.Unboxed    as VU
import           Src.FokkerPlanck.Types
import           Src.Utils.Parallel
import           System.Random

{-# INLINE thetaPlus #-}

thetaPlus :: Double -> Double -> Double
thetaPlus x y
  | z < 0 = z + a
  | z >= a = z - a
  | otherwise = z
  where
    z = x + y
    a = 2 * pi
    
{-# INLINE thetaCheck #-}

thetaCheck :: Double -> Bool
thetaCheck theta
  | theta < (2 * pi) && theta >= 0 = True
  | otherwise = False

{-# INLINE generatePath #-}

generatePath
  :: (RandomGen g)
  => g -> Double -> Int -> ParticleIndex -> (g, VU.Vector ParticleIndex)
generatePath randomGen sigma n init = second VU.fromList . go randomGen n $ init
  where
    go rGen 0 _ = (rGen, [])
    go rGen n (x, y, theta) =
      let (deltaTheta, newGen) = normal' (0, sigma) rGen
          newIndex =
            ((x + cos theta), (y + sin theta), (theta `thetaPlus` deltaTheta))
          (gen, xs) = go newGen (n - 1) newIndex
      in (gen, newIndex : xs)

{-# INLINE generatePathList #-}

generatePathList
  :: (RandomGen g)
  => Int ->  Double -> Int -> ParticleIndex -> g -> [VU.Vector ParticleIndex]
generatePathList 0 _ _ _ _ = []
generatePathList m sigma n init randomGen =
  let (newGen, x) = generatePath randomGen sigma n init
      xs = generatePathList (m - 1) sigma n init newGen
  in x : xs

{-# INLINE count #-}

count :: Int
      -> Int
      -> [VU.Vector ParticleIndex]
      -> (Int, VU.Vector Int)
count len numOrientations xs =
  let maxIdx = len - 1
      shift = maxIdx - (div len 2)
      deltaTheta = 2 * pi / (fromIntegral numOrientations)
      ys =
        L.map
          (VU.filter
             (\(x, y, theta) ->
                if (x <= maxIdx) && (y <= maxIdx) && (x >= 0) && (y >= 0)
                  then True
                  else False) .
           VU.map
             (\(x, y, theta) ->
                ( (round x :: Int) + shift
                , (round y :: Int) + shift
                , (floor $ theta / deltaTheta :: Int))))
          xs
      numTrajectories = L.sum . L.map VU.length $ ys
      arr =
        accumArray (+) 0 ((0, 0, 0), (maxIdx, maxIdx, numOrientations - 1)) .
        L.concatMap (L.map (\idx -> (idx, 1)) . VU.toList) $
        ys
  in (numTrajectories, VU.fromList . elems $ arr)

-- total number of trails = numGen * numTrails

solveMonteCarlo
  :: Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Int
  -> ParticleIndex
  -> IO (R.Array U DIM3 Double)
solveMonteCarlo numGen numTrails numPoints numOrientations sigma len init@(_, _, theta) = do
  unless
    (thetaCheck theta)
    (error $ "solveMonteCarlo: theta out of boundary.\n" L.++ show theta)
  gens <- M.replicateM numGen newStdGen
  let xs = parMap rdeepseq (generatePathList numTrails sigma len init) gens
      (ys, zs) = L.unzip $ parMap rdeepseq (count numPoints numOrientations) xs
      totalNum = fromIntegral $ L.sum ys
      totalNumVec = L.foldl1' (VU.zipWith (+)) zs
  return .
    fromUnboxed (Z :. numPoints :. numPoints :. numOrientations) .
    VU.map (\x -> fromIntegral x / totalNum) $
    totalNumVec
