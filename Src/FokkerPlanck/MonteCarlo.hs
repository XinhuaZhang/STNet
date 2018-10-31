module Src.FokkerPlanck.MonteCarlo
  ( module Src.FokkerPlanck.Types
  , solveMonteCarloR2S1
  , solveMonteCarloR2S1RP
  ) where

import           Control.Arrow
import           Control.Monad                  as M
import           Control.Monad.Parallel         as MP
import           Data.Array                     as Arr
import           Data.Array.Repa                as R
import           Data.List                      as L
import           Data.Random.Normal
import           Data.Vector                    as V
import           Data.Vector.Unboxed            as VU
import           Src.Utils.Coordinates
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

thetaCheck :: Double -> Double
thetaCheck theta =
  if theta < 0
    then theta + 2 * pi
    else if theta >= 2 * pi
           then theta - 2 * pi
           else theta
           
{-# INLINE scalePlus #-}

scalePlus :: Double -> Double -> Double -> Double
scalePlus maxScale x y
  | z <= 0 = 0
  | z >= maxScale = maxScale
  | otherwise = z
  where
    z = x + y

{-# INLINE generatePath #-}

generatePath
  :: (RandomGen g)
  => g -> Double -> Double -> Double -> Int -> ParticleIndex -> (g, VU.Vector ParticleIndex)
generatePath randomGen thetaSigma scaleSigma maxScale numSteps init =
  second VU.fromList . go randomGen numSteps $ init
  where
    go rGen 0 _ = (rGen, [])
    go rGen n (x, y, theta, scale) =
      let (deltaTheta, newGen1) = normal' (0, thetaSigma) rGen
          (deltaScale, newGen2) =
            if scaleSigma == 0
              then (0, newGen1)
              else normal' (0, scaleSigma) newGen1
          newIndex =
            ( (x + scale * cos theta)
            , (y + scale * sin theta)
            , (theta `thetaPlus` deltaTheta)
            , (scalePlus maxScale scale deltaScale))
          (gen, xs) = go newGen2 (n - 1) newIndex
      in (gen, newIndex : xs)

{-# INLINE generatePathList #-}

generatePathList
  :: (RandomGen g)
  => Int
  -> Double
  -> Double
  -> Double
  -> Int
  -> ParticleIndex
  -> g
  -> [VU.Vector ParticleIndex]
generatePathList 0 _ _ _ _ _ _ = []
generatePathList m thetaSigma scaleSigma maxScale numSteps init randomGen =
  let (newGen, x) =
        generatePath randomGen thetaSigma scaleSigma maxScale numSteps init
      xs =
        generatePathList
          (m - 1)
          thetaSigma
          scaleSigma
          maxScale
          numSteps
          init
          newGen
  in x : xs

{-# INLINE countR2S1 #-}

countR2S1 :: Int
      -> Int
      -> [VU.Vector ParticleIndex]
      -> (Int, VU.Vector Int)
countR2S1 len numOrientations xs =
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
             (\(x, y, theta, _) ->
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

solveMonteCarloR2S1
  :: Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Int
  -> ParticleIndex
  -> IO (R.Array U DIM3 Double)
solveMonteCarloR2S1 numGen numTrails numPoints numOrientations thetaSigma numSteps (i,j,o,s) = do
  gens <- M.replicateM numGen newStdGen
  let xs =
        parMap
          rdeepseq
          (generatePathList (div numTrails numGen) thetaSigma 0 10 numSteps (i,j,o,1))
          gens
      (ys, zs) =
        L.unzip $ parMap rdeepseq (countR2S1 numPoints numOrientations) xs
      totalNum = fromIntegral $ L.sum ys
      totalNumVec = L.foldl1' (VU.zipWith (+)) zs
  return .
    fromUnboxed (Z :. numPoints :. numPoints :. numOrientations) .
    VU.map (\x -> fromIntegral x / totalNum) $
    totalNumVec
    


{-# INLINE countR2S1RP #-}

countR2S1RP :: Int
            -> Int
            -> Int
            -> Double
            -> [VU.Vector ParticleIndex]
            -> (Int, VU.Vector Int)
countR2S1RP len numOrientations numScales maxScale xs =
  let maxIdx = len - 1
      shift = maxIdx - (div len 2)
      deltaTheta = 2 * pi / (fromIntegral numOrientations)
      deltaScale = maxScale / (fromIntegral numScales)
      ys =
        L.map
          (VU.filter
             (\(x, y, _, scale) ->
                if (x <= maxIdx) &&
                   (y <= maxIdx) &&
                   (scale < numScales) && (x >= 0) && (y >= 0) && (scale >= 0)
                  then True
                  else False) .
           VU.map
             (\(x, y, theta, scale) ->
                ( (round x :: Int) + shift
                , (round y :: Int) + shift
                , (floor $ theta / deltaTheta :: Int)
                , (floor $ scale / deltaScale :: Int))) .
           VU.filter
             (\(_, _, _, s) ->
                if s == 0
                  then False
                  else True))
          xs
      numTrajectories = L.sum . L.map VU.length $ ys
      arr =
        accumArray
          (+)
          0
          ((0, 0, 0, 0), (maxIdx, maxIdx, numOrientations - 1, numScales - 1)) .
        L.concatMap (L.map (\idx -> (idx, 1)) . VU.toList) $
        ys
  in (numTrajectories, VU.fromList . elems $ arr)

-- total number of trails = numGen * numTrails

solveMonteCarloR2S1RP
  :: Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Int
  -> ParticleIndex
  -> IO (R.Array U DIM4 Double)
solveMonteCarloR2S1RP numGen numTrails numPoints numOrientations numScales thetaSigma scaleSigma maxScale numSteps init = do
  gens <- M.replicateM numGen newStdGen
  let xs =
        parMap
          rdeepseq
          (generatePathList
             (div numTrails numGen)
             thetaSigma
             scaleSigma
             maxScale
             numSteps
             init)
          gens
      (ys, zs) =
        L.unzip $
        parMap
          rdeepseq
          (countR2S1RP numPoints numOrientations numScales maxScale)
          xs
      totalNum = fromIntegral $ L.sum ys
      totalNumVec = L.foldl1' (VU.zipWith (+)) zs
  return .
    fromUnboxed (Z :. numPoints :. numPoints :. numOrientations :. numScales) .
    VU.map (\x -> fromIntegral x / totalNum) $
    totalNumVec
