{-# LANGUAGE BangPatterns #-}
module Src.FokkerPlanck.MonteCarlo
  ( module Src.FokkerPlanck.Types
  , solveMonteCarloR2S1
  -- , solveMonteCarloR2S1'
  -- , solveMonteCarloR2S1''
  -- , solveMonteCarloR2S1RP
  -- , solveMonteCarloR2S1RP''
  ) where

import           Control.Arrow
import           Control.Monad          as M
import           Control.Monad.Parallel as MP
import           Data.Array             as Arr
import           Data.Array.Repa        as R
import           Data.Complex
import           Data.List              as L
import           Data.Random.Normal
import           Data.Vector            as V
import           Data.Vector.Unboxed    as VU
import           Data.Vector    as V
import           Src.Array.UnboxedArray as UA
import           Src.FokkerPlanck.Types
import           Src.Utils.Coordinates
import           Src.Utils.Parallel
import           System.Random          as Random

{-# INLINE thetaPlus #-}

thetaPlus :: Double -> Double -> Double
thetaPlus x y
  | z < 0 =  z + a
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

scalePlus :: Double -> Double -> Double -> (Double, Bool)
scalePlus maxScale x y
  | z < 0 = (-z, True)
  | z >= maxScale = (maxScale, False)
  | otherwise = (z, False)
  where
    z = x + y

-- {-# INLINE generatePath #-}
-- generatePath ::
--      (RandomGen g)
--   => g
--   -> Double
--   -> Double
--   -> Double
--   -> Double
--   -> Int
--   -> ParticleIndex
--   -> (g, VU.Vector ParticleIndex)
-- generatePath randomGen thetaSigma scaleSigma maxScale tao numSteps init =
--   second VU.fromList . go randomGen numSteps $ init
--   where
--     go rGen 0 _ = (rGen, [])
--     go rGen n (!x, !y, !theta, !scale) =
--       let (deltaTheta, newGen1) = normal' (0, thetaSigma) rGen
--           (deltaScale, newGen2) =
--             if scaleSigma == 0
--               then (0, newGen1)
--               else normal' (0, scaleSigma) newGen1
--           (newScale, flag) = scalePlus maxScale scale deltaScale
--           newIndex =
--             ( (x + scale * cos theta)
--             , (y + scale * sin theta)
--             , (if flag
--                  then (theta `thetaPlus` (deltaTheta + pi))
--                  else (theta `thetaPlus` deltaTheta))
--             , newScale)
--           (t, newGen3) = random newGen2
--           (gen, ys) =
--             if t < (1 - exp ((-1) / tao))
--               then go newGen3 0 newIndex
--               else go newGen3 (n - 1) newIndex
--        in (gen, newIndex : ys)
       

{-# INLINE generatePath #-}
generatePath ::
     (RandomGen g)
  => g
  -> Double
  -> Double
  -> Double
  -> Double
  -> Int
  -> ParticleIndex
  -> V.Vector (g,ParticleIndex)
generatePath randomGen thetaSigma scaleSigma maxScale tao numSteps init =
  V.unfoldrN
    numSteps
    (\(!rGen, (!x, !y, !theta, !scale)) ->
       let (deltaTheta, newGen1) = normal' (0, thetaSigma) rGen
           (deltaScale, newGen2) =
             if scaleSigma == 0
               then (0, newGen1)
               else normal' (0, scaleSigma) newGen1
           (newScale, flag) = scalePlus maxScale scale deltaScale
           newIndex =
             ( (x + scale * cos theta)
             , (y + scale * sin theta)
             , (if flag
                  then (theta `thetaPlus` (deltaTheta + pi))
                  else (theta `thetaPlus` deltaTheta))
             , newScale)
           (t, !gen) = random newGen2
        in if t < (1 - exp ((-1) / tao))
             then Nothing
             else Just ((gen, newIndex), (gen, newIndex)))
    (randomGen, init)


-- {-# INLINE generatePathList #-}
-- generatePathList ::
--      (RandomGen g)
--   => Int
--   -> Double
--   -> Double
--   -> Double
--   -> Double
--   -> Int
--   -> ParticleIndex
--   -> g
--   -> [VU.Vector ParticleIndex]
-- generatePathList 0 _ _ _ _ _ _ _ = []
-- generatePathList m thetaSigma scaleSigma maxScale tao numSteps init randomGen =
--   let (newGen, x) =
--         generatePath randomGen thetaSigma scaleSigma maxScale tao numSteps init
--       xs =
--         generatePathList
--           (m - 1)
--           thetaSigma
--           scaleSigma
--           maxScale
--           tao
--           numSteps
--           init
--           newGen
--    in x : xs

{-# INLINE generatePathList #-}
generatePathList ::
     (RandomGen g)
  => Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Int
  -> ParticleIndex
  -> g
  -> [VU.Vector ParticleIndex]
generatePathList m thetaSigma scaleSigma maxScale tao numSteps init randomGen =
  V.toList $
  V.unfoldrN
    m
    (\gen ->
       let (newGen, x) =
             V.unzip $
             generatePath gen thetaSigma scaleSigma maxScale tao numSteps init
        in Just (V.convert x, V.last newGen))
    randomGen

{-# INLINE countR2S1 #-}
countR2S1 :: Int -> Int -> [VU.Vector ParticleIndex] -> (Int, VU.Vector Int)
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
                , (floor $ (thetaCheck theta) / deltaTheta :: Int))))
          xs
      numTrajectories = L.sum . L.map VU.length $ ys
      arr =
        UA.accumulate (+) 0 ((0, 0, 0), (maxIdx, maxIdx, numOrientations - 1)) .
        VU.concat . L.map (VU.map (\idx -> (idx, 1))) $
        ys
   in (numTrajectories, toUnboxedVector arr)

-- total number of trails = numGen * numTrails

{-# INLINE solveMonteCarloR2S1 #-}
solveMonteCarloR2S1 ::
     Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Int
  -> ParticleIndex
  -> IO (R.Array U DIM3 Double)
solveMonteCarloR2S1 numGen numTrails numPoints numOrientations thetaSigma tao numSteps (i, j, o, s) = do
  gens <- M.replicateM numGen newStdGen
  let xs =
        parMap
          rdeepseq
          (generatePathList
             (div numTrails numGen)
             thetaSigma
             0
             10
             tao
             numSteps
             (i, j, o, 1))
          gens
      (ys, zs) =
        L.unzip $ parMap rdeepseq (countR2S1 numPoints numOrientations) xs
      totalNum = fromIntegral $ L.sum ys
      totalNumVec = L.foldl1' (VU.zipWith (+)) zs
  print totalNum
  print $ L.length zs
  return .
    fromUnboxed (Z :. numPoints :. numPoints :. numOrientations) .
    VU.map (\x -> fromIntegral x / totalNum) $
    totalNumVec



-- {-# INLINE countR2S1RP #-}
-- countR2S1RP :: Int
--             -> Int
--             -> Int
--             -> Double
--             -> [VU.Vector ParticleIndex]
--             -> (Int, VU.Vector Int)
-- countR2S1RP len numOrientations numScales maxScale xs =
--   let maxIdx = len - 1
--       shift = maxIdx - (div len 2)
--       deltaTheta = 2 * pi / (fromIntegral numOrientations)
--       deltaScale = maxScale / (fromIntegral numScales)
--       ys =
--         L.map
--           (VU.filter
--              (\(x, y, _, scale) ->
--                 if (x <= maxIdx) &&
--                    (y <= maxIdx) &&
--                    (scale < numScales) && (x >= 0) && (y >= 0) && (scale >= 0)
--                   then True
--                   else False) .
--            VU.map
--              (\(x, y, theta, scale) ->
--                 ( (round x :: Int) + shift
--                 , (round y :: Int) + shift
--                 , (floor $ (thetaCheck theta) / deltaTheta :: Int)
--                 , (floor $ scale / deltaScale :: Int))) .
--            VU.filter
--              (\(_, _, _, s) ->
--                 if s == 0
--                   then False
--                   else True))
--           xs
--       numTrajectories = L.sum . L.map VU.length $ ys
--       arr =
--         UA.accumulate
--           (+)
--           0
--           ((0, 0, 0, 0), (maxIdx, maxIdx, numOrientations - 1, numScales - 1)) .
--         VU.concat . L.map (VU.map (\idx -> (idx, 1))) $
--         ys
--    in (numTrajectories, toUnboxedVector arr)

-- -- total number of trails = numGen * numTrails
-- {-# INLINE solveMonteCarloR2S1RP #-}
-- solveMonteCarloR2S1RP ::
--      Int
--   -> Int
--   -> Int
--   -> Int
--   -> Int
--   -> Double
--   -> Double
--   -> Double
--   -> Double
--   -> Int
--   -> ParticleIndex
--   -> IO (R.Array U DIM4 Double)
-- solveMonteCarloR2S1RP numGen numTrails numPoints numOrientations numScales thetaSigma scaleSigma maxScale tao numSteps init = do
--   gens <- M.replicateM numGen newStdGen
--   let xs =
--         parMap
--           rdeepseq
--           (generatePathList
--              (div numTrails numGen)
--              thetaSigma
--              scaleSigma
--              maxScale
--              tao
--              numSteps
--              init)
--           gens
--       (ys, zs) =
--         L.unzip $
--         parMap
--           rdeepseq
--           (countR2S1RP numPoints numOrientations numScales maxScale)
--           xs
--       totalNum = fromIntegral $ L.sum ys
--       totalNumVec = L.foldl1' (VU.zipWith (+)) zs
--   return .
--     fromUnboxed (Z :. numPoints :. numPoints :. numOrientations :. numScales) .
--     VU.map (\x -> fromIntegral x / totalNum) $
--     totalNumVec


-- {-# INLINE generatePathList' #-}
-- generatePathList' ::
--      (RandomGen g)
--   => Int
--   -> Double
--   -> Double
--   -> Double
--   -> Double
--   -> Int
--   -> [ParticleIndex]
--   -> g
--   -> [VU.Vector ParticleIndex]
-- generatePathList' 0 _ _ _ _ _ _ _ = []
-- generatePathList' m thetaSigma scaleSigma maxScale tao numSteps (init:inits) randomGen =
--   let (newGen, x) =
--         generatePath randomGen thetaSigma scaleSigma maxScale tao numSteps init
--       xs =
--         generatePathList'
--           (m - 1)
--           thetaSigma
--           scaleSigma
--           maxScale
--           tao
--           numSteps
--           inits
--           newGen
--    in x : xs

-- {-# INLINE countR2S1' #-}

-- countR2S1' ::
--      Int
--   -> Int
--   -> Int
--   -> [Double]
--   -> [VU.Vector ParticleIndex]
--   -> (Int, VU.Vector (Complex Double))
-- countR2S1' len numOrientations freq ts xs =
--   let maxIdx = len - 1
--       shift = maxIdx - (div len 2)
--       deltaTheta = 2 * pi / (fromIntegral numOrientations)
--       ys =
--         L.map
--           (VU.filter
--              (\(x, y, theta) ->
--                 if (x <= maxIdx) && (y <= maxIdx) && (x >= 0) && (y >= 0)
--                   then True
--                   else False) .
--            VU.map
--              (\(x, y, theta, _) ->
--                 ( (round x :: Int) + shift
--                 , (round y :: Int) + shift
--                 , (floor $ theta / deltaTheta :: Int))))
--           xs
--       numTrajectories = L.sum . L.map VU.length $ ys
--       arr =
--         accumArray (+) 0 ((0, 0, 0), (maxIdx, maxIdx, numOrientations - 1)) .
--         L.concat $
--         L.zipWith
--           (\vec t ->
--              L.map (\idx -> (idx, exp $ 0 :+ fromIntegral freq * t)) . VU.toList $
--              vec)
--           ys
--           ts
--    in (numTrajectories, VU.fromList . elems $ arr)

-- -- total number of trails = numGen * numTrails

-- solveMonteCarloR2S1' ::
--      Int
--   -> Int
--   -> Int
--   -> Int
--   -> Int
--   -> Double
--   -> Int
--   -> ParticleIndex
--   -> IO (R.Array U DIM3 (Complex Double))
-- solveMonteCarloR2S1' numGen numTrails numPoints numOrientations freq thetaSigma numSteps (i, j, o, s) = do
--   gens <- M.replicateM numGen newStdGen
--   oris <-
--     (L.map (L.map (* (2 * pi)))) <$>
--     M.replicateM numGen (M.replicateM (div numTrails numGen) randomIO) :: IO [[Double]]
--   let xs =
--         withStrategy (parList rdeepseq) $
--         L.zipWith
--           (\gen ori ->
--              let idx = L.map (\t -> (i, j, t, 1)) ori
--              in generatePathList'
--                   (div numTrails numGen)
--                   thetaSigma
--                   0
--                   10
--                   numSteps
--                   idx
--                   gen)
--           gens
--           oris
--       (ys, zs) =
--         L.unzip . withStrategy (parList rdeepseq) $
--         L.zipWith
--           (\x ts -> countR2S1' numPoints numOrientations freq ts x)
--           xs
--           oris
--       totalNum = fromIntegral $ L.sum ys
--       totalNumVec = L.foldl1' (VU.zipWith (+)) zs
--   return .
--     fromUnboxed (Z :. numPoints :. numPoints :. numOrientations) .
--     VU.map (/ totalNum) $
--     totalNumVec



-- {-# INLINE countR2S1'' #-}

-- countR2S1'' ::
--      Int
--   -> Int
--   -> Int
--   -> [Double]
--   -> [VU.Vector ParticleIndex]
--   -> (Int, VU.Vector (Complex Double))
-- countR2S1'' len freq1 freq ts xs =
--   let maxIdx = len - 1
--       shift = maxIdx - (div len 2)
--       ys =
--         L.map
--           (VU.toList .
--            VU.filter
--              (\((x, y), _) ->
--                 if (x <= maxIdx) && (y <= maxIdx) && (x >= 0) && (y >= 0)
--                   then True
--                   else False) .
--            VU.map
--              (\(x, y, theta, _) ->
--                 (((round x :: Int) + shift, (round y :: Int) + shift), theta))) $
--         xs
--       numTrajectories = L.sum . L.map L.length $ ys
--       arr =
--         accumArray (+) 0 ((0, 0, -freq1), (maxIdx, maxIdx, freq1)) .
--         L.concatMap
--           (\n ->
--              L.concat $
--              L.zipWith
--                (\x t ->
--                   L.map
--                     (\((a, b), t1) ->
--                        ( (a, b, n)
--                        , exp $
--                          0 :+ (fromIntegral freq * t + fromIntegral n * t1)))
--                     x)
--                ys
--                ts) $
--         [-freq1 .. freq1]
--    in (numTrajectories, VU.fromList . elems $ arr)


-- solveMonteCarloR2S1'' ::
--      Int
--   -> Int
--   -> Int
--   -> Int
--   -> Int
--   -> Double
--   -> Int
--   -> ParticleIndex
--   -> IO (R.Array U DIM3 (Complex Double))
-- solveMonteCarloR2S1'' numGen numTrails numPoints freq1 freq thetaSigma numSteps (i, j, o, s) = do
--   gens <- M.replicateM numGen newStdGen
--   oris <-
--     (L.map (L.map (* (2 * pi)))) <$>
--     M.replicateM numGen (M.replicateM (div numTrails numGen) randomIO) :: IO [[Double]]
--   let xs =
--         withStrategy (parList rdeepseq) $
--         L.zipWith
--           (\gen ori ->
--              let idx = L.map (\t -> (i, j, t, 1)) ori
--               in generatePathList'
--                    (div numTrails numGen)
--                    thetaSigma
--                    0
--                    10
--                    numSteps
--                    idx
--                    gen)
--           gens
--           oris
--       (ys, zs) =
--         L.unzip . withStrategy (parList rdeepseq) $
--         L.zipWith (\x ts -> countR2S1'' numPoints freq1 freq ts x) xs oris
--       totalNum = fromIntegral $ L.sum ys
--       totalNumVec = L.foldl1' (VU.zipWith (+)) zs
--   return .
--     fromUnboxed (Z :. numPoints :. numPoints :. (2 * freq1 + 1)) .
--     VU.map (/ totalNum) $
--     totalNumVec



-- {-# INLINE countR2S1RP'' #-}

-- countR2S1RP'' ::
--      Int
--   -> [Int]
--   -> Int
--   -> [Int]
--   -> Int
--   -> Double
--   -> [Double]
--   -> [Double]
--   -> [VU.Vector ((Int, Int), (Double, Double))]
--   -> (Int, VU.Vector (Complex Double))
-- countR2S1RP'' len angularFreqs angularFreq1 radialFreqs radialFreq1 maxScale ts ss xs =
--   let maxIdx = len - 1
--       shift = maxIdx - (div len 2)
--       -- ys =
--       --   L.map
--       --     (VU.filter
--       --        (\((x, y), (_, s)) ->
--       --           if (x <= maxIdx) &&
--       --              (y <= maxIdx) && (x >= 0) && (y >= 0) && (s > 0)
--       --             then True
--       --             else False) .
--       --      VU.map
--       --        (\(x, y, theta, scale) ->
--       --           ( ((round x :: Int) + shift, (round y :: Int) + shift)
--       --           , (theta, scale)))) $
--       --   xs
--       numTrajectories = L.sum . L.map VU.length $ xs
--       arr =
--         UA.accumulate
--           (+)
--           0
--           ( (0, 0, L.head angularFreqs, L.head radialFreqs)
--           , (maxIdx, maxIdx, L.last angularFreqs, L.last radialFreqs)) .
--         VU.concat .
--         L.map
--           (\(m, n) ->
--              VU.concat $
--              L.zipWith3
--                (\x t s ->
--                   VU.map
--                     (\((a, b), (t1, s1)) ->
--                        ( (a, b, m, n)
--                        , exp $
--                          0 :+
--                          (fromIntegral angularFreq1 * t + fromIntegral m * t1 +
--                           (fromIntegral radialFreq1 * s + fromIntegral n * s1) *
--                           2 *
--                           pi /
--                           maxScale)))
--                     x)
--                xs
--                ts
--                ss) $
--         [(m, n) | m <- angularFreqs, n <- radialFreqs]
--    in (numTrajectories, toUnboxedVector arr)


-- -- {-# INLINE generateRandomNumber #-}
-- -- generateRandomNumber :: (RandomGen g) => Int -> g -> [Double] -> (g, [Double])
-- -- generateRandomNumber 0 !gen !xs = (gen, xs)
-- -- generateRandomNumber n !gen !xs =
-- --   let (!x, !newGen) = random gen
-- --    in generateRandomNumber (n - 1) newGen (x : xs)

-- {-# INLINE generateRandomNumber #-}
-- generateRandomNumber :: (RandomGen g) => Int -> g -> [Double]
-- generateRandomNumber n = L.take n . randoms


-- -- {-# INLINE solveMonteCarloR2S1RP'' #-}
-- -- solveMonteCarloR2S1RP'' ::
-- --      (RandomGen g)
-- --   => (g,g,g)
-- --   -> Int
-- --   -> Int
-- --   -> [Int]
-- --   -> Int
-- --   -> [Int]
-- --   -> Int
-- --   -> Double
-- --   -> Double
-- --   -> Double
-- --   -> Int
-- --   -> ParticleIndex
-- --   -> R.Array U DIM4 (Complex Double)
-- -- solveMonteCarloR2S1RP'' (gen, gen1, gen2) numTrails numPoints angularFreqs angularFreq1 radialFreqs radialFreq1 thetaSigma scaleSigma maxScale numSteps (i, j, _, _) =
-- --   let oris = (L.map (* (2 * pi))) $ generateRandomNumber numTrails gen1
-- --       scales = (L.map (* maxScale)) $ generateRandomNumber numTrails gen2
-- --       xs =
-- --         let idx = L.zipWith (\t s -> (i, j, t, s)) oris scales
-- --          in generatePathList'
-- --               numTrails
-- --               thetaSigma
-- --               scaleSigma
-- --               maxScale
-- --               numSteps
-- --               idx
-- --               gen
-- --       (totalNum, totalNumVec) =
-- --         countR2S1RP''
-- --           numPoints
-- --           angularFreqs
-- --           angularFreq1
-- --           radialFreqs
-- --           radialFreq1
-- --           maxScale
-- --           oris
-- --           scales
-- --           xs
-- --       result =
-- --         fromUnboxed
-- --           (Z :. numPoints :. numPoints :. (L.length angularFreqs) :.
-- --            (L.length radialFreqs)) .
-- --         VU.map (/ fromIntegral totalNum) $
-- --         totalNumVec
-- --    in deepSeqArray result result


-- {-# INLINE solveMonteCarloR2S1RP'' #-}
-- solveMonteCarloR2S1RP'' ::
--      Int
--   -> Int
--   -> Int
--   -> [Int]
--   -> [Int]
--   -> [Int]
--   -> [Int]
--   -> Double
--   -> Double
--   -> Double
--   -> Int
--   -> ParticleIndex
--   -> IO [R.Array U DIM4 (Complex Double)]
-- solveMonteCarloR2S1RP'' numGen numTrails numPoints angularFreqs angularFreqs1 radialFreqs radialFreqs1 thetaSigma scaleSigma maxScale numSteps (i, j, _, _) = do
--   gens <- M.replicateM numGen newStdGen
--   oris <-
--     (L.map (L.map (* (2 * pi)))) <$>
--     M.replicateM numGen (M.replicateM (div numTrails numGen) randomIO) :: IO [[Double]]
--   scales <-
--     (L.map (L.map (* maxScale))) <$>
--     M.replicateM numGen (M.replicateM (div numTrails numGen) randomIO) :: IO [[Double]]
--   let xs =
--         L.concat .
--         parMap
--           rdeepseq
--           (\(gen, ori, scale) ->
--              let idx = L.zipWith (\t s -> (i, j, t, s)) ori scale
--                  maxIdx = numPoints - 1
--                  shift = maxIdx - (div numPoints 2)
--                  as =
--                    generatePathList'
--                      (div numTrails numGen)
--                      thetaSigma
--                      scaleSigma
--                      maxScale
--                      numSteps
--                      idx
--                      gen
--               in L.map
--                    (VU.filter
--                       (\((x, y), (_, s)) ->
--                          if (x <= maxIdx) &&
--                             (y <= maxIdx) && (x >= 0) && (y >= 0) && (s > 0)
--                            then True
--                            else False) .
--                     VU.map
--                       (\(x, y, theta, scale) ->
--                          ( ((round x :: Int) + shift, (round y :: Int) + shift)
--                          , (theta, scale)))) $
--                  as) $
--         L.zip3 gens oris scales
--       ys =
--         parMap
--           rdeepseq
--           (\(af, rf) ->
--              (\(totalNum, totalNumVec) ->
--                 VU.map (/ fromIntegral totalNum) totalNumVec) $
--              countR2S1RP''
--                numPoints
--                angularFreqs
--                af
--                radialFreqs
--                rf
--                maxScale
--                (L.concat oris)
--                (L.concat scales)
--                xs)
--           [(af, rf) | af <- angularFreqs1, rf <- radialFreqs1]
--   return .
--     L.map
--       (fromUnboxed
--          (Z :. numPoints :. numPoints :. (L.length angularFreqs) :.
--           (L.length radialFreqs))) $
--     ys
