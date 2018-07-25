import           App.STHarmonics.WeightSharing
import           Control.Monad                 as M
import           Data.Array                    as Arr
import           Data.Array.Repa               as R
import           Data.Complex
import           Data.List                     as L
import           Data.Vector.Storable          as VS
import           Data.Vector.Unboxed           as VU
import           Src.Array.Transform           (rotate3D)
import           Src.FokkerPlanck.MonteCarlo
import           Src.Image.ImageIO
import           System.Environment
import           System.Random

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

main = do
  args <- getArgs
  let (orientations:size:len:trails:threads:contrastN:_) =
        L.map (\x -> read x :: Int) . L.take 6 $ args
      (sigma:_) = L.map (\x -> read x :: Double) . L.drop 6 $ args
      -- Generate initial distribution: arrRepa
      xs =
        L.map
          (\(t, x, y) -> (floor (t / 360 * fromIntegral orientations), x, y)) .
        L.map (\x -> read x :: (Double, Int, Int)) . L.drop 7 $
        args -- (theta,x,y), theta is in degree, i.e. (0,359)
      arr =
        accumArray
          (+)
          (0 :: Double)
          ((0, 0, 0), (orientations - 1, size - 1, size - 1)) .
        L.map (\idx -> (idx, 1 / (fromIntegral . L.length $ xs))) $
        xs
      arrRepa =
        fromListUnboxed (Z :. orientations :. size :. size) . Arr.elems $ arr
  -- Generate DFT plan
  lock <- getFFTWLock
  vecTemp1 <-
    VS.fromList <$> M.replicateM (orientations * size * size) randomIO :: IO (VS.Vector Double)
  vecTemp2 <-
    VS.fromList <$> M.replicateM (orientations * size * size) randomIO :: IO (VS.Vector Double)
  (dftPlan1, vecTemp3) <-
    dft1dGPlan
      lock
      getEmptyPlan
      [orientations, size, size]
      [1, 2]
      (VS.zipWith (:+) vecTemp1 vecTemp2)
  (dftPlan, _) <-
    idft1dGPlan lock dftPlan1 [orientations, size, size] [1, 2] vecTemp3
  -- Generate Green's function
  arrG <- solveMonteCarlo threads trails size orientations sigma len (0, 0, 0)
  let arrG' =
        fromUnboxed (Z :. (1 :: Int) :. size :. size) .
        reduceContrast contrastN . toUnboxed . R.sumS $
        arrG
  plotImageRepa "green.png" . Image 8 $ arrG'
  -- Source field
  ys <-
    shareWeightST dftPlan CrossCorrelationST (0, 1) arrRepa .
    computeS . rotate3D $
    arrG
  let arr' =
        fromUnboxed (Z :. (1 :: Int) :. size :. size) .
        reduceContrast contrastN . VU.concat $
        ys
  plotImageRepa "test.png" . Image 8 $ arr'
