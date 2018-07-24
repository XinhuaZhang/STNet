module Src.Utils.DFT
  ( module Src.Utils.DFT.Plan
  , importFFTWWisdom
  , exportFFTWWisdom
  ) where

import           Control.Monad       (unless)
import           Src.Utils.DFT.Plan
import           Src.Utils.DFT.Base (exportWisdomString, importWisdomString)


importFFTWWisdom :: FilePath -> IO ()
importFFTWWisdom path = do
  str <- readFile path
  flag <- importWisdomString str
  unless flag (error $ "initializefftw: importWisdomString (" ++ str ++ ")")
  
exportFFTWWisdom :: FilePath -> IO ()
exportFFTWWisdom path = do
  wisdom <- exportWisdomString
  writeFile path wisdom
