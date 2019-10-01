## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(amber)

## ---- eval = FALSE-------------------------------------------------------
#  library(amber)
#  library(doParallel)
#  library(ncdf4)
#  library(raster)
#  
#  long.name <- 'Gross primary productivity'
#  nc.mod <- system.file('extdata/modelRegular', 'gpp_monthly.nc', package = 'amber')
#  nc.ref <- system.file('extdata/referenceRegular', 'gpp_GBAF_128x64.nc', package = 'amber')
#  mod.id <- 'CLASSIC' # define a model experiment ID
#  ref.id <- 'GBAF' # give reference dataset a name
#  unit.conv.mod <- 86400*1000 # optional unit conversion for model data
#  unit.conv.ref <- 86400*1000 # optional unit conversion for reference data
#  variable.unit <- 'gC m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
#  
#  plot.me <- scores.grid.time(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#  unit.conv.ref, variable.unit)
#  # Create a plot in PDF format
#  plotGrid(long.name, plot.me)

## ---- eval = FALSE-------------------------------------------------------
#  library(amber)
#  library(doParallel)
#  library(ncdf4)
#  library(raster)
#  
#  long.name <- 'Gross primary productivity'
#  nc.mod <- system.file('extdata/modelRegular', 'gpp_monthly.nc', package = 'amber')
#  ref.csv <- system.file('extdata/referenceRegular', 'gpp_monthly_fluxnet.csv', package = 'amber')
#  mod.id <- 'CLASSIC' # define a model experiment ID
#  ref.id <- 'FLUXNET' # give reference dataset a name
#  unit.conv.mod <- 86400*1000 # optional unit conversion for model data
#  unit.conv.ref <- 1 # optional unit conversion for reference data
#  variable.unit <- 'gC m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
#  
#  scores.fluxnet.csv(long.name, nc.mod, ref.csv, mod.id, ref.id,
#  unit.conv.mod, unit.conv.ref, variable.unit)

## ---- eval = FALSE-------------------------------------------------------
#  library(amber)
#  library(doParallel)
#  library(ncdf4)
#  library(raster)
#  library(latex2exp)
#  
#  long.name <- 'Streamflow'
#  nc.mod <- system.file('extdata/modelRegular', 'mrro_monthly.nc', package = 'amber')
#  nc.ref <- system.file('extdata/referenceRegular', 'runoff.nc', package = 'amber')
#  nc.basins <- system.file('extdata/referenceRegular', 'basins.nc', package = 'amber')
#  mod.id <- 'CLASSIC' # model name
#  ref.id <- 'GRDC' # give reference dataset a name
#  unit.conv.mod <- 86400 # optional unit conversion for model data
#  unit.conv.ref <- 86400 # optional unit conversion for reference data
#  variable.unit <- 'kg m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
#  score.weights <- c(1,2,1,1,1) # define score weights
#  
#  scores.runoff(long.name, nc.mod, nc.ref, nc.basins, mod.id, ref.id, unit.conv.mod,
#  unit.conv.ref, variable.unit, score.weights)

