library(testthat)
library(amber)

Sys.setenv("R_TESTS" = "")

rm(list = ls())
library(amber)
library(doParallel)
library(ncdf4)
library(raster)
library(xtable)

#-------------------------------------------------------------------------------

# scores.grid.time

#-------------------------------------------------------------------------------
test_that(
  "Test scores.grid.time",
  {
    skip_on_cran()
    long.name <- 'Gross primary productivity'
    nc.mod <- system.file('extdata/modelRegular', 'gpp_monthly.nc', package = 'amber')
    nc.ref <- system.file('extdata/referenceRegular', 'gpp_GBAF_128x64.nc', package = 'amber')
    mod.id <- 'CLASSIC' # define a model experiment ID
    ref.id <- 'GBAF' # give reference dataset a name
    unit.conv.mod <- 86400*1000 # optional unit conversion for model data
    unit.conv.ref <- 86400*1000 # optional unit conversion for reference data
    variable.unit <- 'gC m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
    result <- amber::scores.grid.time(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
                                      unit.conv.ref, variable.unit)
    # Sample result
    myMean <- mean(result[[1]])
    sample <- myMean[155:160]
    computedValue <- sample
    documentedValue <- c(1.746210, 1.743733, 1.492706, 1.904318, 1.927370, 1.886606)

    expect_equal(computedValue, documentedValue, tolerance = 0.001)
  }
)

#-------------------------------------------------------------------------------

# scores.grid.notime

#-------------------------------------------------------------------------------
test_that(
  "Test scores.grid.notime",
  {
    skip_on_cran()
    long.name <- 'Soil Carbon'
    nc.mod <- system.file('extdata/modelRegular', 'cSoil_monthly.nc', package = 'amber')
    nc.ref <- system.file('extdata/referenceRegular', 'soilc_HWSD_128x64.nc', package = 'amber')
    mod.id <- 'CLASSIC' # define a model experiment ID
    ref.id <- 'HWSD' # give reference dataset a name
    unit.conv.mod <- 1 # optional unit conversion for model data
    unit.conv.ref <- 1 # optional unit conversion for reference data
    variable.unit <- 'kgC m$^{-2}$' # unit after conversion (LaTeX notation)
    result <- scores.grid.notime(long.name, nc.mod, nc.ref, mod.id, ref.id,
                                 unit.conv.mod, unit.conv.ref, variable.unit)
    # Sample result
    myMean <- mean(result[[1]])
    sample <- myMean[155:160]
    computedValue <- sample
    documentedValue <- c(6.770183, 6.741967, 1.306321, 9.493277, 10.263707, 10.457574)

    expect_equal(computedValue, documentedValue, tolerance = 0.001)
  }
)

#-------------------------------------------------------------------------------

# scores.fluxnet.csv

# ------------------------------------------------------------------------------

test_that(
  "Test scores.fluxnet.csv",
  {
    skip_on_cran()
    long.name <- 'Gross primary productivity'
    nc.mod <- system.file('extdata/modelRegular', 'gpp_monthly.nc', package = 'amber')
    ref.csv <- system.file('extdata/referenceRegular', 'gpp_monthly_fluxnet.csv', package = 'amber')
    mod.id <- 'CLASSIC' # define a model experiment ID
    ref.id <- 'FLUXNET' # give reference dataset a name
    unit.conv.mod <- 86400*1000 # optional unit conversion for model data
    unit.conv.ref <- 1 # optional unit conversion for reference data
    variable.unit <- 'gC m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
    scores.fluxnet.csv(long.name, nc.mod, ref.csv, mod.id, ref.id,
                       unit.conv.mod, unit.conv.ref, variable.unit, outputDir = tempdir())

    result <- read.table(paste(tempdir(), "scorevalues_GPP_CLASSIC_vs_FLUXNET_from_2000-01_to_2002-12", sep = "/"))
    computedValue <- result$S_overall
    documentedValue <- c(0.6712645, 0.6230724)

    expect_equal(computedValue, documentedValue, tolerance = 0.001)
  }
)

#-------------------------------------------------------------------------------

# scores.runoff

# ------------------------------------------------------------------------------

test_that(
  "Test scores.runoff",
  {
    skip_on_cran()
    long.name <- 'Streamflow'
    nc.mod <- system.file('extdata/modelRegular', 'mrro_monthly.nc', package = 'amber')
    nc.ref <- system.file('extdata/referenceRegular', 'runoff.nc', package = 'amber')
    nc.basins <- system.file('extdata/referenceRegular', 'basins.nc', package = 'amber')
    mod.id <- 'CLASSIC' # model name
    ref.id <- 'GRDC' # give reference dataset a name
    unit.conv.mod <- 86400 # optional unit conversion for model data
    unit.conv.ref <- 86400 # optional unit conversion for reference data
    variable.unit <- 'kg m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
    scores.runoff(long.name, nc.mod, nc.ref, nc.basins, mod.id, ref.id, unit.conv.mod,
                  unit.conv.ref, variable.unit, outputDir = tempdir())

    result <- read.table(paste(tempdir(), "scorevalues_MRRO_CLASSIC_vs_GRDC_from_2000-01_to_2002-12", sep = "/"))
    computedValue <- result$S_overall
    documentedValue <- c(0.5142253, 0.4490813)

    expect_equal(computedValue, documentedValue, tolerance = 0.001)
  }
)
