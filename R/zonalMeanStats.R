################################################################################
#' Zonal mean plots of AMBER results (bias, bias scores, etc)
#' @description This function computes zonal mean values of model and reference data and the zonal mean bias, centralized root-mean-squre error, phase, inter-annual variability, and corresponding scores. The computation is based on the NetCDF files produced by \link{scores.grid.time}.
#' @param inputDir A string that gives the location of NetCDF files produced by \link{scores.grid.time}, e.g. '/home/project/study'.
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @return A list with two tables. The first table gives the zonal mean values of centralized root-mean-squre error, phase, inter-annual variability, and corresponding
#' scores for each variable and globally gridded reference data set. The second table gives the physical units of each variable in LaTeX notation (e.g. 'W m$^{-2}$').
#' Both tables are written to two text files (zonalMeanStats and zonalMeanStatsUnits) if the user specifies an output directory.
#'
#' @examples
#'
#' library(amber)
#' library(foreach)
#' library(ncdf4)
#' library(raster)
#'
#' inputDir <- paste(system.file('extdata', package = 'amber'), 'zonalMeanStats', sep = '/')
#' zonalMeanStats(inputDir, outputDir = FALSE)
#'
#' @export
zonalMeanStats <- function(inputDir, outputDir = FALSE) {

    # netcdf files
    nc.mod.mean <- list.files(path = inputDir, pattern = "mod-mean.nc")
    nc.ref.mean <- list.files(path = inputDir, pattern = "ref-mean.nc")
    nc.bias <- list.files(path = inputDir, pattern = "-bias.nc")
    nc.crmse <- list.files(path = inputDir, pattern = "-crmse.nc")
    nc.phase <- list.files(path = inputDir, pattern = "-phase.nc")
    nc.iav <- list.files(path = inputDir, pattern = "-iav.nc")
    nc.score <- list.files(path = inputDir, pattern = "-score.nc")
    #---------------------------------------------------------------------------
    # get physical units of biases
    bias.files <- list(nc.bias)
    bias.files <- unlist(bias.files)
    names <- bias.files
    names <- gsub(".nc", "", names)
    names <- gsub("-", ".", names)
    #
    unitList <- foreach::foreach(i = 1:length(bias.files)) %do% {
        nc.file <- unlist(bias.files[i])
        nc.file <- paste(inputDir, nc.file, sep = "/")
        # get units
        nc <- ncdf4::nc_open(nc.file)
        v2 <- nc$var[[2]]
        units <- v2$units
        variable <- gsub(pattern = ".bias", replacement = "", names[i])
        units <- data.frame(variable, units)
        colnames(units) <- c("variables", "unit")
        assign(paste("units", names[i], sep = "."), units)
    }
    units <- Reduce(rbind, unitList)
    #---------------------------------------------------------------------------
    # get data from all variables
    all.files <- list(nc.mod.mean, nc.ref.mean, nc.bias, nc.crmse, nc.phase, nc.iav, nc.score)
    all.files <- unlist(all.files)
    names <- all.files
    names <- gsub(".nc", "", names)
    names <- gsub("-", ".", names)
    #
    myList <- foreach::foreach(i = 1:length(all.files)) %do% {
        nc.file <- unlist(all.files[i])
        # get the correponding bias file to create a mask that omits all values that mod and ref do not have in common
        parts <- unlist(strsplit(nc.file, "-"))
        bias.file.name <- paste(parts[1], parts[2], "bias.nc", sep = "-")
        nc.bias <- paste(inputDir, bias.file.name, sep = "/")
        mask <- raster::raster(nc.bias)
        mask <- mask - mask + 1
        # get zone
        nc.file <- paste(inputDir, nc.file, sep = "/")
        data <- raster::raster(nc.file)
        data <- data * mask  # mask out all grid cells that mod and ref do not have in commmon
        z <- data
        xy <- sp::coordinates(data)
        lat <- xy[, 2]
        z[] <- lat
        # compute zonal mean values
        zonal.mean <- raster::zonal(data, z, "mean", digits = 5)
        # convert to data frame
        zonal.mean <- data.frame(zonal.mean)
        colnames(zonal.mean) <- c("zone", "mean")
        assign(paste("zonal.mean", names[i], sep = "."), zonal.mean)
    }
    data <- Reduce(function(x, y) merge(x, y, by = "zone"), myList)
    colnames(data) <- c("zone", names)
    #---------------------------------------------------------------------------
    if (outputDir != FALSE) {
        utils::write.table(data, paste(outputDir, "zonalMeanStats", sep = "/"))
        utils::write.table(units, paste(outputDir, "zonalMeanStatsUnits", sep = "/"))
    }
    return(list(data, units))
}

utils::globalVariables("%do%")
