################################################################################
#' Correlation matrix for statistical metrics computed by AMBER
#' @description This function produces a correlation matrices for \eqn{bias}, \eqn{rmse},
#' \eqn{crmse}, \eqn{phase}, and corresponding scores. The input data consist of netCDF files
#' produced by \link{scores.grid.time} and \link{scores.grid.notime}.
#' @param metric A string that indicates for what statistical metric the correlation matrix should be computed.
#' Options are 'bias', 'crmse', 'phase', 'bias-score', 'crmse-score', 'phase-score', or 'iav-score'.
#' @param inputDir A string that gives the location of NetCDF files produced by \link{scores.grid.time}, e.g. '/home/project/study'.
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @param significanceLevel A number that gives the desired significance level of a correlation, e.g. 0.01
#' @param plot.width  A number that gives the plot width, e.g. 8
#' @param plot.height  A number that gives the plot height, e.g. 6.8
#' @param plot.margin An R object that gives the plot margin, e.g. c(10, 10, 1, 4)
#' @return A list with the Spearman correlation coefficient and corresponding p-values, and a Figure of the correlation matrix
#'
#' @examples
#'
#' library(amber)
#' library(foreach)
#' library(Hmisc)
#' library(ncdf4)
#' library(raster)
#'
#' inputDir <- paste(system.file('extdata', package = 'amber'), 'zonalMeanStats', sep = '/')
#' correlationMatrix(metric = 'bias', inputDir)
#'
#' @export
correlationMatrix <- function(metric, inputDir, outputDir = FALSE, significanceLevel = 0.01, plot.width = 8, plot.height = 6.8,
    plot.margin = c(10, 10, 1, 4)) {

    # netcdf files
    if (metric == "bias") {
        nc.files <- list.files(path = inputDir, pattern = "-bias.nc")
    }
    if (metric == "crmse") {
        nc.files <- list.files(path = inputDir, pattern = "-crmse.nc")
    }
    if (metric == "phase") {
        nc.files <- list.files(path = inputDir, pattern = "-phase.nc")
    }
    # if (metric == 'iav') { nc.files <- list.files(path = inputDir, pattern = '-iav.nc') }
    if (metric == "bias-score") {
        nc.files <- list.files(path = inputDir, pattern = "-bias-score.nc")
    }
    if (metric == "rmse-score") {
        nc.files <- list.files(path = inputDir, pattern = "-rmse-score.nc")
    }
    if (metric == "phase-score") {
        nc.files <- list.files(path = inputDir, pattern = "-phase-score.nc")
    }
    if (metric == "iav-score") {
        nc.files <- list.files(path = inputDir, pattern = "-iav-score.nc")
    }
    nc.files <- unlist(nc.files)
    #---------------------------------------------------------------------------
    myList <- foreach::foreach(i = 1:length(nc.files)) %do% {
        # convert file names into variable names (e.g. 'AGLBIO-GEOCARBON-bias.nc' -> 'AGLBIO.GEOCARBON')
        names <- nc.files
        names <- gsub(".nc", "", names)
        names <- gsub(paste("-", metric, sep = ""), "", names)
        names <- gsub("-", ".", names)
        # get file
        nc.mod <- unlist(nc.files[i])
        nc.mod <- paste(inputDir, nc.mod, sep = "/")
        # get values
        data <- raster::raster(nc.mod)
        lonLat <- sp::coordinates(data)
        lon <- lonLat[, 1]
        lat <- lonLat[, 2]
        lon <- round(lon, 3)
        lat <- round(lat, 3)
        data <- raster::values(data)
        # convert to data frame
        data <- data.frame(lon, lat, data)
        colnames(data) <- c("lon", "lat", names[i])
        assign(paste("dataValues", names[i], sep = "."), data)
    }
    data <- Reduce(function(x, y) merge(x, y, by = c("lon", "lat")), myList)
    data <- subset(data, select = -c(lon, lat))  # drop lon and lat
    myLabel <- colnames(data)  # get lables
    myLabel <- gsub(".bias", "", myLabel)
    myLabel <- gsub(".crmse", "", myLabel)
    myLabel <- gsub(".phase", "", myLabel)
    # myLabel <- gsub('.iav', '', myLabel)
    myLabel <- gsub(".bias.score", "", myLabel)
    myLabel <- gsub(".crmse.score", "", myLabel)
    myLabel <- gsub(".phase.score", "", myLabel)
    myLabel <- gsub(".iav.score", "", myLabel)

    #---------------------------------------------------------------------------
    # correlation matrix
    #---------------------------------------------------------------------------
    cmatrix <- Hmisc::rcorr(as.matrix(data), type = "spearman")  # type can be pearson or spearman
    CorrCoef <- cmatrix$r
    CorrCoef[upper.tri(CorrCoef, diag = TRUE)] <- NA  # omit upper triangle
    # n <- cmatrix$n
    pValue <- cmatrix$P
    pValue[upper.tri(pValue)] <- NA  # omit upper triangle
    # rename CorrCoref and pValue to return them at the end of this function
    CorrCorefTable <- CorrCoef
    pValueTable <- pValue
    # convert tables to raster
    CorrCoef <- raster::raster(CorrCoef)
    pValue <- raster::raster(pValue)
    # set not significant correlations to NA
    fun <- function(x) {
        x[x > significanceLevel] <- NA
        return(x)
    }
    pValue <- raster::calc(pValue, fun)  #
    significant <- pValue - pValue + 1  # give significant corrlations a value of 1
    # mask out insignificant correlations
    CorrCoef.sig <- CorrCoef * significant
    n <- length(myLabel)
    raster::extent(CorrCoef) <- c(0, n, 0, n)
    raster::extent(CorrCoef.sig) <- c(0, n, 0, n)
    data <- CorrCoef
    data.sig <- CorrCoef.sig
    # prepare plot inputs
    my.breaks <- seq(-1, 1, 0.1)
    my.labels <- seq(-1, 1, 0.2)
    my.col <- viridis::plasma(n = length(my.breaks)/2 - 1, direction = 1)
    my.col <- c(my.col, rev(my.col))
    my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
    my.legend.args <- list(text = paste("Spearman Correlation Coefficient (", metric, ")", sep = ""), side = 2, font = 1,
        line = 1, cex = 1)
    #---------------------------------------------------------------------------
    # plot
    #---------------------------------------------------------------------------
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", "correlationMatrix-", metric, ".pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(1, 1), font.main = 1, mar = plot.margin, lwd = 1, cex = 1)
    raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = NA, axes = FALSE, box = FALSE)
    raster::text(data.sig, digits = 1, cex = 0.7)
    graphics::axis(side = 2, at = rev(seq(0.5, n, 1)), labels = myLabel, las = 2)
    graphics::axis(side = 1, at = seq(0.5, n, 1), labels = myLabel, las = 2)
    raster::plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args,
        legend.width = 1.5, legend.shrink = 1, font = 1)
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
    return(list(CorrCorefTable, pValueTable))
}
utils::globalVariables("%do%")
