################################################################################
#' Correlation matrix for statistical metrics computed by AMBER for FLUXNET data
#' @description This function produces a correlation matrices for mean values, \eqn{bias},
#' \eqn{crmse}, \eqn{phase}, and corresponding scores. The input data consist of text files
#' produced by \link{scores.fluxnet.csv}.
#' @param metric A string that indicates for what statistical metric the correlation matrix should be computed.
#' Options are 'mod.mean', ref.mean', bias', 'crmse', 'phase', 'bias.score', 'crmse.score', 'phase.score', or 'iav.score'.
#' @param inputDir A string that gives the location of text files produced by \link{scores.fluxnet.csv}, e.g. '/home/project/study'.
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @param fileNames An object of strings that give the filenames that should be included. The default is c('GPP_FLUXNET', 'HFLS_FLUXNET', 'HFSS_FLUXNET', 'NEE_FLUXNET', 'RECO_FLUXNET', 'RNS_FLUXNET')
#' @param significanceLevel A number that gives the desired significance level of a correlation, e.g. 0.01
#' @param plot.width  A number that gives the plot width, e.g. 8
#' @param plot.height  A number that gives the plot height, e.g. 6.8
#' @param plot.margin An R object that gives the plot margin, e.g. c(10, 10, 1, 4)
#' @return A list with the Spearman correlation coefficient and corresponding p-values, and a Figure of the correlation matrix
#'
#' @examples
#' library(amber)
#' library(classInt)
#' library(doParallel)
#' library(foreach)
#' library(Hmisc)
#' library(latex2exp)
#' library(ncdf4)
#' library(parallel)
#' library(raster)
#' library(rgdal)
#' library(rgeos)
#' library(scico)
#' library(sp)
#' library(stats)
#' library(utils)
#' library(viridis)
#' library(xtable)
#'
#' inputDir <- paste(system.file('extdata', package = 'amber'), 'scores', sep = '/')
#' correlationMatrixFluxnet(metric = 'bias', inputDir)
#'
#' @export
correlationMatrixFluxnet <- function(metric, inputDir, outputDir = FALSE, fileNames = c("GPP_FLUXNET", "HFLS_FLUXNET", "HFSS_FLUXNET",
    "NEE_FLUXNET", "RECO_FLUXNET", "RNS_FLUXNET"), significanceLevel = 0.01, plot.width = 8, plot.height = 6.8, plot.margin = c(10,
    10, 1, 4)) {
    #---------------------------------------------------------------------------
    myList <- foreach::foreach(i = 1:length(fileNames)) %do% {
        fileName <- fileNames[i]
        data <- utils::read.table(paste(inputDir, fileName, sep = "/"))
        fileName <- gsub("_", ".", fileName)
        # get file
        lon <- data$lon
        lat <- data$lat
        lon <- round(lon, 3)
        lat <- round(lat, 3)

        if (metric == "mod.mean") {
            data <- data$mod.mean
        }
        if (metric == "ref.mean") {
            data <- data$ref.mean
        }
        if (metric == "bias") {
            data <- data$bias
        }
        if (metric == "crmse") {
            data <- data$crmse
        }
        if (metric == "phase") {
            data <- data$phase
        }
        if (metric == "bias.score") {
            data <- data$bias.score
        }
        if (metric == "rmse.score") {
            data <- data$rmse.score
        }
        if (metric == "phase.score") {
            data <- data$phase.score
        }
        if (metric == "iav.score") {
            data <- data$iav.score
        }

        data <- data.frame(lon, lat, data)
        colnames(data) <- c("lon", "lat", fileName)
        assign(paste("dataValues", fileName, sep = "."), data)
    }
    data <- Reduce(function(x, y) merge(x, y, by = c("lon", "lat")), myList)
    data <- subset(data, select = -c(lon, lat))  # drop lon and lat
    myLabel <- colnames(data)  # get lables

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
    my.legend.args <- list(text = paste("Spearman Correlation Coefficient (", metric, ")", sep = ""), side = 2, font = 1, line = 1,
        cex = 1)
    #---------------------------------------------------------------------------
    # plot
    #---------------------------------------------------------------------------
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", "correlationMatrixFluxnet-", metric, ".pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(1, 1), font.main = 1, mar = plot.margin, lwd = 1, cex = 1)
    raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = NA, axes = FALSE, box = FALSE)
    raster::text(data.sig, digits = 1, cex = 0.7)
    # raster::text(data.sig, digits = 1, fun = function(x) { abs(x) < 0.5 }, cex = 1) raster::text(data.sig, digits = 1, fun =
    # function(x) { abs(x) >= 0.5 }, col = 'white', cex = 1)
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
