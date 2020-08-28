################################################################################
#' Plot a matrix that shows the difference between two correlation
#' coefficient matrices computed by \link{correlationMatrix}.
#' @description This function plots a matrix that shows the difference between two correlation
#' coefficient matrices computed by \link{correlationMatrix}. This is useful for assessing how well
#' model data reproduces correlations that are evident in reference data. The difference is computed
#' as the absolute value of the first correlation matrix minus the absolute value of the second
#' correlation matrix.
#' @param cm.one An R objective that gives a correlation matrix computed by \link{correlationMatrix}.
#' @param cm.two An R objective that gives another correlation matrix computed by \link{correlationMatrix}.
#' @param myRows Optional: the user can highlight relations between variables by
#' specifying variable names along the rows and columns of a matrix. Those relations
#' are then highlighted by plotting a polygon around the corresponding grid cells
#' and by plotting the corresponding value.
#' @param myColumns Optional: Same as myRows but for variable names listed along the columns of the matrix.
#' @param inputDir A string that gives the location of NetCDF files produced by \link{scores.grid.time}, e.g. '/home/project/study'.
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @param ofileName A string that gives the output file name, e.g. 'myOutput.pdf'
#' @param plot.width  A number that gives the plot width, e.g. 8
#' @param plot.height  A number that gives the plot height, e.g. 6.8
#' @param plot.margin An R object that gives the plot margin, e.g. c(10, 10, 1, 4)
#' @return A Figure of a matrix that shows the difference between two correlation
#' coefficient matrices.
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
#' inputDir <- paste(system.file('extdata', package = 'amber'), 'zonalMeanStats', sep = '/')
#' cm.one <- correlationMatrix(metric = 'bias', inputDir = inputDir)
#' cm.two <- correlationMatrix(metric = 'bias', inputDir = inputDir)
#'
#' correlationMatrixDiff(cm.one, cm.two, inputDir = inputDir)
#'
#' # You can specify certain relationships to highlight them in your correlation matrix
#' myRows <- c('LAI.AVHRR', 'LAI.AVHRR')
#' myColumns <- c('ALBS.CERES', 'GPP.FluxCom')
#'
#' correlationMatrixDiff(cm.one, cm.two, inputDir = inputDir, myRows = myRows, myColumns = myColumns)
#'
#' @export
correlationMatrixDiff <- function(cm.one, cm.two, myRows = NA, myColumns = NA, inputDir, outputDir = FALSE, ofileName = "correlationMatrixDiff.pdf",
    plot.width = 8, plot.height = 6.8, plot.margin = c(10, 10, 1, 4)) {

    cm.mod <- data.frame(cm.one[1])
    cm.ref <- data.frame(cm.two[1])

    delta <- abs(cm.mod) - abs(cm.ref)
    delta <- round(delta, 2)

    myLabel <- colnames(delta)
    n <- length(myLabel)

    #---------------------------------------------------------------------------
    # Create a mask that marks user-specified relations in the correlation matrix I use the mask to determine what values will be
    # printed
    mask <- delta
    # the loop sets all usewr-defined relationships to the value two
    if (length(myRows) > 1) {
        for (i in 1:length(myRows)) {
            mask[myRows[i], myColumns[i]] <- 2
        }
    }
    # set all other values to NA
    mask[mask < 2] <- NA
    # convert matrix to raster
    mask <- raster::raster(as.matrix(mask))
    mask <- mask/mask  # unity mask
    #---------------------------------------------------------------------------

    delta <- raster::raster(as.matrix(delta))

    raster::extent(delta) <- c(0, n, 0, n)
    raster::extent(mask) <- c(0, n, 0, n)

    # prepare plot inputs
    my.breaks <- seq(-1, 1, 0.1)
    my.labels <- seq(-1, 1, 0.2)
    my.col <- scico::scico(n = length(my.breaks), direction = 1, palette = "roma")
    my.col <- c(my.col, rev(my.col))
    my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
    my.legend.args <- list(text = "Difference in absolute R", side = 2, font = 1, line = 1, cex = 1)

    #---------------------------------------------------------------------------
    # plot
    #---------------------------------------------------------------------------
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", ofileName, sep = ""), width = plot.width, height = plot.height)
    }

    graphics::par(mfrow = c(1, 1), font.main = 1, mar = plot.margin, lwd = 1, cex = 1)
    raster::plot(delta, col = my.col, breaks = my.breaks, legend = FALSE, main = NA, axes = FALSE, box = FALSE)
    values4text <- mask * delta
    # values4text <- delta

    if (length(myRows) > 1) {

        myMax <- max(abs(raster::values(values4text)), na.rm = TRUE)
        myMin <- min(abs(raster::values(values4text)), na.rm = TRUE)

        maskPolygon <- raster::rasterToPolygons(mask, dissolve = TRUE)
        raster::plot(maskPolygon, lwd = 1, add = TRUE)
        if (myMin < 0.5) {
            raster::text(values4text, digits = 1, fun = function(x) {
                abs(x) < 0.5
            }, col = "black", cex = 0.7)
        }
        if (myMax >= 0.5) {
            raster::text(values4text, digits = 1, fun = function(x) {
                abs(x) >= 0.5
            }, col = "white", cex = 0.7)
        }
    }
    graphics::axis(side = 2, at = rev(seq(0.5, n, 1)), labels = myLabel, las = 2)
    graphics::axis(side = 1, at = seq(0.5, n, 1), labels = myLabel, las = 2)
    raster::plot(delta, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args,
        legend.width = 1.5, legend.shrink = 1, font = 1)
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}
utils::globalVariables("%do%")
