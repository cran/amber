################################################################################
#' Ensemble mean plots of AMBER results (bias, bias scores, etc)
#' @description This function plots ensemble mean, minimum, and maximum values of a statistical
#'  metric computed by \link{scores.grid.time} and \link{scores.grid.notime}.
#' @param long.name A string that gives the full name of the variable, e.g. 'Gross primary productivity'
#' @param metric A string that specifies what statistical metric should be plotted.
#' This includes for instance 'bias', 'crmse', 'phase', 'iav', 'bias-score', 'rmse-score', 'phase-score', and 'iav-score'.
#' @param mod.path.list A List of directories where AMBER output is stored for different model runs,
#' e.g. list(mod01.path, mod02.path, mod03.path)
#' @param modelIDs An R object with the different model run IDs, e.g. c('CLASSIC.CRUJRAv2', 'CLASSIC.GSWP3W5E5', 'CLASSIC.CRUNCEP')
#' @param myVariables  An R object with the variable names of interest, e.g. c('GPP.FluxCom', 'RECO.FluxCom').
#' @param shp.filename A string that gives the coastline shapefile
#' @param my.xlim An R object that gives the longitude range that you wish to plot, e.g. c(-180, 180)
#' @param my.ylim An R object that gives the longitude range that you wish to plot, e.g. c(-90, 90)
#' @param plot.width Number that gives the plot width, e.g. 8
#' @param plot.height Number that gives the plot height, e.g. 4
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @param subcaption A string that defines the subcaption of the figure, e.g. '(a)'.
#' @return Figures in PDF format.
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
#' long.name <- 'Gross Primary Productivity'
#' metric <- 'mod-mean'
#'
#' mod01.path <- paste(system.file('extdata', package = 'amber'), 'model01', sep = '/')
#' mod02.path <- paste(system.file('extdata', package = 'amber'), 'model02', sep = '/')
#' mod.path.list <- list(mod01.path, mod02.path)
#' modelIDs <- c('CLASSIC.CRUJRAv2', 'CLASSIC.GSWP3W5E5')
#'
#' myVariables <- c('GPP-GOSIF', 'GPP-MODIS')
#'
#' plotEnsembleMean(long.name, metric, mod.path.list, modelIDs, myVariables,
#'  plot.width = 5, plot.height = 5.5)
#'
#' @export
plotEnsembleMean <- function(long.name, metric, mod.path.list, modelIDs, myVariables, shp.filename = system.file("extdata/ne_110m_land/ne_110m_land.shp",
    package = "amber"), my.xlim = c(-180, 180), my.ylim = c(-60, 85), plot.width = 5, plot.height = 7, outputDir = FALSE, subcaption = "") {

    # coastline
    my.projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #"+proj=longlat +ellps=WGS84"
    land <- intFun.coast(my.xlim, my.ylim, my.projection, shp.filename)
    # get meta data
    filePath <- file.path(mod.path.list[[1]], myVariables[1])
    fileName <- paste(filePath, "-", metric, ".nc", sep = "")
    nc <- ncdf4::nc_open(fileName)
    units <- nc$var[[2]]$units

    units <- gsub("\\", "", units, fixed = TRUE)

    units.modelMean <- latex2exp::TeX(units)

    if (metric == "phase") {
        units <- "months"
    }
    if (metric == "bias-score") {
        units <- "(-)"
    }
    if (metric == "rmse-score") {
        units <- "(-)"
    }
    if (metric == "phase-score") {
        units <- "(-)"
    }
    if (metric == "iav-score") {
        units <- "(-)"
    }

    vname <- nc$var[[2]]$name
    legend.bar.text <- latex2exp::TeX(units)

    # number of ensemble members and evaluations per ensemble member
    nmod <- length(mod.path.list)
    nvar <- length(myVariables)

    # Loop each ensemble member
    eachModel <- foreach::foreach(m = 1:nmod) %do% {

        # Loop each evaluation
        eachVariable <- foreach::foreach(v = 1:nvar) %do% {

            filePath <- file.path(mod.path.list[[m]], myVariables[v])
            fileName <- paste(filePath, "-", metric, ".nc", sep = "")

            data <- raster::raster(fileName)
            suppressWarnings(raster::projection(data) <- my.projection)
            return(data)
        }
        data <- do.call(raster::stack, eachVariable)
    }
    data <- do.call(raster::stack, eachModel)
    # compute mean value
    ensembleMean <- raster::mean(data, na.rm = TRUE)
    ensembleMin <- min(data, na.rm = TRUE)
    ensembleMax <- max(data, na.rm = TRUE)

    # Get the ensemble mean value

    # Loop each ensemble member
    eachModel <- foreach::foreach(m = 1:nmod) %do% {

        # Loop each evaluation
        eachVariable <- foreach::foreach(v = 1:nvar) %do% {

            filePath <- file.path(mod.path.list[[m]], myVariables[v])
            fileName <- paste(filePath, "-mod-mean.nc", sep = "")

            data <- raster::raster(fileName)
            suppressWarnings(raster::projection(data) <- my.projection)
            return(data)
        }
        data <- do.call(raster::stack, eachVariable)
    }
    data <- do.call(raster::stack, eachModel)
    modelMean <- raster::mean(data, na.rm = TRUE)

    data <- raster::stack(modelMean, ensembleMean, ensembleMin, ensembleMax)

    # plot data
    meanEnsembleMetric <- paste("mean ensemble", metric, sep = " ")
    minEnsembleMetric <- paste("min ensemble", metric, sep = " ")
    maxEnsembleMetric <- paste("max ensemble", metric, sep = " ")

    my.title <- c("model ensemble mean", meanEnsembleMetric, minEnsembleMetric, maxEnsembleMetric)

    # color legend inputs for mean, min, and max Recall: If I use 4 colors I need 5 breaks and 3 labels
    legendInputs <- foreach::foreach(i = 1:raster::nlayers(data)) %do% {
        min.max.int <- intFun.min.max.int.ext(data[[i]])
        if (metric == "bias" && i > 1) {
            min.max.int <- intFun.min.max.int.bias(data[[i]])
        }
        if (metric == "bias" && i > 2) {
            min.max.int <- intFun.min.max.int.bias(data[[3:4]])  # this ensures that min bias and max bias have the same legend
        }
        min <- min.max.int[1]
        max <- min.max.int[2]
        interval <- min.max.int[3]
        my.breaks <- round(seq(min, max, interval), 3)
        # add two more breaks to account for the outlier colors that are added further below
        my.breaks <- c(min - interval, my.breaks, max + interval)
        my.labels <- round(seq(min, max, interval), 3)  # location of labels
        my.col <- viridis::viridis(n = length(my.breaks) - 3, direction = -1)
        my.col <- c("magenta", my.col, "black")  # add colors for outliers
        my.col.bias <- scico::scico(n = length(my.breaks) - 3, palette = "vik")
        my.col.bias <- c("magenta", my.col.bias, "black")  # add colors for outliers
        if (metric == "bias" && i > 1) {
            my.col <- my.col.bias
        }
        my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
        list(my.col, my.breaks, my.labels, my.axis.args, min, max, interval)
    }

    my.legend.args <- list(text = legend.bar.text, side = 2, font = 1, line = -5, cex = 0.75)
    my.filename <- paste(vname, metric, "ensemble_mean.pdf", sep = "_")

    # plot
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))

    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(4, 1), font.main = 1, oma = c(2, 0, 2, 0), mar = c(0, 3, 0, 1), lwd = 1, cex = 1)
    my.xlim <- c(raster::extent(land)[1], raster::extent(land)[2])
    my.ylim <- c(raster::extent(land)[3], raster::extent(land)[4])

    for (i in 1:raster::nlayers(data)) {

        my.col <- legendInputs[[i]][[1]]
        my.breaks <- legendInputs[[i]][[2]]
        my.labels <- legendInputs[[i]][[3]]
        my.axis.args <- legendInputs[[i]][[4]]

        min <- legendInputs[[i]][[5]]
        max <- legendInputs[[i]][[6]]
        interval <- legendInputs[[i]][[7]]

        plotMe <- data[[i]]
        # Change values of outliers for the purpose of the legend
        plotMe[plotMe > max] <- max + interval/2
        plotMe[plotMe < min] <- min - interval/2

        raster::plot(data[[1]], col = NA, legend = FALSE, xlim = my.xlim, ylim = my.ylim, main = NA, axes = FALSE)

        title <- paste(subcaption, " ", long.name, " (evaluation ensemble ", metric, ", N = ", nmod * nvar, ")", sep = "")

        # ticks
        if (i == 1) {
            graphics::axis(3, labels = title, at = 1, tcl = 0.3)
        }
        if (i > 1) {
            graphics::axis(3, labels = FALSE, tcl = 0.3)
        }
        graphics::axis(2, labels = TRUE, tcl = 0.3, las = 2)
        graphics::axis(3, labels = FALSE, tcl = 0.3)
        graphics::axis(4, labels = FALSE, tcl = 0.3)
        if (i == 4) {
            graphics::axis(1, labels = TRUE, tcl = 0.3)
        }

        graphics::par(new = TRUE)
        raster::plot(land, col = "grey", border = NA, add = TRUE)
        raster::plot(plotMe, col = my.col, breaks = my.breaks, legend = FALSE, add = TRUE)
        raster::plot(land, add = TRUE)
        # graphics::box()

        if (i == 1) {
            my.legend.args.modelMean <- list(text = units.modelMean, side = 2, font = 1, line = -5, cex = 0.75)
            raster::plot(plotMe, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args.modelMean,
                legend.width = 1, legend.shrink = 0.9, font = 1)
        }

        if (i > 1) {
            raster::plot(plotMe, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args,
                legend.width = 1, legend.shrink = 0.9, font = 1)
        }


        graphics::legend(15, -40, my.title[i], bty = "n", cex = 0.9)
        if (i == 1) {
            graphics::legend("bottomleft", c(myVariables, modelIDs), bty = "n", cex = 0.6)
        }
    }
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}

utils::globalVariables("%do%")
