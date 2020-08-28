################################################################################
#' Hovmoeller Diagram for model ensemble
#' @description This function plots Hovmoeller diagrams of monthly climatological
#' mean values and biases computed by \link{scores.grid.time} for multiple ensemble members
#' @param mod.path.list A List of directories where AMBER output is stored for different model runs,
#' e.g. list(mod01.path, mod02.path, mod03.path)
#' @param modelIDs An R object with the different model run IDs, e.g. c('CLASSIC.CRUJRAv2', 'CLASSIC.GSWP3W5E5', 'CLASSIC.CRUNCEP')
#' @param myVariables  An R object with the variable names of interest, e.g. c('GPP.FluxCom', 'RECO.FluxCom').
#' @param myBin An integer number that defines the latitudinal range used for computing the zonal mean.
#' For instance, a value of 10 implies that a zonal mean is computed for every 10 degrees latitude.
#' @param gridCellWidth A number that is used as a factor to adjust the width of grid cells, e.g. 1.
#' @param my.ylim An R object that gives the longitude range that you wish to plot, e.g. c(-90, 90)
#' @param plot.width Number that gives the plot width, e.g. 8
#' @param plot.height Number that gives the plot height, e.g. 4
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @return Figures in PDF format.
#' @examples
#'
#' \donttest{
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
#' mod01.path <- paste(system.file('extdata', package = 'amber'), 'model01', sep = '/')
#' mod02.path <- paste(system.file('extdata', package = 'amber'), 'model02', sep = '/')
#' mod.path.list <- list(mod01.path, mod02.path)
#' modelIDs <- c('CLASSIC.CRUJRAv2', 'CLASSIC.GSWP3W5E5')
#'
#' myVariables <- c('GPP-MODIS', 'GPP-GOSIF')
#'
#' plotEnsembleHovmoeller(mod.path.list = mod.path.list,
#' modelIDs = modelIDs, myVariables = myVariables, myBin = 20, gridCellWidth = 2,
#' my.ylim = c(-100, 100), plot.width = 8.4, plot.height = 5.0)
#' } #donttest
#' @export
plotEnsembleHovmoeller <- function(mod.path.list = mod.path.list, modelIDs = modelIDs, myVariables = myVariables, myBin = 20,
    gridCellWidth = 2, my.ylim = c(-100, 100), plot.width = 8.4, plot.height = 5, outputDir = FALSE) {

    # input files look like this: GPP-GOSIF-bias-clim-mly.nc GPP-GOSIF-mod-clim-mly.nc GPP-GOSIF-ref-clim-mly.nc

    # get meta data
    filePath <- file.path(mod.path.list[[1]], myVariables[1])
    fileName <- paste(filePath, "-mod-clim-mly.nc", sep = "")
    nc <- ncdf4::nc_open(fileName)
    units <- nc$var[[2]]$units
    vname <- nc$var[[2]]$name

    units <- gsub("\\", "", units, fixed = TRUE)
    legend.bar.text <- latex2exp::TeX(units)

    # number of ensemble members and evaluations per ensemble member
    nmod <- length(mod.path.list)
    nvar <- length(myVariables)

    # (a) model Loop each ensemble member
    eachModel <- foreach::foreach(m = 1:nmod) %do% {

        # Loop each evaluation
        eachVariable <- foreach::foreach(v = 1:nvar) %do% {

            filePath <- file.path(mod.path.list[[m]], myVariables[v])
            fileName <- paste(filePath, "-mod-clim-mly.nc", sep = "")
            data <- raster::stack(fileName)
        }
    }
    eachModel <- unlist(eachModel)
    mod.mean <- do.call(raster::overlay, c(eachModel, fun = mean))
    mod.min <- do.call(raster::overlay, c(eachModel, fun = min))
    mod.max <- do.call(raster::overlay, c(eachModel, fun = max))

    # (b) reference Loop each ensemble member
    eachModel <- foreach::foreach(m = 1:nmod) %do% {

        # Loop each evaluation
        eachVariable <- foreach::foreach(v = 1:nvar) %do% {

            filePath <- file.path(mod.path.list[[m]], myVariables[v])
            fileName <- paste(filePath, "-ref-clim-mly.nc", sep = "")
            data <- raster::stack(fileName)
        }
    }
    eachModel <- unlist(eachModel)
    ref.mean <- do.call(raster::overlay, c(eachModel, fun = mean))
    ref.min <- do.call(raster::overlay, c(eachModel, fun = min))
    ref.max <- do.call(raster::overlay, c(eachModel, fun = max))

    # (c) bias Loop each ensemble member
    eachModel <- foreach::foreach(m = 1:nmod) %do% {

        # Loop each evaluation
        eachVariable <- foreach::foreach(v = 1:nvar) %do% {

            filePath <- file.path(mod.path.list[[m]], myVariables[v])
            fileName <- paste(filePath, "-bias-clim-mly.nc", sep = "")
            data <- raster::stack(fileName)
        }
    }
    eachModel <- unlist(eachModel)
    bias.mean <- do.call(raster::overlay, c(eachModel, fun = mean))
    bias.min <- do.call(raster::overlay, c(eachModel, fun = min))
    bias.max <- do.call(raster::overlay, c(eachModel, fun = max))

    myDataList <- list(mod.mean, mod.min, mod.max, ref.mean, ref.min, ref.max, bias.mean, bias.min, bias.max)

    # make Hovmoeller diagram
    eachElement <- foreach::foreach(i = 1:length(myDataList)) %do% {
        data <- myDataList[[i]]
        xy <- sp::coordinates(data)
        lat <- xy[, 2]
        lat <- sort(unique(lat), decreasing = FALSE)

        latBreaks <- seq(-180, 180, myBin)
        myCut <- cut(x = lat, breaks = latBreaks)

        lat.floor <- floor(lat/myBin) * myBin  # new
        lat.ceiling <- ceiling(lat/myBin) * myBin  # new
        latMin <- min(lat.floor)
        latMax <- max(lat.ceiling)

        data <- raster::as.array(data)
        data <- apply(data, c(1, 3), mean, na.rm = TRUE)  # zonal mean

        # compute means across larger steps in latitude (e.g. every 10 degrees lat)
        data <- data.frame(myCut, data)
        data <- stats::aggregate(data[, 2:13], by = list(data$myCut), mean, na.rm = TRUE)

        data <- raster::raster(as.matrix(data[, 2:13]))

        myExtent <- c(0, gridCellWidth * ncol(data), latMin/10, latMax/10)
        raster::extent(data) <- myExtent
        return(data)
    }

    # mean, min, max
    mod <- raster::stack(eachElement[1:3])
    ref <- raster::stack(eachElement[4:6])
    bias <- raster::stack(eachElement[7:9])

    # plot inputs
    my.col.lab <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
    lat.lab <- seq(-80, 80, 20)  # axis label for vertical axis (latitude)
    my.legend.args <- list(text = legend.bar.text, side = 2, font = 1, line = -5, cex = 0.75)

    # mod and ref
    min.max.int <- intFun.min.max.int.mod.ref(mod, ref)
    min <- min.max.int[1]
    max <- min.max.int[2]
    interval <- min.max.int[3]
    my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
    my.labels <- round(seq(min, max, interval), 3)  # location of labels
    rm(min, max)
    my.col <- viridis::viridis(n = length(my.labels) - 1, direction = -1)
    my.col.mod.ref <- my.col
    my.breaks.mod.ref <- my.breaks
    my.axis.args.mod.ref <- list(at = my.labels, labels = my.labels, cex.axis = 1)

    # bias
    min.max.int <- intFun.min.max.int.diff(bias)
    min <- min.max.int[1]
    max <- min.max.int[2]
    interval <- min.max.int[3]
    my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
    my.labels <- round(seq(min, max, interval), 3)  # location of labels
    rm(min, max)
    my.col <- scico::scico(n = length(my.breaks) - 1, palette = "vik")
    my.col.bias <- my.col
    my.breaks.bias <- my.breaks
    my.axis.args.bias <- list(at = my.labels, labels = my.labels, cex.axis = 1)

    # plot
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))

    my.filename <- paste(vname, "EnsembleHovmoeller", sep = "_")

    ensembleInfo <- c("Mean", "Min", "Max")

    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, ".pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(3, 3), font.main = 1, oma = c(2, 3, 2, 3), mar = c(0, 0, 0, 2), lwd = 1, cex = 1)

    for (i in c(1, 2, 3)) {

        raster::plot(mod[[i]], col = my.col.mod.ref, breaks = my.breaks.mod.ref, ylim = my.ylim/10, legend = FALSE, axes = FALSE,
            ylab = NA, cex.lab = 0.8)
        graphics::legend("topleft", ensembleInfo[i], bty = "n")
        if (i == 1) {
            graphics::mtext(paste("(a) ", vname, " Model Ensemble (N = ", nmod, ")"), side = 3, line = 1, cex = 0.75)
        }

        if (i == 2) {
            graphics::mtext("Degrees Latitude", side = 2, line = 3, cex = 0.75)
        }

        graphics::axis(side = 2, at = lat.lab/10, labels = lat.lab, las = 2)

        if (i == 3) {

            graphics::axis(side = 1, at = seq(gridCellWidth/2, gridCellWidth * ncol(data), gridCellWidth), labels = my.col.lab,
                cex.axis = 0.7, las = 1)
        }



        raster::plot(ref[[i]], col = my.col.mod.ref, breaks = my.breaks.mod.ref, ylim = my.ylim/10, legend = FALSE, axes = FALSE,
            ylab = NA, cex.lab = 0.8)
        graphics::legend("topleft", ensembleInfo[i], bty = "n")
        if (i == 1) {
            graphics::mtext(paste("(b) ", vname, " Reference Ensemble (N = ", nvar, ")"), side = 3, line = 1, cex = 0.75)
        }

        if (i == 3) {
            graphics::axis(side = 1, at = seq(gridCellWidth/2, gridCellWidth * ncol(data), gridCellWidth), labels = my.col.lab,
                cex.axis = 0.7, las = 1)
        }


        raster::plot(bias[[i]], col = my.col.bias, breaks = my.breaks.bias, ylim = my.ylim/10, legend = FALSE, axes = FALSE,
            ylab = NA, cex.lab = 0.8)
        graphics::legend("topleft", ensembleInfo[i], bty = "n")
        if (i == 1) {
            graphics::mtext(paste("(c) ", vname, "Mean Ensemble Bias (N = ", nmod * nvar, ")"), side = 3, line = 1, cex = 0.75)
        }

        if (i == 1) {
            raster::plot(mod[[1]], legend.only = TRUE, col = my.col.mod.ref, breaks = my.breaks.mod.ref, axis.args = my.axis.args.mod.ref,
                legend.args = my.legend.args, legend.width = 2, legend.shrink = 1, font = 1, horizontal = FALSE)
        }

        if (i == 3) {
            raster::plot(bias[[1]], legend.only = TRUE, col = my.col.bias, breaks = my.breaks.bias, axis.args = my.axis.args.bias,
                legend.args = my.legend.args, legend.width = 2, legend.shrink = 1, font = 1, horizontal = FALSE)
        }

        if (i == 3) {
            graphics::axis(side = 1, at = seq(gridCellWidth/2, gridCellWidth * ncol(data), gridCellWidth), labels = my.col.lab,
                cex.axis = 0.7, las = 1)
        }

    }

    if (outputDir != FALSE) {
        grDevices::dev.off()
    }


}

utils::globalVariables("%do%")
