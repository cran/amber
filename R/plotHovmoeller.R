################################################################################
#' Plot Hovmoeller diagrams that show monthly climatological mean values and biases
#' @description This function plots Hovmoeller diagrams of monthly climatological
#' mean values and biases computed by \link{scores.grid.time}.
#' @param plot.me A list that is produced by \link{scores.grid.time}
#' @param long.name A string that gives the full name of the variable, e.g. 'Gross primary productivity'
#' @param mod.id A string that identifies the source of the reference data set, e.g. 'CLASSIC'
#' @param ref.id A string that identifies the source of the reference data set, e.g. 'MODIS'
#' @param myBin An integer number that defines the latitudinal range used for computing the zonal mean.
#' For instance, a value of 10 implies that a zonal mean is computed for every 10 degrees latitude.
#' @param gridCellWidth A number that is used as a factor to adjust the width of grid cells, e.g. 1.
#' @param plot.width Number that gives the plot width, e.g. 8
#' @param plot.height Number that gives the plot height, e.g. 8
#' @param my.ylim An R object with the latitudinal range that should be plotted, e.g. c(-40, 65).
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @return Figures in PDF format.
#' @examples
#'
#' \donttest{
#' # Global plots on a regular grid
#'
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
#' long.name <- 'Gross primary productivity'
#' nc.mod <- system.file('extdata/modelRegular', 'gpp_monthly.nc', package = 'amber')
#' nc.ref <- system.file('extdata/referenceRegular', 'gpp_GBAF_128x64.nc', package = 'amber')
#' mod.id <- 'CLASSIC' # define a model experiment ID
#' ref.id <- 'GBAF' # give reference dataset a name
#' unit.conv.mod <- 86400*1000 # optional unit conversion for model data
#' unit.conv.ref <- 86400*1000 # optional unit conversion for reference data
#' variable.unit <- 'gC m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
#'
#' # Short version using default settings:
#' plot.me <- scores.grid.time(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit)
#' plotHovmoeller(plot.me, long.name, mod.id, ref.id)
#'
#' } #donttest
#' @export
plotHovmoeller <- function(plot.me, long.name, mod.id, ref.id, myBin = 20, gridCellWidth = 2, plot.width = 4, plot.height = 5.2,
    my.ylim = c(-100, 100), outputDir = FALSE) {

    #---------------------------------------------------------------------------
    # bias, model, ref
    #---------------------------------------------------------------------------
    for (i in c(6, 7, 5)) {
        # 5 = bias, 6 = mod, 7 = ref

        # get data
        data <- plot.me[[i]]

        # get metadata
        meta <- raster::metadata(data)
        my.filename <- unlist(meta[1])
        # id <- unlist(meta[2]) id <- gsub('Monthly_', '', id) my.title <- gsub('_', ' ', id)

        myUnit <- meta[[4]]
        myUnit <- gsub("\\", "", myUnit, fixed = TRUE)
        legend.bar.text <- latex2exp::TeX(myUnit)

        # legend.bar.text <- latex2exp::TeX(meta[[4]])

        # make zonal mean
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
        if (i == 5) {
            bias <- data
        }
        if (i == 6) {
            mod <- data
        }
        if (i == 7) {
            ref <- data
        }
    }

    # Show model, reference, and bias in a single Figure.  plot name
    my.filename <- gsub("_", "-", my.filename)
    my.filename <- gsub(".", "-", my.filename, fixed = TRUE)

    # plot inputs
    my.col.lab <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
    lat.lab <- seq(-80, 80, 20)  # axis label for vertical axis (latitude)
    my.legend.args <- list(text = legend.bar.text, side = 2, font = 1, line = 0.1, cex = 0.75)

    # mod and ref
    min.max.int <- intFun.min.max.int.mod.ref(mod, ref)
    min <- min.max.int[1]
    max <- min.max.int[2]
    interval <- min.max.int[3]
    my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
    my.labels <- round(seq(min, max, interval), 3)  # location of labels
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
    my.col <- scico::scico(n = length(my.breaks) - 1, palette = "vik")
    my.col.bias <- my.col
    my.breaks.bias <- my.breaks
    my.axis.args.bias <- list(at = my.labels, labels = my.labels, cex.axis = 1)

    # plot
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))

    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, ".pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(3, 1), font.main = 1, oma = c(2, 0, 2, 0), mar = c(0, 5, 0, 2), lwd = 1, cex = 1, xpd = NA)
    raster::plot(mod, col = my.col.mod.ref, breaks = my.breaks.mod.ref, ylim = my.ylim/10, legend = FALSE, axes = FALSE, ylab = NA,
        cex.lab = 0.8)
    graphics::axis(side = 2, at = lat.lab/10, labels = lat.lab, las = 2)
    graphics::mtext(long.name, side = 3, line = 1, cex = 0.75)
    raster::plot(mod, legend.only = TRUE, col = my.col.mod.ref, breaks = my.breaks.mod.ref, axis.args = my.axis.args.mod.ref,
        legend.args = my.legend.args, legend.width = 1, legend.shrink = 0.9, font = 1, horizontal = FALSE)
    graphics::mtext(mod.id, side = 4, line = -0.5, cex = 0.75)

    #---------------------------------------------------------------------------

    raster::plot(ref, col = my.col.mod.ref, breaks = my.breaks.mod.ref, ylim = my.ylim/10, legend = FALSE, axes = FALSE, ylab = "Degrees Latitude")
    graphics::axis(side = 2, at = lat.lab/10, labels = lat.lab, las = 2)
    raster::plot(ref, legend.only = TRUE, col = my.col.mod.ref, breaks = my.breaks.mod.ref, axis.args = my.axis.args.mod.ref,
        legend.args = my.legend.args, legend.width = 1, legend.shrink = 0.9, font = 1, horizontal = FALSE)
    graphics::mtext(ref.id, side = 4, line = -0.5, cex = 0.75)

    #---------------------------------------------------------------------------

    raster::plot(bias, col = my.col.bias, breaks = my.breaks.bias, ylim = my.ylim/10, legend = FALSE, axes = FALSE, ylab = NA)
    graphics::axis(side = 1, at = seq(gridCellWidth/2, gridCellWidth * ncol(data), gridCellWidth), labels = my.col.lab, cex.axis = 0.7,
        las = 1)
    graphics::axis(side = 2, at = lat.lab/10, labels = lat.lab, las = 2)
    raster::plot(bias, legend.only = TRUE, col = my.col.bias, breaks = my.breaks.bias, axis.args = my.axis.args.bias, legend.args = my.legend.args,
        legend.width = 1, legend.shrink = 0.9, font = 1, horizontal = FALSE)
    graphics::mtext("Bias", side = 4, line = -0.5, cex = 0.75)

    #---------------------------------------------------------------------------

    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}
if (getRversion() >= "2.15.1") utils::globalVariables(c("axis", "dev.off", "plot"))
