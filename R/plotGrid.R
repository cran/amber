################################################################################ 
#' Plots raster layers of a raster stack object
#' @description This function plots the results from \link{scores.grid.time} and
#' \link{scores.grid.notime}.
#' @param long.name A string that gives the full name of the variable, e.g. 'Gross primary productivity'
#' @param plot.me A list that is produced by \link{scores.grid.time} or
#' \link{scores.grid.notime}.
#' @param irregular logical: TRUE if data is on an irregular grid and FALSE if
#' data is on a regular grid
#' @param my.projection A string that gives the projection of the irregular grid
#' @param shp.filename A string that gives the coastline shapefile
#' @param my.xlim An R object that gives the longitude range that you wish to
#' plot, e.g. c(-180, 180)
#' @param my.ylim An R object that gives the longitude range that you wish to
#' plot, e.g. c(-90, 90)
#' @param plot.width Number that gives the plot width, e.g. 8
#' @param plot.height Number that gives the plot height, e.g. 4
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @return Figures in PDF format.
#' This may include the model data
#' (mean, \eqn{mod.mean}; interannual-variability, \eqn{mod.iav}; month of
#' annual cycle maximum, \eqn{mod.max.month}),
#' the reference data
#' (mean, \eqn{ref.mean}; interannual-variability, \eqn{ref.iav}; month of
#' annual cycle maximum, \eqn{ref.max.month}),
#' statistical metrics
#' (bias, \eqn{bias}; root mean square error, \eqn{rmse}; time difference of the
#'  annual cycle maximum, \eqn{phase}),
#' and scores
#' (bias score, \eqn{bias.score}; root mean square error score,
#' \eqn{rmse.score}; inter-annual variability score \eqn{iav.score};
#' annual cycle score (\eqn{phase.score}).

#' @examples
#'
#' \donttest{
#' # (1) Global plots on a regular grid
#'
#' library(amber)
#' library(doParallel)
#' library(foreach)
#' library(ncdf4)
#' library(parallel)
#' library(raster)
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
#' plotGrid(long.name, plot.me)
#'
#' # Additional parameters:
#' score.weights <- c(1,2,1,1,1) # score weights of S_bias, S_rmse, S_phase, S_iav, S_dist
#' outlier.factor <- 1
#' irregular <- FALSE
#' my.projection <- '+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.'
#'
#' plot.me <- scores.grid.time(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit, score.weights, outlier.factor, irregular,
#' my.projection)
#' plotGrid(long.name, plot.me)
#'
#' # (2) Regional plots on a rotated grid
#' long.name <- 'Gross primary productivity'
#' nc.mod <- system.file('extdata/modelRotated', 'gpp_monthly.nc', package = 'amber')
#' nc.ref <- system.file('extdata/referenceRotated', 'gpp_GBAF_rotated.nc', package = 'amber')
#' mod.id <- 'CLASSIC' # define a model experiment ID
#' ref.id <- 'GBAF' # give reference dataset a name
#' unit.conv.mod <- 86400*1000 # optional unit conversion for model data
#' unit.conv.ref <- 86400*1000 # optional unit conversion for reference data
#' variable.unit <- 'gC m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
#' score.weights <- c(1,2,1,1,1) # score weights of S_bias, S_rmse, S_phase, S_iav, S_dist
#' outlier.factor <- 10
#' irregular <- TRUE
#' my.projection <- '+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.'
#'
#' plot.me <- scores.grid.time(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit, score.weights, outlier.factor, irregular,
#' my.projection)
#'
#' # Plot results:
#' irregular <- TRUE # data is on an irregular grid
#' my.projection <- '+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.'
#' shp.filename <- system.file('extdata/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp',
#'  package = 'amber')
#' my.xlim <- c(-171, 0) # longitude range that you wish to plot
#' my.ylim <- c(32, 78) # latitude range that you wish to plot
#' plot.width <- 7 # plot width
#' plot.height <- 3.8 # plot height
#'
#' plotGrid(long.name, plot.me, irregular, my.projection,
#' shp.filename, my.xlim, my.ylim, plot.width, plot.height)
#' }
#' @export
plotGrid <- function(long.name, plot.me, irregular = FALSE, my.projection = "+proj=longlat +ellps=WGS84", shp.filename = system.file("extdata/ne_110m_land/ne_110m_land.shp", 
    package = "amber"), my.xlim = c(-180, 180), my.ylim = c(-60, 85), plot.width = 8, plot.height = 3.8, outputDir = FALSE) {
    
    land <- intFun.coast(my.xlim, my.ylim, my.projection, shp.filename)  # reproject coastline
    # loop through all layers of the raster stack
    
    mod.outlier.points <- plot.me[[2]]
    ref.outlier.points <- plot.me[[3]]
    if (length(plot.me) > 3) 
        bias.significance <- plot.me[[4]]
    
    plot.me <- plot.me[[1]]
    for (i in 1:raster::nlayers(plot.me)) {
        data <- plot.me[[i:i]]
        # get metadata
        meta <- raster::metadata(data)
        my.filename <- unlist(meta[1])
        id <- unlist(meta[2])
        my.title <- gsub("_", " ", id)
        min.max.int <- unlist(meta[3])
        legend.bar.text <- latex2exp::TeX(meta[[4]])
        
        # for legend
        min <- min.max.int[1]
        max <- min.max.int[2]
        interval <- min.max.int[3]
        my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
        my.labels <- round(seq(min, max, interval), 3)  # location of labels
        my.col <- viridis::viridis(n = length(my.breaks) - 1, direction = -1)
        # color options: 'magma', 'inferno', 'plasma', 'viridis', 'cividis'
        my.col.bias <- scico::scico(n = length(my.breaks) - 1, palette = "vik")
        my.col.phase <- grDevices::rainbow(n = length(my.breaks) - 1)
        if (i == 3) 
            {
                my.col <- my.col.bias
            }  # divergent color scheme for bias plots
        if (i == 7 || i == 8) 
            {
                my.col <- my.col.phase
            }  # circular color scheme for phase plots
        my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
        my.legend.args <- list(text = legend.bar.text, side = 2, font = 1, line = 1, cex = 1)
        # plot name
        my.filename <- gsub("_", "-", my.filename)
        my.filename <- gsub(".", "-", my.filename, fixed = TRUE)
        # mod.mean
        dummy <- stats::runif(360 * 180, min = min, max = max)
        dummy <- matrix(dummy, nrow = 180)
        dummy <- raster::raster(dummy)
        
        if (irregular == FALSE) {
            myExtent <- c(-180, 180, -90, 90)
        }
        if (irregular == TRUE) {
            myExtent <- c(0.8553854, 2.041861, -0.1171118, 0.4934048)
        }
        
        raster::extent(dummy) <- myExtent
        # plot
        oldpar <- graphics::par(mfrow = c(1, 2))
        on.exit(graphics::par(oldpar))
        
        if (outputDir != FALSE) {
            grDevices::pdf(paste(outputDir, "/", my.filename, ".pdf", sep = ""), width = plot.width, height = plot.height)
        }
        graphics::par(font.main = 1, mar = c(3, 3, 3, 4), lwd = 1, cex = 1)
        my.xlim <- c(raster::extent(land)[1], raster::extent(land)[2])
        my.ylim <- c(raster::extent(land)[3], raster::extent(land)[4])
        raster::plot(dummy, col = NA, legend = FALSE, xlim = my.xlim, ylim = my.ylim, main = paste(long.name, my.title, 
            sep = "\n"), axes = FALSE)
        raster::plot(land, col = "grey", border = NA, add = TRUE)
        raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, add = TRUE)
        # mark grid cells with extreme outliers
        if (i == 1) 
            {
                graphics::points(mod.outlier.points, pch = 16, cex = 1, col = "red")
            }  # mod.mean
        if (i == 2) 
            {
                graphics::points(ref.outlier.points, pch = 16, cex = 1, col = "red")
            }  # ref.mean
        # mark grid cells where bias is not statistically significant
        if (i == 3 & exists("bias.significance") == TRUE) 
            {
                bias.significance <- raster::rasterToPoints(bias.significance)
                if (irregular == FALSE) {
                  graphics::points(bias.significance, pch = 16, cex = 0.2)
                }
                if (irregular == TRUE) {
                  graphics::points(bias.significance, pch = 16, cex = 0.1)
                }
            }  # non-significant differences are marked with dots
        
        plot(land, add = TRUE)
        if (irregular == FALSE) {
            graphics::axis(1, labels = TRUE, tcl = 0.3)
        }
        if (irregular == FALSE) {
            graphics::axis(2, labels = TRUE, tcl = 0.3, las = 2)
        }
        plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args, 
            legend.width = 1.5, legend.shrink = 1, font = 1)
        if (outputDir != FALSE) {
            grDevices::dev.off()
        }
    }
}
if (getRversion() >= "2.15.1") utils::globalVariables(c("axis", "dev.off", "plot"))
