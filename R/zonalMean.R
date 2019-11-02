################################################################################ 
#' Zonal mean plots of model and reference data
#' @description This function plots zonal mean values and corresponding inter-quartile
#' ranges of model and reference data.
#' @param long.name A string that gives the full name of the variable, e.g. 'Gross primary productivity'
#' @param nc.mod A string that gives the path and name of the netcdf file that contains the model output, e.g. '/home/model_gpp.nc'
#' @param nc.ref A string that gives the path and name of the netcdf file that contains the reference data output, e.g. '/home/reference_gpp.nc'
#' @param mod.id A string that identifies the source of the reference data set, e.g. 'CLASSIC'
#' @param ref.id A string that identifies the source of the reference data set, e.g. 'MODIS'
#' @param unit.conv.mod A number that is used as a factor to convert the unit of the model data, e.g. 86400
#' @param unit.conv.ref A number that is used as a factor to convert the unit of the reference data, e.g. 86400
#' @param variable.unit A string that gives the final units using LaTeX notation, e.g. 'gC m$^{-2}$ day$^{-1}$'
#' @param outlier.factor A number that is used to define outliers, e.g. 10.
#'  Plotting raster objects that contain extreme outliers lead to figures where
#'  most grid cells are presented by a single color since the color legend covers
#'  the entire range of values. To avoid this, the user may define outliers that
#'  will be masked out and marked with a red dot. Outliers are all values that
#'  exceed the interquartile range multiplied by the outlier factor defined here.
#' @param plot.width Number that gives the plot width, e.g. 8
#' @param plot.height Number that gives the plot height, e.g. 4
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @return A figure in PDF format that gives the zonal mean values of model and
#' reference data. The bold line presents the mean and the shaded area the
#' corresponding interquartile range (IQR). The IQR presents the inter-annual variability
#' and longitudinal variability, combined.
#'
#' @examples
#'
#' library(amber)
#' library(latex2exp)
#' library(ncdf4)
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
#' outlier.factor <- 1000
#'
#' zonalMean(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit, outlier.factor)
#'
#' @export
zonalMean <- function(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod, unit.conv.ref, variable.unit, outlier.factor = 1000, 
    plot.width = 8, plot.height = 3.8, outputDir = FALSE) {
    
    # Data preparation
    nc <- ncdf4::nc_open(nc.mod)
    variable.name <- names(nc[["var"]])
    ncdf4::nc_close(nc)
    variable.name <- variable.name[length(variable.name)]  # take the last variable (relevant for CanESM5)
    variable.name <- ifelse(variable.name == "burntFractionAll", "burnt", variable.name)  # rename burntFractionAll to shorter name
    variable.name <- toupper(variable.name)  # make variable name upper-case
    mod <- raster::brick(nc.mod)
    dates.mod <- raster::getZ(mod)
    dates.mod <- format(as.Date(dates.mod), "%Y-%m")  # only year and month
    start.date.mod <- min(dates.mod)
    end.date.mod <- max(dates.mod)
    
    # reference data
    ref <- raster::brick(nc.ref)
    dates.ref <- raster::getZ(ref)
    dates.ref <- format(as.Date(dates.ref), "%Y-%m")  # only year and month
    start.date.ref <- min(dates.ref)
    end.date.ref <- max(dates.ref)
    
    # find common time period
    start.date <- max(start.date.mod, start.date.ref)
    end.date <- min(end.date.mod, end.date.ref)
    
    # subset common time period
    mod <- mod[[which(format(as.Date(raster::getZ(mod)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(mod)), "%Y-%m") <= 
        end.date)]]
    ref <- ref[[which(format(as.Date(raster::getZ(ref)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(ref)), "%Y-%m") <= 
        end.date)]]
    
    # get layer names
    mod.names <- names(mod)
    ref.names <- names(ref)
    
    # unit conversion if appropriate
    mod <- mod * unit.conv.mod
    ref <- ref * unit.conv.ref
    
    # Extreme outliers are set to NA in the grid
    
    mod.mean <- raster::mean(mod, na.rm = TRUE)  # time mean
    mod.outlier_range <- intFun.grid.define.outlier(mod.mean, outlier.factor)  # define outlier range
    outlier.neg <- mod.outlier_range[1]
    outlier.pos <- mod.outlier_range[2]
    mod.mask_outliers <- intFun.grid.outliers.na(mod.mean, outlier.neg, outlier.pos)
    mod.mask_outliers <- mod.mask_outliers - mod.mask_outliers + 1
    mod <- mod * mod.mask_outliers
    names(mod) <- mod.names
    
    # reference data
    ref.mean <- raster::mean(ref, na.rm = TRUE)  # time mean
    ref.outlier_range <- intFun.grid.define.outlier(ref.mean, outlier.factor)  # define outlier range
    outlier.neg <- ref.outlier_range[1]
    outlier.pos <- ref.outlier_range[2]
    ref.mask_outliers <- intFun.grid.outliers.na(ref.mean, outlier.neg, outlier.pos)
    ref.mask_outliers <- ref.mask_outliers - ref.mask_outliers + 1
    ref <- ref * ref.mask_outliers
    names(ref) <- ref.names
    
    # Compute zonal mean
    
    data.list <- list(mod, ref)
    names <- c("mod", "ref")
    for (i in 1:length(data.list)) {
        data <- raster::stack(data.list[i])
        data.mean <- raster::mean(data, na.rm = TRUE)  # time mean
        # data.sd <- raster::calc(data, sd) # time standard deviation min and max cover the interquartile range
        data.min <- raster::calc(data, fun = function(x) {
            stats::quantile(x, probs = c(0.25), na.rm = TRUE)
        })
        data.max <- raster::calc(data, fun = function(x) {
            stats::quantile(x, probs = c(0.75), na.rm = TRUE)
        })
        
        # get 'z'
        z <- data.mean
        xy <- sp::coordinates(data.mean)
        lat <- xy[, 2]
        z[] <- lat
        # compute zonal mean values
        zonal.mean.mean <- raster::zonal(data.mean, z, "mean", digits = 5)
        zonal.mean.min <- raster::zonal(data.min, z, "mean", digits = 5)
        zonal.mean.max <- raster::zonal(data.max, z, "mean", digits = 5)
        zone <- zonal.mean.mean[, 1]
        # combine values
        rbind(zonal.mean.min, zonal.mean.max)
        zonal.mean <- data.frame(zone, zonal.mean.min[, 2], zonal.mean.max[, 2], zonal.mean.mean[, 2])
        zonal.mean <- stats::na.omit(zonal.mean)  # omit rows with NA
        assign(paste("zonal.mean", names[i], sep = "."), zonal.mean)
    }
    
    # prepare plot file name
    
    my.filename <- paste(variable.name, ref.id, "zonalMean", sep = "-")
    # plot title
    my.title <- paste("Zonal mean", variable.name, mod.id, "vs", ref.id, "from", start.date, "to", end.date)
    # colors
    my.col.mod <- "black"
    my.col.ref <- "red"
    my.col.mod.range <- grDevices::adjustcolor(my.col.mod, alpha = 0.25)
    my.col.ref.range <- grDevices::adjustcolor(my.col.ref, alpha = 0.25)
    
    # limits
    my.xlim <- c(min(zonal.mean.mod[1], zonal.mean.ref[1]), max(zonal.mean.mod[1], zonal.mean.ref[1]))
    # my.xllim <- c(-90,90)
    my.ylim <- c(min(zonal.mean.mod[2], zonal.mean.ref[2]), max(zonal.mean.mod[3], zonal.mean.ref[3]))
    
    # legend bar text
    legend.bar.text <- latex2exp::TeX(variable.unit)
    
    # model data polygons for uncertainty range
    zonal.mean <- zonal.mean.mod
    zone <- zonal.mean$zone
    zonal.mean.min <- zonal.mean$zonal.mean.min
    zonal.mean.max <- zonal.mean$zonal.mean.max
    zonal.mean.mean <- zonal.mean$zonal.mean.mean
    poly.x.mod <- c(zone, rev(zone))
    poly.y.mod <- c(zonal.mean.min, rev(zonal.mean.max))
    zone.mod <- zone
    
    # reference data polygons for uncertainty range
    zonal.mean <- zonal.mean.ref
    zone <- zonal.mean$zone
    zonal.mean.min <- zonal.mean$zonal.mean.min
    zonal.mean.max <- zonal.mean$zonal.mean.max
    zonal.mean.mean <- zonal.mean$zonal.mean.mean
    poly.x.ref <- c(zone, rev(zone))
    poly.y.ref <- c(zonal.mean.min, rev(zonal.mean.max))
    zone.ref <- zone
    
    # plot
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, ".pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(font.main = 1, mar = c(4, 5, 3, 1), lwd = 1, cex = 1, tcl = 0.3)
    # plot
    graphics::plot(zonal.mean.ref$zone, zonal.mean.ref$zonal.mean.max, main = paste(long.name, my.title, sep = "\n"), type = "l", 
        xlab = "degrees latitude", ylab = legend.bar.text, xlim = my.xlim, ylim = my.ylim, col = NA, las = 1)
    # ref
    graphics::polygon(poly.x.ref, poly.y.ref, col = my.col.ref.range, border = NA)
    graphics::lines(zone.ref, zonal.mean.ref$zonal.mean.mean, col = my.col.ref, lwd = 2)
    # mod
    graphics::polygon(poly.x.mod, poly.y.mod, col = my.col.mod.range, border = NA)
    graphics::lines(zone.mod, zonal.mean.mod$zonal.mean.mean, col = my.col.mod, lwd = 2)
    # legend
    graphics::legend("topright", c("model mean", "model IQR", "reference mean", "reference IQR"), col = c(my.col.mod, my.col.mod.range, 
        my.col.ref, my.col.ref.range), lty = c(1, NA, 1, NA), lwd = 2, pch = c(NA, 15, NA, 15), bty = "n")
    # ticks
    graphics::axis(1, at = seq(-90, 90, 10), labels = FALSE, tcl = 0.3)
    graphics::axis(3, at = seq(-90, 90, 10), labels = FALSE, tcl = 0.3)
    graphics::axis(4, labels = FALSE, tcl = 0.3)
    
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}
if (getRversion() >= "2.15.1") utils::globalVariables(c("zonal.mean.mod", "zonal.mean.ref"))
