################################################################################ 
#' Plots the time-mean of a variable stored in NetCDF model output on a regular grid
#' @description This function plots the time-mean, spatial-mean, zonal mean, and
#' seasonal cycle of variable stored in NetCDF model output. The function expects
#' model data to be on a regular grid.
#' @param long.name A string that gives the full name of the variable, e.g. 'Gross primary productivity'
#' @param nc.mod A string that gives the path and name of the netcdf file that contains the model output, e.g. '/home/model_gpp.nc'
#' @param mod.id A string that identifies the source of the reference data set, e.g. 'CLASSIC'
#' @param unit.conv.mod A number that is used as a factor to convert the unit of the model data, e.g. 86400
#' @param variable.unit A string that gives the final units using LaTeX notation, e.g. 'gC m$^{-2}$ day$^{-1}$'
#' @param timePeriod A string that gies the time period over which to average the data, e.g. c('1980-01', '2017-12')
#' @param outlier.factor A number that is used to define outliers, e.g. 10.
#'  Plotting raster objects that contain extreme outliers lead to figures where
#'  most grid cells are presented by a single color since the color legend covers
#'  the entire range of values. To avoid this, the user may define outliers that
#'  will be masked out and marked with a red dot. Outliers are all values that
#'  exceed the interquartile range multiplied by the outlier factor defined here.
#' @param my.xlim An R object that gives the longitude range that you wish to plot, e.g. c(-180, 180)
#' @param my.ylim An R object that gives the longitude range that you wish to plot, e.g. c(-90, 90)
#' @param plot.width Number that gives the plot width, e.g. 8
#' @param plot.height Number that gives the plot height, e.g. 4
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory
#'
#' @return Figures in PDF format that show the time-mean, spatial-mean, zonal mean, and seasonal cycle.
#'
#' @examples
#'
#' library(amber)
#' library(ncdf4)
#' library(raster)
#'
#' long.name <- 'Gross primary productivity'
#' nc.mod <- system.file('extdata/modelRegular', 'gpp_monthly.nc', package = 'amber')
#' mod.id <- 'CLASSIC' # define a model experiment ID
#' unit.conv.mod <- 86400*1000 # optional unit conversion for model data
#' variable.unit <- 'gC m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
#' timePeriod <- c('1980-01', '2017-12')
#' outlier.factor <- 1
#'
#' plotNc(long.name, nc.mod, mod.id, unit.conv.mod, variable.unit, timePeriod, outlier.factor)
#'
#' @export
plotNc <- function(long.name, nc.mod, mod.id, unit.conv.mod, variable.unit, timePeriod, outlier.factor = 1000, my.xlim = c(-180, 
    180), my.ylim = c(-60, 85), plot.width = 8, plot.height = 3.8, outputDir = FALSE) {
    
    shp.filename <- system.file("extdata/ne_110m_land/ne_110m_land.shp", package = "amber")
    land <- raster::shapefile(shp.filename)
    
    nc <- ncdf4::nc_open(nc.mod)
    variable.name <- names(nc[["var"]])
    ncdf4::nc_close(nc)
    
    variable.name <- variable.name[length(variable.name)]  # take the last variable (relevant for CanESM5)
    variable.name <- ifelse(variable.name == "burntFractionAll", "burnt", variable.name)  # rename burntFractionAll to shorter name
    variable.name <- toupper(variable.name)  # make variable name upper-case
    mod <- raster::brick(nc.mod)
    mod <- raster::rotate(mod)
    
    dates.mod <- raster::getZ(mod)
    dates.mod <- format(as.Date(dates.mod), "%Y-%m")  # only year and month
    start.date.mod <- min(dates.mod)
    end.date.mod <- max(dates.mod)
    
    # find common time period
    start.date <- max(start.date.mod, timePeriod[1])
    end.date <- min(end.date.mod, timePeriod[2])
    
    # subset time period
    mod <- mod[[which(format(as.Date(raster::getZ(mod)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(mod)), "%Y-%m") <= 
        end.date)]]
    # get layer names
    mod.names <- names(mod)
    
    # unit conversion if appropriate
    mod <- mod * unit.conv.mod
    
    # Make a string that summarizes metadata
    meta.data.mod <- paste(variable.name, mod.id, "from", start.date, "to", end.date, sep = "_")
    
    # all extreme outliers are set to NA in the grid they can be marked as a dot in the plot model data
    mod.mean <- raster::mean(mod, na.rm = TRUE)  # time mean
    mod.outlier_range <- intFun.grid.define.outlier(mod.mean, outlier.factor)  # define outlier range
    outlier.neg <- mod.outlier_range[1]
    outlier.pos <- mod.outlier_range[2]
    mod.mask_outliers <- intFun.grid.outliers.na(mod.mean, outlier.neg, outlier.pos)
    mod.mask_outliers <- mod.mask_outliers - mod.mask_outliers + 1
    mod <- mod * mod.mask_outliers
    names(mod) <- mod.names
    mod.outlier.points <- intFun.grid.outliers.points(mod.mean, outlier.neg, outlier.pos)
    
    #-----------------------------------------------------------------------------
    
    # (1) Time mean value
    
    #-----------------------------------------------------------------------------
    
    mod.mean <- raster::mean(mod, na.rm = TRUE)  # time mean
    mmi.mean <- intFun.min.max.int.mod.ref(mod.mean, mod.mean)
    raster::metadata(mod.mean) <- list(paste(variable.name, mod.id, "mod_mean", sep = "_"), paste("Mean", meta.data.mod, 
        sep = "_"), mmi.mean, variable.unit)
    data <- mod.mean
    
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
    
    my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
    my.legend.args <- list(text = legend.bar.text, side = 2, font = 1, line = 1, cex = 1)
    
    # plot name
    my.filename <- gsub("_", "-", my.filename)
    my.filename <- gsub(".", "-", my.filename, fixed = TRUE)
    
    # plot
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, "-plotNc.pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(3, 3, 3, 4), lwd = 1, cex = 1)
    raster::plot(data, col = NA, legend = FALSE, xlim = my.xlim, ylim = my.ylim, main = paste(long.name, my.title, sep = "\n"), 
        axes = FALSE)
    raster::plot(land, col = "grey", border = NA, add = TRUE)
    raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, add = TRUE)
    graphics::axis(1, labels = TRUE, tcl = 0.3)
    graphics::axis(2, labels = TRUE, tcl = 0.3, las = 2)
    # mark grid cells with extreme outliers
    graphics::points(mod.outlier.points, pch = 16, cex = 1, col = "red")
    plot(land, add = TRUE)
    # legend
    plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args, 
        legend.width = 1.5, legend.shrink = 1, font = 1)
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
    
    #-----------------------------------------------------------------------------
    
    # Zonal mean
    
    #-----------------------------------------------------------------------------
    
    data <- mod
    data.mean <- raster::mean(data, na.rm = TRUE)  # time mean
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
    # combine values rbind(zonal.mean.min, zonal.mean.max)
    zonal.mean <- data.frame(zone, zonal.mean.min[, 2], zonal.mean.max[, 2], zonal.mean.mean[, 2])
    zonal.mean <- stats::na.omit(zonal.mean)  # omit rows with NA
    colnames(zonal.mean) <- c("lat", "q25", "q75", "mean")
    zonal.mean.mod <- zonal.mean
    
    # prepare plot file name
    
    my.filename <- paste(variable.name, mod.id, "zonalMean-plotNc", sep = "-")
    # plot title
    my.title <- paste("Zonal mean", variable.name, mod.id, "from", start.date, "to", end.date)
    # colors
    my.col.mod <- "black"
    my.col.mod.range <- grDevices::adjustcolor(my.col.mod, alpha = 0.25)
    
    # limits
    my.xlim <- c(min(zonal.mean.mod[1]), max(zonal.mean.mod[1]))
    # my.xllim <- c(-90,90)
    my.ylim <- c(min(zonal.mean.mod[2]), max(zonal.mean.mod[3]))
    
    # legend bar text
    legend.bar.text <- latex2exp::TeX(variable.unit)
    
    # model data polygons for uncertainty range
    zonal.mean <- zonal.mean.mod
    zone <- zonal.mean$lat
    zonal.mean.min <- zonal.mean$q25
    zonal.mean.max <- zonal.mean$q75
    zonal.mean.mean <- zonal.mean$mean
    poly.x.mod <- c(zone, rev(zone))
    poly.y.mod <- c(zonal.mean.min, rev(zonal.mean.max))
    zone.mod <- zone
    
    # plot
    
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, ".pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(font.main = 1, mar = c(4, 5, 3, 1), lwd = 1, cex = 1, tcl = 0.3)
    # plot
    graphics::plot(zonal.mean$lat, zonal.mean$q75, main = paste(long.name, my.title, sep = "\n"), type = "l", xlab = "degrees latitude", 
        ylab = legend.bar.text, xlim = my.xlim, ylim = my.ylim, col = NA, las = 1)
    # mod
    graphics::polygon(poly.x.mod, poly.y.mod, col = my.col.mod.range, border = NA)
    graphics::lines(zone, zonal.mean.mod$mean, col = my.col.mod, lwd = 2)
    # legend
    graphics::legend("topright", c("model mean and IQR"), col = c(my.col.mod), lty = 1, lwd = 2, bty = "n")
    # ticks
    graphics::axis(1, at = seq(-90, 90, 10), labels = FALSE, tcl = 0.3)
    graphics::axis(3, at = seq(-90, 90, 10), labels = FALSE, tcl = 0.3)
    graphics::axis(4, labels = FALSE, tcl = 0.3)
    
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
    
    #---------------------------------------------------------------------------
    
    # Seasonal cyle
    
    #---------------------------------------------------------------------------
    
    data <- mod
    
    # date
    date <- as.Date(names(data), format = "X%Y.%m.%d")
    index <- format(as.Date(names(mod), format = "X%Y.%m.%d"), format = "%m")
    month <- as.numeric(index)
    
    # time series
    time.series <- raster::cellStats(data, stat = "mean", na.rm = TRUE)
    time.series <- data.frame(date, month, time.series)
    colnames(time.series) <- c("date", "month", "data")
    time.series.mod <- time.series
    
    # seasonal cycle
    index <- list(time.series$month)
    seasonMean <- tapply(time.series$data, index, mean, na.rm = TRUE)
    seasonMin <- tapply(time.series$data, index, min, na.rm = TRUE)
    seasonMax <- tapply(time.series$data, index, max, na.rm = TRUE)
    
    months <- seq(1, 12, 1)
    
    # make polygon inputs for plot
    poly.x <- c(months, rev(months))
    poly.y <- c(seasonMin, rev(seasonMax))
    
    
    seasonMean.mod <- seasonMean
    seasonMin.mod <- seasonMin
    seasonMax.mod <- seasonMax
    
    poly.x.mod <- poly.x
    poly.y.mod <- poly.y
    
    # prepare plot inputs
    
    my.filename <- paste(variable.name, mod.id, "seasonalCycle", sep = "-")
    
    # colors
    my.col.mod <- "black"
    my.col.mod.range <- grDevices::adjustcolor(my.col.mod, alpha = 0.25)
    
    # legend bar text
    legend.bar.text <- latex2exp::TeX(variable.unit)
    
    # plot
    
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, "-timeseries-plotNc.pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(4, 5, 3, 1), lwd = 1, cex = 1, tcl = 0.3)
    
    # (a) time series plot
    
    # plot title
    
    my.title <- paste("Global mean", variable.name, mod.id, "from", start.date, "to", end.date)
    
    # time series
    
    my.ylim <- c(min(time.series.mod$data, na.rm = TRUE), max(time.series.mod$data, na.rm = TRUE))
    
    graphics::plot(as.Date(time.series.mod$date), time.series.mod$data, col = "red", main = paste(long.name, my.title, sep = "\n"), 
        type = "l", xlab = NA, ylab = legend.bar.text, ylim = my.ylim, las = 1)
    graphics::lines(as.Date(time.series.mod$date), time.series.mod$data, col = "black")
    # ticks
    graphics::axis(4, labels = FALSE, tcl = 0.3)
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
    
    # (b) climatological mean seasonal cycle
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, "-clim-plotNc.pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(4, 5, 3, 1), lwd = 1, cex = 1, tcl = 0.3)
    
    my.ylim <- c(min(poly.y.mod), max(poly.y.mod))
    
    graphics::plot(seq(1, 12, 1), seasonMean.mod, main = paste(long.name, my.title, sep = "\n"), type = "l", xaxt = "n", 
        xlab = "", ylab = legend.bar.text, ylim = my.ylim, col = NA, las = 1)
    
    months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    axis(1, at = 1:12, labels = months)
    
    # mod
    graphics::polygon(poly.x.mod, poly.y.mod, col = my.col.mod.range, border = NA)
    graphics::lines(seq(1, 12, 1), seasonMean.mod, col = my.col.mod, lwd = 2)
    
    # legend
    graphics::legend("topleft", "model mean and total range", col = c(my.col.mod), lty = 1, lwd = 2, bty = "n")
    # ticks
    graphics::axis(4, labels = FALSE, tcl = 0.3)
    
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
    
}

