################################################################################
#' Zonal mean plots of model and reference data
#' @description This function plots the mean seasonal cycle and corresponding inter-quartile
#' range of model and reference data.
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
#' 180 degrees and FALSE if you want longitudes to range from 0 to 360 degrees
#' @param plot.width Number that gives the plot width, e.g. 8
#' @param plot.height Number that gives the plot height, e.g. 4
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#'
#' @return A plot in PDF format that gives the (a) global mean time series and
#' (b) the climatological mean seasonal cycle of model and reference data.
#' The bold line presents the spatial mean values and the shaded area gives the
#' total range caused by interannual variability.
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
#' seasonalCycle(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit, outlier.factor)
#'
#' @export
seasonalCycle <- function(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod, unit.conv.ref, variable.unit, outlier.factor = 1000,
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
    mod <- mod[[which(format(as.Date(raster::getZ(mod)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(mod)), "%Y-%m") <= end.date)]]
    ref <- ref[[which(format(as.Date(raster::getZ(ref)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(ref)), "%Y-%m") <= end.date)]]

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

    # Compute seasonal cycle

    data.list <- list(mod, ref)
    names <- c("mod", "ref")

    for (i in 1:length(data.list)) {
        data <- raster::stack(data.list[i])

        # date
        date <- as.Date(names(data), format = "X%Y.%m.%d")
        index <- format(as.Date(names(mod), format = "X%Y.%m.%d"), format = "%m")
        month <- as.numeric(index)

        # time series
        time.series <- raster::cellStats(data, stat = "mean", na.rm = TRUE)
        time.series <- data.frame(date, month, time.series)
        colnames(time.series) <- c("date", "month", "data")

        # seasonal cycle
        index <- list(time.series$month)
        seasonMean <- tapply(time.series$data, index, mean, na.rm = TRUE)
        seasonMin <- tapply(time.series$data, index, min, na.rm = TRUE)
        seasonMax <- tapply(time.series$data, index, max, na.rm = TRUE)

        months <- seq(1, 12, 1)

        # make polygon inputs for plot
        poly.x <- c(months, rev(months))
        poly.y <- c(seasonMin, rev(seasonMax))

        # assign names
        assign(paste("time.series", names[i], sep = "."), time.series)
        assign(paste("seasonMean", names[i], sep = "."), seasonMean)
        assign(paste("poly.y", names[i], sep = "."), poly.y)
    }

    # prepare plot inputs

    my.filename <- paste(variable.name, ref.id, "seasonalCycle", sep = "-")

    # colors
    my.col.mod <- "black"
    my.col.ref <- "red"
    my.col.mod.range <- grDevices::adjustcolor(my.col.mod, alpha = 0.25)
    my.col.ref.range <- grDevices::adjustcolor(my.col.ref, alpha = 0.25)

    # legend bar text
    legend.bar.text <- latex2exp::TeX(variable.unit)

    # plot
    oldpar <- graphics::par(mfrow = c(1,2))
    on.exit(graphics::par(oldpar))
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, "-timeseries.pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(4, 5, 3, 1), lwd = 1, cex = 1, tcl = 0.3)

    # (a) time series plot

    # plot title

    my.title <- paste("Global mean", variable.name, mod.id, "vs", ref.id, "from", start.date, "to", end.date)

    # time series

    my.ylim <- c(min(time.series.mod$data, time.series.ref$data, na.rm = TRUE), max(time.series.mod$data, time.series.ref$data, na.rm = TRUE))

    graphics::plot(as.Date(time.series.ref$date), time.series.ref$data, col = "red", main = paste(long.name, my.title, sep = "\n"),
        type = "l", xlab = NA, ylab = legend.bar.text, ylim = my.ylim, las = 1)
    graphics::lines(as.Date(time.series.mod$date), time.series.mod$data, col = "black")
    graphics::legend("topleft", c("model", "reference"), col = c("black", "red"), lty = 1, bty = "n")
    # ticks
    graphics::axis(4, labels = FALSE, tcl = 0.3)
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }

    # (b) climatological mean seasonal cycle
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, "-clim.pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(4, 5, 3, 1), lwd = 1, cex = 1, tcl = 0.3)

    my.ylim <- c(min(poly.y.mod, poly.y.ref), max(poly.y.mod, poly.y.ref))

    graphics::plot(seq(1, 12, 1), seasonMean.mod, main = paste(long.name, my.title, sep = "\n"), type = "l", xaxt = "n", xlab = "",
        ylab = legend.bar.text, ylim = my.ylim, col = NA, las = 1)

    months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    axis(1, at = 1:12, labels = months)

    # ref
    months <- seq(1, 12, 1)
    graphics::polygon(poly.x, poly.y.ref, col = my.col.ref.range, border = NA)
    graphics::lines(months, seasonMean.ref, col = my.col.ref, lwd = 2)

    # mod
    graphics::polygon(poly.x, poly.y.mod, col = my.col.mod.range, border = NA)
    graphics::lines(months, seasonMean.mod, col = my.col.mod, lwd = 2)

    # legend
    graphics::legend("topleft", c("model (mean and total range)", "reference (mean and total range)"), col = c(my.col.mod, my.col.ref),
        lty = 1, lwd = 2, bty = "n")
    # ticks
    graphics::axis(4, labels = FALSE, tcl = 0.3)

    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}
if (getRversion() >= "2.15.1") utils::globalVariables(c("time.series.mod", "time.series.ref", "time.series.mod", "poly.y.mod", "poly.y.ref",
    "seasonMean.mod", "seasonMean.ref", "start.date.mod", "start.date.ref"))

