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
#' @param outlier.factor A number that is used to define outliers, e.g. 10. Outliers are all values that
#' exceed the interquartile range multiplied by the outlier factor defined here.
#' 180 degrees and FALSE if you want longitudes to range from 0 to 360 degrees
#' @param plot.width Number that gives the plot width, e.g. 8
#' @param plot.height Number that gives the plot height, e.g. 4
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @param myLevel A number that determines what level of the output netCDF file to use.
#' This is relevant for files with multiple levels, which applies to soil data.
#' By default, myLevel is set to 1.
#' @param globalValues Either 'globalMean' or 'globalSum'. If set to 'globalMean',
#' values are averaged across all grid cells. If set to 'globalSum', values are
#' summed up. The sum is weighted by grid cell area.
#' @return Three plots in PDF format that give the (a) monthly time series, (b)
#' annual time series, and (c) the climatological mean seasonal cycle of model
#' and reference data. Concerning the latter, thee bold line presents the mean
#' and the shaded area gives the total range caused by interannual variability.
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
#' unit.conv.ref, variable.unit)
#' @export
seasonalCycle <- function(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod, unit.conv.ref, variable.unit, outlier.factor = 1000,
    plot.width = 6, plot.height = 5, outputDir = FALSE, myLevel = 1, globalValues = "globalMean") {

    # Data preparation
    nc <- ncdf4::nc_open(nc.mod)
    variable.name <- names(nc[["var"]])
    ncdf4::nc_close(nc)
    variable.name <- variable.name[length(variable.name)]  # take the last variable (relevant for CanESM5)
    variable.name <- ifelse(variable.name == "burntFractionAll", "burnt", variable.name)  # rename burntFractionAll to shorter name
    variable.name <- toupper(variable.name)  # make variable name upper-case
    mod <- raster::brick(nc.mod, level = myLevel)
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

    # create a mask to excludes all grid cells that the model and reference data do not have in common.  This mask varies in
    # time.
    mask <- (mod * ref)
    mask <- mask - mask + 1
    mod <- mod * mask
    names(mod) <- mod.names  # this adds the corresponding dates
    ref <- ref * mask
    names(ref) <- ref.names  # this adds the corresponding dates
    # now mod and ref are based on the same grid cells

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
        if (globalValues == "globalMean") {
            time.series <- raster::cellStats(data, stat = "mean", na.rm = TRUE)  # spatial mean
            time.series <- data.frame(date, month, time.series)
            colnames(time.series) <- c("date", "month", "data")
            year <- format(as.Date(time.series$date), "%Y")
            time.series <- data.frame(year, time.series)
            annualValues <- tapply(time.series$data, time.series$year, mean)
            annualValues <- data.frame(annualValues)
            year <- rownames(annualValues)
            year <- as.numeric(year)
            annualValues <- data.frame(year, annualValues)
            time.series <- merge(time.series, annualValues, by = "year")
            time.series[time.series$month != 1, "annualValues"] <- NA
            time.series.yr <- stats::na.omit(time.series)
            time.series.yr <- subset(time.series.yr, select = -c(data, month))
        }

        if (globalValues == "globalSum") {
            area <- raster::area(data) * 1000 * 1000  # area of grid cell in m^2
            dataXarea <- data * area
            time.series <- raster::cellStats(dataXarea, stat = "sum", na.rm = TRUE)  # spatial sum
            time.series <- data.frame(date, month, time.series)
            colnames(time.series) <- c("date", "month", "data")
            year <- format(as.Date(time.series$date), "%Y")
            time.series <- data.frame(year, time.series)
            annualValues <- tapply(time.series$data, time.series$year, sum)
            annualValues <- data.frame(annualValues)
            year <- rownames(annualValues)
            year <- as.numeric(year)
            annualValues <- data.frame(year, annualValues)
            time.series <- merge(time.series, annualValues, by = "year")
            time.series[time.series$month != 1, "annualValues"] <- NA
            time.series.yr <- stats::na.omit(time.series)
            time.series.yr <- subset(time.series.yr, select = -c(data, month))
        }

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
        assign(paste("time.series.yr", names[i], sep = "."), time.series.yr)
        assign(paste("seasonMean", names[i], sep = "."), seasonMean)
        assign(paste("poly.y", names[i], sep = "."), poly.y)

    }

    # prepare plot inputs

    my.filename <- paste(variable.name, ref.id, "seasonalCycle", sep = "-")

    # colors
    myPalette <- "viridis"
    myCol <- viridis::viridis(n = 2, end = 0.75, option = myPalette)
    my.col.mod <- myCol[1]
    my.col.ref <- myCol[2]
    my.col.mod.range <- grDevices::adjustcolor(my.col.mod, alpha = 0.25)
    my.col.ref.range <- grDevices::adjustcolor(my.col.ref, alpha = 0.25)

    # legend bar text
    legend.bar.text <- latex2exp::TeX(variable.unit)

    # plot

    # (a) monthly time series plot
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, "-timeseries.pdf", sep = ""), width = plot.width * 1.5, height = plot.height)
    }
    graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(4, 5, 3, 2), lwd = 1, cex = 1, tcl = 0.3)

    my.title <- paste(variable.name, mod.id, "vs", ref.id, "from", start.date, "to", end.date)

    my.ylim <- c(min(time.series.mod$data, time.series.ref$data, na.rm = TRUE), max(time.series.mod$data, time.series.ref$data,
        na.rm = TRUE))

    graphics::plot(as.Date(time.series.ref$date), time.series.ref$data, col = my.col.ref, main = paste(long.name, my.title, sep = "\n"),
        type = "l", xlab = NA, ylab = legend.bar.text, ylim = my.ylim, las = 1, cex.main = 1)
    graphics::lines(as.Date(time.series.mod$date), time.series.mod$data, col = my.col.mod)
    graphics::legend("topleft", c("model", "reference"), col = c(my.col.mod, my.col.ref), pch = 16, bty = "n")
    graphics::axis(4, labels = FALSE, tcl = 0.3)

    if (outputDir != FALSE) {
        grDevices::dev.off()
    }

    # (b) annual time series plot

    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, "-timeseries-YLY.pdf", sep = ""), width = plot.width, height = plot.height)
    }

    graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(4, 5, 3, 2), lwd = 1, cex = 1, tcl = 0.3)

    my.title <- paste(variable.name, mod.id, "vs", ref.id, "from", start.date, "to", end.date)

    my.ylim <- c(min(time.series.yr.mod$annualValues, time.series.yr.ref$annualValues, na.rm = TRUE), max(time.series.yr.mod$annualValues,
        time.series.yr.ref$annualValues, na.rm = TRUE))

    graphics::plot(as.Date(time.series.yr.ref$date), time.series.yr.ref$annualValues, col = my.col.ref, main = paste(paste(long.name,
        "(annual)", sep = " "), my.title, sep = "\n"), type = "l", xlab = NA, ylab = legend.bar.text, ylim = my.ylim, las = 1,
        cex.main = 0.8)
    graphics::lines(as.Date(time.series.yr.mod$date), time.series.yr.mod$annualValues, col = my.col.mod)
    graphics::legend("topleft", c("model", "reference"), col = c(my.col.mod, my.col.ref), pch = 16, bty = "n")
    graphics::axis(4, labels = FALSE, tcl = 0.3)

    if (outputDir != FALSE) {
        grDevices::dev.off()
    }

    # (c) climatological mean seasonal cycle
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, "-clim.pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(4, 5, 3, 2), lwd = 1, cex = 1, tcl = 0.3)

    my.title <- paste("Monthly", variable.name, mod.id, "vs", ref.id, "from", start.date, "to", end.date)

    my.ylim <- c(min(poly.y.mod, poly.y.ref), max(poly.y.mod, poly.y.ref))

    graphics::plot(seq(1, 12, 1), seasonMean.mod, main = paste(long.name, my.title, sep = "\n"), type = "l", xaxt = "n", xlab = "",
        ylab = legend.bar.text, ylim = my.ylim, col = NA, las = 1, cex.main = 0.8)

    months <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
    axis(1, at = 1:12, labels = months)
    axis(3, at = 1:12, labels = FALSE)

    # ref
    months <- seq(1, 12, 1)
    graphics::polygon(poly.x, poly.y.ref, col = my.col.ref.range, border = NA)
    graphics::lines(months, seasonMean.ref, col = my.col.ref, lwd = 2)

    # mod
    graphics::polygon(poly.x, poly.y.mod, col = my.col.mod.range, border = NA)
    graphics::lines(months, seasonMean.mod, col = my.col.mod, lwd = 2)

    # legend
    graphics::legend("topleft", c("model", "reference"), col = c(my.col.mod, my.col.ref), pch = 16, bty = "n")
    # ticks
    graphics::axis(4, labels = FALSE, tcl = 0.3)

    if (outputDir != FALSE) {
        grDevices::dev.off()
    }

}
if (getRversion() >= "2.15.1") utils::globalVariables(c("time.series.mod", "time.series.ref", "time.series.yr.mod", "time.series.yr.ref",
    "poly.y.mod", "poly.y.ref", "seasonMean.mod", "seasonMean.ref", "start.date.mod", "start.date.ref"))

