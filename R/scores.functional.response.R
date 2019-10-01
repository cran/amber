################################################################################
#' Response of a variable to its forcing
#' @description This function conducts a relationship analysis by assessing the functional
#' response of two variables. The variables consist of a dependent variable \eqn{y}
#' (e.g. GPP) and an independent variable \eqn{x} (e.g. temperature).
#' Usually, there are 4 datasets involved:
#' (1) modeled \eqn{x}, (2) modeled \eqn{y}, (3) reference \eqn{x}, (4) reference \eqn{y},
#' where \eqn{reference} refers to observation-based datasets. The time period of
#' analysis should cover the period that all data sets have in common.
#' When the model is forced with observations, (1) and (3) are identical.
#' The common time period is then defined by (2) and (4).
#' The performance of a model is expressed through a score value that ranges
#' from zero to one, where increasing values imply better performance:
#'
#' \eqn{$\varepsilon_{func}^u=\sqrt{\frac{\int(f_{mod}(u)-f_{ref}(u))^2du}{\int f_{ref}(u)^2du}}$}
#'
#' where \eqn{f_{mod}(u)} and \eqn{f_{ref}(u)} are the binned time means of
#' the model and reference data, respectively.
#'
#' @param nc.mod.x A string that gives the path and name of the netcdf file that contains the model forcing data, e.g. '/home/model_tas.nc'
#' @param nc.mod.y A string that gives the path and name of the netcdf file that contains the model output data, e.g. '/home/model_gpp.nc'
#' @param nc.ref.x A string that gives the path and name of the netcdf file that contains the reference forcing data, e.g. '/home/ref_tas.nc'
#' @param nc.ref.y A string that gives the path and name of the netcdf file that contains the reference data, e.g. '/home/ref_gpp.nc'
#' @param mod.id A string that identifies the source of the reference data set, e.g. 'CLASSIC'
#' @param ref.id A string that identifies the source of the reference data set, e.g. 'MODIS'
#' @param unit.conv.mod.x A number that is used as a factor to convert the unit of the model forcing data, e.g. 86400
#' @param unit.conv.mod.y A number that is used as a factor to convert the unit of the model output data, e.g. 86400
#' @param unit.conv.ref.x A number that is used as a factor to convert the unit of the reference forcing data, e.g. 86400
#' @param unit.conv.ref.y A number that is used as a factor to convert the unit of the reference data, e.g. 86400
#' @param x.bin A number that gives the size of the x-bin, e.g. 2
#' @param y.bin A number that gives the size of the y-bin, e.g. 2
#' @param my.xlab A string that serves as the label for the x-axis, e.g. expression(paste('near surface air temperature (',degree,'C)')) # (R notation)
#' @param my.ylab A string that serves as the label for the y-axis, e.g. expression('GPP (g C m'^{-2}~'day'^{-1}~')') # (R notation)
#' @param legendA A string that defines the legend location for the upper subplot. The default value is 'topleft'.
#' @param legendB A string that defines the legend location for the lower subplot. The default value is 'topleft'.
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @return A figure in PDF format that shows the functional relation
#' ship between two variables, the frequency of grid cells for each bin,
#' and the score value.
#' @examples
#'
#' library(amber)
#' library(ncdf4)
#' library(raster)
#' nc.mod.x <- system.file('extdata/modelRegular', 'tas_monthly.nc', package = 'amber')
#' nc.mod.y <- system.file('extdata/modelRegular', 'gpp_monthly.nc', package = 'amber')
#' nc.ref.x <- system.file('extdata/modelRegular', 'tas_monthly.nc', package = 'amber')
#' nc.ref.y <- system.file('extdata/referenceRegular', 'gpp_GBAF_128x64.nc', package = 'amber')
#' mod.id <- 'CLASSIC' # define a model experiment ID
#' ref.id <- 'GBAF' # give reference dataset a name
#' unit.conv.mod.x <- 1
#' unit.conv.mod.y <- 86400*1000
#' unit.conv.ref.x <- 1
#' unit.conv.ref.y <- 86400*1000
#' x.bin <- 2.5 # adjust as required
#' y.bin <- 1 # adjust as required
#' my.xlab <- expression(paste('near surface air temperature (',degree,'C)')) # (R notation)
#' my.ylab <- expression('GPP (g C m'^{-2}~'day'^{-1}~')') # (R notation)
#'
#' scores.functional.response(nc.mod.x, nc.mod.y, nc.ref.x, nc.ref.y, mod.id,
#' ref.id, unit.conv.mod.x, unit.conv.mod.y, unit.conv.ref.x, unit.conv.ref.y,
#' x.bin, y.bin, my.xlab, my.ylab)
#'
#' @export
scores.functional.response <- function(nc.mod.x, nc.mod.y, nc.ref.x, nc.ref.y, mod.id, ref.id, unit.conv.mod.x, unit.conv.mod.y, unit.conv.ref.x,
    unit.conv.ref.y, x.bin, y.bin, my.xlab, my.ylab, legendA = "topleft", legendB = "topleft", outputDir = FALSE) {
    # (1) Find a common time period and subset the data accordingly

    # variable names (x)
    nc <- ncdf4::nc_open(nc.mod.x)
    variable.name <- names(nc[["var"]])
    variable.name <- variable.name[length(variable.name)]  # take the last variable (relevant for CanESM5)
    variable.name.x <- variable.name
    variable.name.x <- toupper(variable.name.x)
    ncdf4::nc_close(nc)

    # variable names (y)
    nc <- ncdf4::nc_open(nc.mod.y)
    variable.name <- names(nc[["var"]])
    variable.name <- variable.name[length(variable.name)]  # take the last variable (relevant for CanESM5)
    variable.name.y <- variable.name
    variable.name.y <- toupper(variable.name.y)  # make variable name upper-case
    ncdf4::nc_close(nc)
    # model data
    mod.x <- raster::brick(nc.mod.x)
    mod.y <- raster::brick(nc.mod.y)
    dates.mod <- raster::getZ(mod.y)
    dates.mod <- format(as.Date(dates.mod), "%Y-%m")  # only year and month
    start.date.mod <- min(dates.mod)
    end.date.mod <- max(dates.mod)
    # reference data
    ref.x <- raster::brick(nc.ref.x)
    ref.y <- raster::brick(nc.ref.y)
    dates.ref <- raster::getZ(ref.y)
    dates.ref <- format(as.Date(dates.ref), "%Y-%m")  # only year and month
    start.date.ref <- min(dates.ref)
    end.date.ref <- max(dates.ref)

    # find common time period
    start.date <- max(start.date.mod, start.date.ref)
    end.date <- min(end.date.mod, end.date.ref)

    # subset common time period
    mod.x <- mod.x[[which(format(as.Date(raster::getZ(mod.x)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(mod.x)), "%Y-%m") <=
        end.date)]]
    mod.y <- mod.y[[which(format(as.Date(raster::getZ(mod.y)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(mod.y)), "%Y-%m") <=
        end.date)]]
    ref.x <- ref.x[[which(format(as.Date(raster::getZ(ref.x)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(ref.x)), "%Y-%m") <=
        end.date)]]
    ref.y <- ref.y[[which(format(as.Date(raster::getZ(ref.y)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(ref.y)), "%Y-%m") <=
        end.date)]]

    # unit conversion if appropriate
    mod.x <- mod.x * unit.conv.mod.x
    mod.y <- mod.y * unit.conv.mod.y
    ref.x <- ref.x * unit.conv.ref.x
    ref.y <- ref.y * unit.conv.ref.y

    my.filename <- paste(variable.name.y, variable.name.x, ref.id, "functional_response", sep = "_")

    # (2) Convert data into bins model data
    mod.x.mean <- raster::mean(mod.x, na.rm = TRUE)
    mod.y.mean <- raster::mean(mod.y, na.rm = TRUE)
    ref.x.mean <- raster::mean(ref.x, na.rm = TRUE)
    ref.y.mean <- raster::mean(ref.y, na.rm = TRUE)

    # make sure all data is based on the same grid cells
    mask <- (mod.x.mean * mod.y.mean * ref.x.mean * ref.y.mean)
    mask <- mask - mask + 1
    mod.x.mean <- mod.x.mean/mask
    mod.y.mean <- mod.y.mean/mask
    ref.x.mean <- ref.x.mean/mask
    ref.y.mean <- ref.y.mean/mask

    # convert to bins
    mod.x.bin <- intFun.bin(mod.x.mean, x.bin)
    mod.y.bin <- intFun.bin(mod.y.mean, y.bin)
    ref.x.bin <- intFun.bin(ref.x.mean, x.bin)
    ref.y.bin <- intFun.bin(ref.y.mean, y.bin)

    # (3) calculate the mean value and corresponding variability of gpp for each temperature interval model
    x <- raster::values(mod.x.bin)
    y <- raster::values(mod.y.bin)
    relation.mod <- data.frame(x, y)
    relation.mod <- stats::na.omit(relation.mod)  # omit rows with NA
    index <- list(relation.mod$x)
    y.mean <- data.frame(tapply(relation.mod$y, index, mean))  # mean value for x bin
    y.sd <- data.frame(tapply(relation.mod$y, index, sd))
    y.sd[is.na(y.sd)] = 0  # NA standard deviation (sd) is set to zero. This is necessary for plotting sd polygons
    my.table <- table(relation.mod$x, relation.mod$y)
    freq.x.mod <- margin.table(my.table, 1)  # number of grid cells per x bin
    freq.x.mod <- as.data.frame(freq.x.mod)
    colnames(freq.x.mod) <- c(variable.name.x, "n")
    relation.mod <- cbind(y.mean, y.sd)
    x <- as.numeric(rownames(relation.mod))
    relation.mod <- cbind(x, relation.mod)  # This adds the x-axis values as a column
    colnames(relation.mod) <- c(variable.name.x, variable.name.y, "sd")
    relation.mod <- merge(relation.mod, freq.x.mod, by = variable.name.x)
    rownames(relation.mod) <- c()

    # reference
    x <- raster::values(ref.x.bin)
    y <- raster::values(ref.y.bin)
    relation.ref <- data.frame(x, y)
    relation.ref <- stats::na.omit(relation.ref)  # omit rows with NA
    index <- list(relation.ref$x)
    y.mean <- data.frame(tapply(relation.ref$y, index, mean))  # mean value for x bin
    y.sd <- data.frame(tapply(relation.ref$y, index, sd))
    y.sd[is.na(y.sd)] = 0  # NA standard deviation (sd) is set to zero. This is necessary for plotting sd polygons
    my.table <- table(relation.ref$x, relation.ref$y)
    freq.x.ref <- margin.table(my.table, 1)  # number of grid cells per x bin
    freq.x.ref <- as.data.frame(freq.x.ref)
    colnames(freq.x.ref) <- c(variable.name.x, "n")
    relation.ref <- cbind(y.mean, y.sd)
    x <- as.numeric(rownames(relation.ref))
    relation.ref <- cbind(x, relation.ref)  # This adds the x-axis values as a column
    colnames(relation.ref) <- c(variable.name.x, variable.name.y, "sd")
    relation.ref <- merge(relation.ref, freq.x.ref, by = variable.name.x)
    rownames(relation.ref) <- c()
    relation.all <- merge(relation.mod, relation.ref, by = variable.name.x, all = TRUE)
    relation <- stats::na.omit(relation.all)  # omit rows with NA

    # compute the score
    a <- relation[, 2]
    b <- relation[, 4]
    epsilon.f <- intFun.rel.error(a, b)
    S_funresp <- exp(-epsilon.f)

    # colors
    my.col.mod <- "black"
    my.col.ref <- "red"
    my.col.mod.range <- grDevices::adjustcolor(my.col.mod, alpha = 0.25)
    my.col.ref.range <- grDevices::adjustcolor(my.col.ref, alpha = 0.25)
    # polygon coordinates
    poly.x.mod <- c(relation.mod[, 1], rev(relation.mod[, 1]))
    poly.y.mod <- c(relation.mod[, 2] + relation.mod[, 3], rev(relation.mod[, 2] - relation.mod[, 3]))
    poly.x.ref <- c(relation.ref[, 1], rev(relation.ref[, 1]))
    poly.y.ref <- c(relation.ref[, 2] + relation.ref[, 3], rev(relation.ref[, 2] - relation.ref[, 3]))

    # limits
    my.xlim <- c(min(relation.ref[, 1]), max(relation.ref[, 1]))
    my.ylim <- c(min(poly.y.mod, poly.y.ref), max(poly.y.mod, poly.y.ref))

    # plot
    oldpar <- graphics::par(mfrow = c(1,2))
    on.exit(graphics::par(oldpar))
    my.filename <- gsub("_", "-", my.filename)
    my.filename <- gsub(".", "-", my.filename, fixed = TRUE)
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, ".pdf", sep = ""), width = 3.5, height = 4)
    }
    graphics::par(font.main = 1, mar = c(0, 5, 0, 0), lwd = 1, cex = 1, cex.axis = 1.25, oma = c(5, 1, 3, 1))

    my.matrix <- matrix(c(1, 1, 2, 1, 1, 2), ncol = 2)

    # defines location of subplots
    graphics::layout(my.matrix)  # defines location of subplots
    plot(relation.mod[, 1], relation.mod[, 2], xlim = my.xlim, ylim = my.ylim, ylab = my.ylab, col = NA, las = 1, axes = FALSE)
    graphics::polygon(poly.x.mod, poly.y.mod, col = my.col.mod.range, border = NA)
    graphics::polygon(poly.x.ref, poly.y.ref, col = my.col.ref.range, border = NA)
    graphics::lines(relation.ref[, 1], relation.ref[, 2], col = my.col.ref, lwd = 2)
    graphics::lines(relation.mod[, 1], relation.mod[, 2], col = my.col.mod, lwd = 2)
    graphics::points(relation.ref[, 1], relation.ref[, 2], col = my.col.ref, pch = 17)
    graphics::points(relation.mod[, 1], relation.mod[, 2], col = my.col.mod, pch = 16)
    graphics::box()
    graphics::legend(legendA, c(expression("model mean", paste("model ", sigma, sep = ""), "reference mean", paste("reference ", sigma,
        sep = ""))), col = c(my.col.mod, my.col.mod.range, my.col.ref, my.col.ref.range), lty = c(1, NA, 1, NA), lwd = 2, pch = c(16,
        15, 17, 15), bty = "n")
    # my.legend <- substitute(paste('score = ', score), list(score=round(S_funresp,2))) legend('bottomright', legend=my.legend,
    # bty='n') Ticks
    graphics::axis(1, labels = FALSE, tcl = 0.3)
    graphics::axis(2, labels = TRUE, tcl = 0.3, las = 2)
    graphics::axis(3, labels = FALSE, tcl = 0.3)
    graphics::axis(4, labels = FALSE, tcl = 0.3)
    graphics::mtext(paste(mod.id, " vs ", ref.id, "(", start.date, " to ", end.date, ")", sep = ""), side = 3, line = 1, cex = 0.75)
    # grid cell frequency
    my.ylim <- c(0, max(relation.mod[, 4], relation.ref[, 4]))

    plot(x = relation.mod[, 1] - x.bin/8, y = relation.mod[, 4], xlim = my.xlim, ylim = my.ylim, type = "h", col = my.col.mod, xlab = my.xlab,
        ylab = "grid cell frequency", lwd = 2, las = 1, axes = FALSE)
    graphics::lines(x = relation.ref[, 1] + x.bin/8, y = relation.ref[, 4], ylim = my.ylim, type = "h", col = my.col.ref, xlab = NA,
        ylab = NA, lwd = 2)
    graphics::legend(legendB, lwd = 2, col = c(my.col.mod, my.col.ref), c("model", "reference"), bty = "n")
    graphics::box()
    graphics::mtext(my.xlab, side = 1, line = 3, cex = 0.75)

    graphics::axis(1, labels = TRUE, tcl = 0.3)
    graphics::axis(2, labels = TRUE, tcl = 0.3, las = 2)
    graphics::axis(4, labels = FALSE, tcl = 0.3)

    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}
