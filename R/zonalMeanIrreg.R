################################################################################ 
#' Zonal mean plots of model and reference data on an irregular grid
#' @description This function plots zonal mean values and corresponding inter-quartile
#' ranges of model and reference data. The function expects data to be on an irregular grid.
#' @param long.name A string that gives the full name of the variable, e.g. 'Gross primary productivity'
#' @param nc.mod A string that gives the path and name of the netcdf file that contains the model output, e.g. '/home/model_gpp.nc'
#' @param nc.ref A string that gives the path and name of the netcdf file that contains the reference data output, e.g. '/home/reference_gpp.nc'
#' @param mod.id A string that identifies the source of the reference data set, e.g. 'CLASSIC'
#' @param ref.id A string that identifies the source of the reference data set, e.g. 'MODIS'
#' @param unit.conv.mod A number that is used as a factor to convert the unit of the model data, e.g. 86400
#' @param unit.conv.ref A number that is used as a factor to convert the unit of the reference data, e.g. 86400
#' @param variable.unit A string that gives the final units using LaTeX notation, e.g. 'gC m$^{-2}$ day$^{-1}$'
#' @param my.projection A string that defines the projection of the irregular grid
#' @param outlier.factor A number that is used to define outliers, e.g. 10.
#'  Plotting raster objects that contain extreme outliers lead to figures where
#'  most grid cells are presented by a single color since the color legend covers
#'  the entire range of values. To avoid this, the user may define outliers that
#'  will be masked out and marked with a red dot. Outliers are all values that
#'  exceed the interquartile range multiplied by the outlier factor defined here.
#' @param plot.width Number that gives the plot width, e.g. 8
#' @param plot.height Number that gives the plot height, e.g. 4
#' @param numCores An integer that defines the number of cores, e.g. 2
#' @param timeInt A string that gives the time interval of the model data, e.g. 'month' or 'year'
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.

#' @return A figure in PDF format that gives the zonal mean values of model and
#' reference data. The bold line presents the mean and the shaded area the
#' corresponding interquartile range (IQR). The IQR presents the inter-annual variability
#' and longitudinal variability, combined.
#'
#' @examples
#' \donttest{
#'
#' library(amber)
#' library(doParallel)
#' library(foreach)
#' library(latex2exp)
#' library(ncdf4)
#' library(parallel)
#' library(raster)
#'
#' long.name <- 'Gross primary productivity'
#' nc.mod <- system.file('extdata/modelRotated', 'gpp_monthly.nc', package = 'amber')
#' nc.ref <- system.file('extdata/referenceRotated', 'gpp_GBAF_rotated.nc', package = 'amber')
#' mod.id <- 'CLASSIC' # define a model experiment ID
#' ref.id <- 'GBAF' # give reference dataset a name
#' unit.conv.mod <- 86400*1000 # optional unit conversion for model data
#' unit.conv.ref <- 86400*1000 # optional unit conversion for reference data
#' variable.unit <- 'gC m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
#' outlier.factor <- 1000
#'
#' zonalMeanIrreg(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit, outlier.factor, numCores = 2)
#' }
#' @export
zonalMeanIrreg <- function(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod, unit.conv.ref, variable.unit, outlier.factor = 1000, 
    my.projection = "+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.", plot.width = 8, plot.height = 3.8, 
    numCores = 2, timeInt = "month", outputDir = FALSE) {
    
    regular <- "+proj=longlat +ellps=WGS84"
    rotated <- my.projection
    # get variable name
    nc <- ncdf4::nc_open(nc.mod)
    variable.name <- base::names(nc[["var"]])
    variable.name <- variable.name[1]
    variable.name <- ifelse(variable.name == "burntFractionAll", "burnt", variable.name)  # rename burntFractionAll to shorter name
    variable.name <- toupper(variable.name)  # make variable name upper-case
    ncdf4::nc_close(nc)
    # get data for both, model and reference data
    modRef <- c(nc.mod, nc.ref)
    for (id in 1:length(modRef)) {
        my.nc <- modRef[id]
        nc <- ncdf4::nc_open(my.nc)
        data <- ncdf4::ncvar_get(nc)  # get values
        lon <- ncdf4::ncvar_get(nc, "lon")
        lat <- ncdf4::ncvar_get(nc, "lat")
        time <- ncdf4::ncvar_get(nc, "time")
        # 
        nCol <- base::length(lon[, 1])
        nRow <- base::length(lon[1, ])
        nTime <- base::length(time)
        # compute dates
        origin <- ncdf4::ncatt_get(nc, "time", attname = "units")[2]
        origin <- base::strsplit(origin$value, " ")
        origin <- base::unlist(origin)[3]
        start.date <- base::as.Date(origin) + time[1] + 30  # the result is consitent with cdo sinfo
        dates <- base::seq(as.Date(start.date), by = timeInt, length = nTime)
        start.date <- min(dates)
        end.date <- max(dates)
        start.date <- format(as.Date(start.date), "%Y-%m")
        end.date <- format(as.Date(end.date), "%Y-%m")
        
        lon <- base::matrix(lon, ncol = 1)
        lat <- base::matrix(lat, ncol = 1)
        
        lonLat <- base::data.frame(lon, lat)
        sp::coordinates(lonLat) <- ~lon + lat
        raster::projection(lonLat) <- regular
        lonLat <- sp::spTransform(lonLat, sp::CRS(rotated))
        myExtent <- raster::extent(lonLat)
        # create an empty raster
        r <- raster::raster(ncols = nCol, nrows = nRow)  # empty raster
        # create a raster by looping through all time steps
        cl <- parallel::makePSOCKcluster(numCores)
        doParallel::registerDoParallel(cl)
        myRaster <- foreach::foreach(i = 1:nTime) %dopar% {
            myValues <- data[, , i]
            myValues <- base::apply(base::t(myValues), 2, rev)  # rotate values
            r <- raster::setValues(r, myValues)  # assign values
            raster::extent(r) <- myExtent  # extent using the rotated projection
            r <- r * 1  # this is necessary for the base::do.call function below
        }
        
        myStack <- base::do.call(raster::stack, myRaster)
        parallel::stopCluster(cl)
        myStack <- raster::setZ(myStack, dates, name = "time")
        names(myStack) <- dates
        assign(paste(c("mod", "ref")[id], sep = ""), myStack)
        assign(paste("start.date", c("mod", "ref")[id], sep = "."), start.date)
        assign(paste("end.date", c("mod", "ref")[id], sep = "."), end.date)
        
    }
    
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
    
    # Make a string that summarizes metadata. This will be added to each netcdf file (longname).
    
    # The string can then be accessed like this: names(raster(file.nc))
    
    meta.data.mod <- paste(variable.name, mod.id, "from", start.date, "to", end.date, sep = "_")
    meta.data.ref <- paste(variable.name, ref.id, "from", start.date, "to", end.date, sep = "_")
    meta.data.com <- paste(variable.name, mod.id, "vs", ref.id, "from", start.date, "to", end.date, sep = "_")
    # all extreme outliers are set to NA in the grid they can be marked as a dot in the plot model data
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
    
    
    # Ensure that both data sets use the same grid cells
    mask <- mod.mean * ref.mean
    mask <- mask - mask + 1
    
    mod <- mod * mask
    ref <- ref * mask
    
    
    # Compute zonal mean
    
    data.list <- list(mod, ref)
    names <- c("mod", "ref")
    for (i in 1:length(data.list)) {
        data <- raster::stack(data.list[i])
        data.mean <- raster::mean(data, na.rm = TRUE)  # time mean
        values <- raster::getValues(data.mean)
        xy <- sp::coordinates(data.mean)
        data <- data.frame(xy, values)
        
        sp::coordinates(data) <- ~x + y
        raster::projection(data) <- my.projection
        data <- sp::spTransform(data, sp::CRS(regular))
        lonLat <- sp::coordinates(data)
        data <- data.frame(lonLat, data$values)
        colnames(data) <- c("lon", "lat", "values")
        data <- stats::na.omit(data)  # omit rows with NA
        
        zone <- round(data$lat, 0)
        data <- data.frame(zone, data)
        index <- list(data$zone)
        
        # define function that computes quantiles
        fun.q25 <- function(x) {
            q25 <- stats::quantile(x, probs = c(0.25), na.rm = TRUE)
            return(q25)
        }
        
        fun.q75 <- function(x) {
            q75 <- stats::quantile(x, probs = c(0.75), na.rm = TRUE)
            return(q75)
        }
        
        zonal.mean.mean <- tapply(data$values, index, mean)
        zonal.mean.min <- tapply(data$values, index, fun.q25)
        zonal.mean.max <- tapply(data$values, index, fun.q75)
        
        zone <- as.numeric(rownames(zonal.mean.mean))
        
        # combine values rbind(zonal.mean.min, zonal.mean.max)
        zonal.mean <- data.frame(zone, zonal.mean.min, zonal.mean.max, zonal.mean.mean)
        zonal.mean <- stats::na.omit(zonal.mean)  # omit rows with NA
        colnames(zonal.mean) <- c("lat", "q25", "q75", "mean")
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
    my.ylim <- c(min(zonal.mean.mod[2], zonal.mean.ref[2]), max(zonal.mean.mod[3], zonal.mean.ref[3]))
    
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
    
    # reference data polygons for uncertainty range
    zonal.mean <- zonal.mean.ref
    zone <- zonal.mean$lat
    zonal.mean.min <- zonal.mean$q25
    zonal.mean.max <- zonal.mean$q75
    zonal.mean.mean <- zonal.mean$mean
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
    graphics::plot(zonal.mean.ref$lat, zonal.mean.ref$q75, main = paste(long.name, my.title, sep = "\n"), type = "l", xlab = "degrees latitude", 
        ylab = legend.bar.text, xlim = my.xlim, ylim = my.ylim, col = NA, las = 1)
    # ref
    graphics::polygon(poly.x.ref, poly.y.ref, col = my.col.ref.range, border = NA)
    graphics::lines(zone.ref, zonal.mean.ref$mean, col = my.col.ref, lwd = 2)
    # mod
    graphics::polygon(poly.x.mod, poly.y.mod, col = my.col.mod.range, border = NA)
    graphics::lines(zone.mod, zonal.mean.mod$mean, col = my.col.mod, lwd = 2)
    # legend
    graphics::legend("topright", c("model mean and IQR", "reference mean and IQR"), col = c(my.col.mod, my.col.ref), lwd = 2, 
        bty = "n")
    # ticks
    graphics::axis(1, at = seq(-90, 90, 10), labels = FALSE, tcl = 0.3)
    graphics::axis(3, at = seq(-90, 90, 10), labels = FALSE, tcl = 0.3)
    graphics::axis(4, labels = FALSE, tcl = 0.3)
    
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}
if (getRversion() >= "2.15.1") utils::globalVariables(c("zonal.mean.mod", "zonal.mean.ref", "start.date.mod", "start.date.ref", 
    "end.date.mod", "end.date.ref"))


