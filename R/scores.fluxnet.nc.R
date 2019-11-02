################################################################################ 
#' Scores for FLUXNET reference data in NetCDF format
#' @description This function compares model output against
#' FLUXNET measurements in NetCDF format. The performance of a model is
#' expressed through scores that range from zero to one, where increasing values
#' imply better performance. These scores are computed in five steps:
#' \eqn{(i)} computation of a statistical metric,
#' \eqn{(ii)} nondimensionalization,
#' \eqn{(iii)} conversion to unit interval,
#' \eqn{(iv)} spatial integration, and
#' \eqn{(v)} averaging scores computed from different statistical metrics.
#' The latter includes the bias, root-mean-square error, phase shift,
#' inter-annual variability, and spatial distribution. The corresponding equations
#' are documented in \code{\link{amber-package}}.
#'
#' @param long.name A string that gives the full name of the variable, e.g. 'Gross primary productivity'
#' @param nc.mod A string that gives the path and name of the netcdf file that contains the model output, e.g. '/home/model_gpp.nc'
#' @param nc.ref A string that gives the path and name of the netcdf file that contains the reference data output, e.g. '/home/reference_gpp.nc'
#' @param mod.id A string that identifies the source of the reference data set, e.g. 'CanESM2'
#' @param ref.id A string that identifies the source of the reference data set, e.g. 'MODIS'
#' @param unit.conv.mod A number that is used as a factor to convert the unit of the model data, e.g. 86400
#' @param unit.conv.ref A number that is used as a factor to convert the unit of the reference data, e.g. 86400
#' @param variable.unit A string that gives the final units using LaTeX notation, e.g. 'gC m$^{-2}$ day$^{-1}$'
#' @param score.weights R object that gives the weights of each score (\eqn{S_{bias}}, \eqn{S_{rmse}}, \eqn{S_{phase}}, \eqn{S_{iav}}, \eqn{S_{dist}})
#' that are used for computing the overall score, e.g. c(1,2,1,1,1)
#' @param rotate.me logical: TRUE if you want longitudes to range from -180 to
#' 180 degrees and FALSE if you want longitudes to range from 0 to 360 degrees
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
#' @param numCores An integer that defines the number of cores, e.g. 2
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @return (1) Figures in PDF format that show global maps of
#' the model data at the location of FLUXNET sites
#' (mean, \eqn{mod.mean}; interannual-variability, \eqn{mod.iav}; month of annual cycle maximum, \eqn{mod.max.month}),
#' the reference data
#' (mean, \eqn{ref.mean}; interannual-variability, \eqn{ref.iav}; month of annual cycle maximum, \eqn{ref.max.month}),
#' statistical metrics
#' (bias, \eqn{bias}; centralized root mean square error, \eqn{crmse}; time difference of the annual cycle maximum, \eqn{phase}),
#' and scores
#' (bias score, \eqn{bias.score}; root mean square error score, \eqn{rmse.score}; inter-annual variability score \eqn{iav.score}; annual cycle score (\eqn{phase.score}).
#'
#' (2) Four text files: (i) score values and (ii) score inputs for each individual
#' site, and (iii) score values and (iv) score inputs averaged across sites.
#' when averaging over all station.
#' @examples
#'
#' library(amber)
#' library(doParallel)
#' library(foreach)
#' library(ncdf4)
#' library(parallel)
#' library(raster)

#' long.name <- 'Gross primary productivity'
#' nc.mod <- system.file('extdata/modelRegular', 'gpp_monthly.nc', package = 'amber')
#' nc.ref <- system.file('extdata/referenceRegular', 'gpp_FLUXNET.nc', package = 'amber')
#' mod.id <- 'CLASSIC' # define a model experiment ID
#' ref.id <- 'FLUXNET' # give reference dataset a name
#' unit.conv.mod <- 86400*1000 # optional unit conversion for model data
#' unit.conv.ref <- 86400*1000 # optional unit conversion for reference data
#' variable.unit <- 'gC m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
#' score.weights <- c(1,2,1,1,1) # score weights of S_bias, S_rmse, S_phase, S_iav, S_dist
#' # global plot
#' scores.fluxnet.nc(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#'  unit.conv.ref, variable.unit, score.weights)
#'
#' \donttest{
#' # regional plot
#' scores.fluxnet.nc(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit, score.weights, my.xlim = c(-150, -60), my.ylim = c(20, 60),
#' plot.width = 6, plot.height = 3.8)
#' scores.fluxnet.nc(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit, score.weights)
#'
#' # (2) Example for data on a rotated grid
#' nc.mod <- system.file('extdata/modelRotated', 'gpp_monthly.nc', package = 'amber')
#' nc.ref <- system.file('extdata/referenceRegular', 'gpp_FLUXNET.nc', package = 'amber')
#' mod.id <- 'CLASSIC' # define a model experiment ID
#' ref.id <- 'FLUXNET' # give reference dataset a name
#' unit.conv.mod <- 86400*1000 # optional unit conversion for model data
#' unit.conv.ref <- 86400*1000 # optional unit conversion for reference data
#' variable.unit <- 'gC m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
#' score.weights <- c(1,2,1,1,1) # score weights of S_bias, S_rmse, S_phase, S_iav, S_dist
#' rotate.me <- FALSE
#' irregular <- TRUE
#' my.projection <-'+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.'
#' shp.filename <- system.file('extdata/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp',
#'  package = 'amber')
#' my.xlim <- c(-171, 0) # longitude range that you wish to plot
#' my.ylim <- c(32, 78) # latitude range that you wish to plot
#' plot.width <- 7 # plot width
#' plot.height <- 3.8 # plot height
#' numCores = 2
#'
#'scores.fluxnet.nc(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit, score.weights, rotate.me, irregular,
#' my.projection,
#' shp.filename, my.xlim, my.ylim)
#' }
#' @export
scores.fluxnet.nc <- function(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod, unit.conv.ref, variable.unit, score.weights, 
    rotate.me = TRUE, irregular = FALSE, my.projection = "+proj=longlat +ellps=WGS84", shp.filename = system.file("extdata/ne_110m_land/ne_110m_land.shp", 
        package = "amber"), my.xlim = c(-180, 180), my.ylim = c(-60, 85), plot.width = 8, plot.height = 3.8, numCores = 2, 
    outputDir = FALSE) {
    
    # (I) Data preparation -----------------------------------------------------
    
    # (1) Reproject data from an irregular to a regular grid if 'irregular=TRUE'
    
    if (irregular == TRUE) {
        regular <- "+proj=longlat +ellps=WGS84"
        rotated <- my.projection
        # get variable name
        nc <- ncdf4::nc_open(nc.mod)
        variable.name <- base::names(nc[["var"]])
        variable.name <- variable.name[1]
        variable.name <- ifelse(variable.name == "burntFractionAll", "burnt", variable.name)  # rename burntFractionAll to shorter name
        variable.name <- toupper(variable.name)  # make variable name upper-case
        ncdf4::nc_close(nc)
        
        # process model data
        my.nc <- nc.mod
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
        dates <- base::seq(as.Date(start.date), by = "month", length = nTime)
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
        mod <- myStack
        start.date.mod <- start.date
        end.date.mod <- end.date
        
        # process reference data
        nc <- ncdf4::nc_open(nc.ref)
        variable.name.ref <- names(nc[["var"]])
        variable.name.ref <- variable.name.ref[length(variable.name.ref)]  # take the last variable
        # get date
        time <- nc$dim[[1]]  # days since 1850-01-01 00:00:00
        n <- length(time$vals)  # number of months in file
        origin <- strsplit(time$units, " ")
        origin <- unlist(origin)[3]
        start.date.ref <- as.Date(origin) + time$vals[1] + 30  # the result is consitent with cdo sinfo
        dates.ref <- seq(as.Date(start.date.ref), by = "month", length = n)
        dates.ref <- format(as.Date(dates.ref), "%Y-%m")  # only year and month
        start.date.ref <- min(dates.ref)
        end.date.ref <- max(dates.ref)
        # convert data to data frame
        lon <- ncdf4::ncvar_get(nc, varid = "lon")
        lat <- ncdf4::ncvar_get(nc, varid = "lat")
        
        lon <- base::matrix(lon, ncol = 1)
        lat <- base::matrix(lat, ncol = 1)
        
        # reproject coordinates to rotated grid
        lonLat <- base::data.frame(lon, lat)
        sp::coordinates(lonLat) <- ~lon + lat
        raster::projection(lonLat) <- regular
        lonLat <- sp::spTransform(lonLat, sp::CRS(rotated))
        lon <- sp::coordinates(lonLat)[, 1]
        lat <- sp::coordinates(lonLat)[, 2]
        data <- ncdf4::ncvar_get(nc, varid = variable.name.ref)
        ncdf4::nc_close(nc)
        data <- data.frame(lon, lat, data * unit.conv.ref)
        colnames(data) <- c("lon", "lat", dates.ref)
        site.id <- paste("site", seq(1, nrow(data), 1), sep = "_")
        rownames(data) <- site.id
        ref <- data
        
    } else {
        
        # (2) Process data if 'irregular=FALSE' --------------------------------
        
        # model data
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
        myExtent <- c(-180, 180, -90, 90)  # used later
        # reference data
        nc <- ncdf4::nc_open(nc.ref)
        variable.name.ref <- names(nc[["var"]])
        variable.name.ref <- variable.name.ref[length(variable.name.ref)]  # take the last variable
        # get date
        time <- nc$dim[[1]]  # days since 1850-01-01 00:00:00
        n <- length(time$vals)  # number of months in file
        origin <- strsplit(time$units, " ")
        origin <- unlist(origin)[3]
        start.date.ref <- as.Date(origin) + time$vals[1] + 30  # the result is consitent with cdo sinfo
        dates.ref <- seq(as.Date(start.date.ref), by = "month", length = n)
        dates.ref <- format(as.Date(dates.ref), "%Y-%m")  # only year and month
        start.date.ref <- min(dates.ref)
        end.date.ref <- max(dates.ref)
        # convert data to data frame
        lon <- ncdf4::ncvar_get(nc, varid = "lon")
        lat <- ncdf4::ncvar_get(nc, varid = "lat")
        data <- ncdf4::ncvar_get(nc, varid = variable.name.ref)
        ncdf4::nc_close(nc)
        data <- data.frame(lon, lat, data * unit.conv.ref)
        colnames(data) <- c("lon", "lat", dates.ref)
        site.id <- paste("site", seq(1, nrow(data), 1), sep = "_")
        rownames(data) <- site.id
        ref <- data
        
    }
    
    #---------------------------------------------------------------------------
    
    # find common time period
    start.date <- max(start.date.mod, start.date.ref)
    end.date <- min(end.date.mod, end.date.ref)
    # subset common time period for model data
    mod <- mod[[which(format(as.Date(raster::getZ(mod)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(mod)), "%Y-%m") <= 
        end.date)]]
    # subset common time period for reference data get a sequence of common dates
    dates <- seq(from = as.Date(paste(start.date, "15", sep = "-")), to = as.Date(paste(end.date, "15", sep = "-")), by = "month")
    dates <- format(as.Date(dates), "%Y-%m")
    common.dates <- data.frame(dates)
    colnames(common.dates) <- "dates"
    ## add an index to each date
    index <- seq(1, length(dates.ref), 1)
    dates.ref <- data.frame(index, dates.ref)
    colnames(dates.ref) <- c("index", "dates")
    ## merge common dates and reference dates the index will be used to subset the reference data set if the time covered by
    ## the reference data exceeds the time covered by the model data
    column.id <- merge(dates.ref, common.dates, by = "dates")
    column.id <- column.id$index
    ref <- data.frame(ref[2 + column.id])  # 2+ because column 1 and 2 are are lon and lat
    colnames(ref) <- c(dates)
    ref <- intFun.site.points(lon, lat, ref)
    # compute statistics for each site by extracting model data for each site
    if (rotate.me == TRUE) {
        mod <- raster::rotate(mod)
    }
    mod <- mod * unit.conv.mod
    mod <- raster::extract(mod, ref, method = "bilinear")
    mod <- data.frame(mod)
    rownames(mod) <- site.id
    ref <- as.data.frame(ref)
    ref <- ref[-c(1, 2)]  # drop lon lat column
    # check whether number of columns and rows match between mod and ref
    stopifnot(ncol(ref) == ncol(mod))
    stopifnot(nrow(ref) == nrow(mod))
    # make a string that summarizes metadata
    meta.data.mod <- paste(variable.name, mod.id, "from", start.date, "to", end.date, sep = "_")
    meta.data.ref <- paste(variable.name, ref.id, "from", start.date, "to", end.date, sep = "_")
    meta.data.com <- paste(variable.name, mod.id, "vs", ref.id, "from", start.date, "to", end.date, sep = "_")
    # statistical analysis the time period of the reference data varies among sites some sites are located in the model
    # ocean exclude data in one data set that is not available in the other data set
    mask.mod <- (mod - mod + 1)
    mask.ref <- (ref - ref + 1)
    mask <- mask.mod * mask.ref
    mod <- mod * mask
    ref <- ref * mask
    # get number of years of data coverage of each station
    n.months <- apply(mask, 1, sum, na.rm = TRUE)
    n.years <- round(n.months/12)
    n.years <- intFun.site.points(lon, lat, n.years)
    # transpose data, where rows are dates and columns are sites this is necessary for using mapply below
    mod <- data.frame(t(mod))
    ref <- data.frame(t(ref))
    # (1) bias
    mod.mean <- apply(mod, 2, mean, na.rm = TRUE)  # time mean
    ref.mean <- apply(ref, 2, mean, na.rm = TRUE)  # time mean
    weights <- ref.mean  # weights used for spatial integral
    bias <- mod.mean - ref.mean  # time mean
    ref.sd <- apply(ref, 2, sd, na.rm = TRUE)  # standard deviation of reference data
    epsilon_bias <- abs(bias)/ref.sd  # relative error
    bias.score <- exp(-epsilon_bias)  # bias score as a function of space
    S_bias_not.weighted <- mean(bias.score, na.rm = TRUE)  # scalar score (not weighted)
    # calculate the weighted scalar score
    a <- bias.score * weights
    b <- sum(a, na.rm = TRUE)  # this is a scalar, the sum of all values
    S_bias_weighted <- b/sum(weights, na.rm = TRUE)  # scalar score (weighted)
    # compute global mean values of score input(s)
    mod.mean.scalar <- mean(mod.mean, na.rm = TRUE)  # global mean value
    ref.mean.scalar <- mean(ref.mean, na.rm = TRUE)  # global mean value
    bias.scalar <- mean(bias, na.rm = TRUE)  # global mean value
    ref.sd.scalar <- mean(ref.sd, na.rm = TRUE)
    epsilon_bias.scalar <- mean(epsilon_bias, na.rm = TRUE)
    # (2) root mean square error (rmse)
    rmse <- mapply(intFun.rmse, mod, ref)  # compute rmse
    mod.anom <- data.frame(apply(mod, 2, intFun.anom))  # compute anomalies
    ref.anom <- data.frame(apply(ref, 2, intFun.anom))  # compute anomalies
    crmse <- mapply(intFun.crmse, mod.anom, ref.anom)
    epsilon_rmse <- crmse/ref.sd  # relative error
    rmse.score <- exp(-epsilon_rmse)  # rmse score as a function of space
    S_rmse_not.weighted <- mean(rmse.score, na.rm = TRUE)  # scalar score (not weighted)
    # calculate the weighted scalar score
    a <- rmse.score * weights  # this is a raster
    b <- sum(a, na.rm = TRUE)  # this is a scalar, the sum up all values
    S_rmse_weighted <- b/sum(weights, na.rm = TRUE)  # scalar score (weighted)
    # compute global mean values of score input(s)
    rmse.scalar <- mean(rmse, na.rm = TRUE)  # global mean value
    crmse.scalar <- mean(crmse, na.rm = TRUE)  # global mean value
    epsilon_rmse.scalar <- mean(epsilon_rmse, na.rm = TRUE)
    # (3) phase shift make monthly means month index is the same for mod and ref, given the temporal subset made above
    month <- format(as.Date(paste(dates, "15", sep = "-")), format = "%m")
    month <- as.numeric(month)
    mod <- data.frame(month, mod)
    ref <- data.frame(month, ref)
    # compute climatological mean monthly values
    index <- list(mod$month)
    mod.clim.mly <- apply(mod, 2, function(x) {
        tapply(x, index, mean, na.rm = TRUE)
    })
    ref.clim.mly <- apply(ref, 2, function(x) {
        tapply(x, index, mean, na.rm = TRUE)
    })
    mod.clim.mly.mm <- mod.clim.mly  # will be used when computing iav further below
    ref.clim.mly.mm <- ref.clim.mly  # will be used when computing iav further below
    # drop 'month' column
    mod.clim.mly <- subset(mod.clim.mly, select = -c(month))
    ref.clim.mly <- subset(ref.clim.mly, select = -c(month))
    # find month with max value
    mod.max.month <- apply(mod.clim.mly, 2, which.max)
    ref.max.month <- apply(ref.clim.mly, 2, which.max)
    mod.max.month <- as.numeric(mod.max.month)
    mod.max.month <- data.frame(as.numeric(mod.max.month))
    ref.max.month <- data.frame(as.numeric(ref.max.month))
    # get shortest time distance between these months
    abs.diff <- abs(mod.max.month - ref.max.month)  # absolute difference from 0 to 12 months
    phase <- apply(abs.diff, 2, intFun.theta)  # shortest distance from 0 to 6 months (theta)
    phase.score <- 0.5 * (1 + cos(2 * pi * phase/12))  # score from 0 (6 months) to 1 (0 months)
    S_phase_not.weighted <- mean(phase.score, na.rm = TRUE)  # scalar score (not weighted)
    # calculate the weighted scalar score
    a <- phase.score * weights  # this is spatial data
    b <- sum(a, na.rm = TRUE)  # this is a scalar, the sum up all values
    S_phase_weighted <- b/sum(weights, na.rm = TRUE)  # scalar score (weighted)
    # compute global mean values of score input(s)
    phase.scalar <- mean(phase, na.rm = TRUE)  # global mean value
    mod.max.month.scalar <- mean(mod.max.month[, 1], na.rm = TRUE)  # global mean value
    ref.max.month.scalar <- mean(ref.max.month[, 1], na.rm = TRUE)  # global mean value
    # (4) interannual variability
    
    # This approach assumes that all data start in Jan all months after the last Dec will be dropped if data does not end in
    # Dec
    mod.anom <- intFun.anom.mly(mod, mod.clim.mly.mm)
    ref.anom <- intFun.anom.mly(ref, ref.clim.mly.mm)
    mod.iav <- apply(mod.anom, 2, intFun.iav)
    ref.iav <- apply(ref.anom, 2, intFun.iav)
    
    # set values close to zero to NA
    ref.iav.na <- ref.iav
    ref.iav.na[ref.iav.na < 10^(-5)] <- NA
    
    epsilon_iav <- abs((mod.iav - ref.iav))/ref.iav.na  # I changed Eq. 26 so that epsilon_iav > =  0
    iav.score <- exp(-epsilon_iav)  # iav score as a function of space
    S_iav_not.weighted <- mean(iav.score, na.rm = TRUE)  # scalar score (not weighted)
    # calculate the weighted scalar score
    a <- iav.score * weights  # this is a raster
    b <- sum(a, na.rm = TRUE)  # this is a scalar, the sum up all values
    S_iav_weighted <- b/sum(weights, na.rm = TRUE)  # scalar score (weighted)
    # compute global mean values of score input(s)
    mod.iav.scalar <- mean(mod.iav, na.rm = TRUE)  # global mean value
    ref.iav.scalar <- mean(ref.iav, na.rm = TRUE)  # global mean value
    epsilon_iav.scalar <- mean(epsilon_iav, na.rm = TRUE)  # global mean value
    
    # (5) dist
    mod.sigma.scalar <- sd(mod.mean, na.rm = TRUE)  # standard deviation of period mean data
    ref.sigma.scalar <- sd(ref.mean, na.rm = TRUE)  # standard deviation of period mean data
    sigma <- mod.sigma.scalar/ref.sigma.scalar
    y <- mod.mean
    x <- ref.mean
    reg <- stats::lm(y ~ x)
    R <- sqrt(summary(reg)$r.squared)
    S_dist <- 2 * (1 + R)/(sigma + 1/sigma)^2  # weighting does not apply
    # scores
    w.bias <- score.weights[1]
    w.rmse <- score.weights[2]
    w.phase <- score.weights[3]
    w.iav <- score.weights[4]
    w.dist <- score.weights[5]
    # not weighted
    S_bias <- S_bias_not.weighted
    S_rmse <- S_rmse_not.weighted
    S_phase <- S_phase_not.weighted
    S_iav <- S_iav_not.weighted
    # weight importance of statisitcal metrics and compute overall score
    S_overall <- (w.bias * S_bias + w.rmse * S_rmse + w.phase * S_phase + w.iav * S_iav + w.dist * S_dist)/(w.bias + w.rmse + 
        w.phase + w.iav + w.dist)
    scores <- data.frame(variable.name, ref.id, S_bias, S_rmse, S_phase, S_iav, S_dist, S_overall)
    scores_not.weighted <- scores
    # weighted (except for S_dist)
    S_bias <- S_bias_weighted
    S_rmse <- S_rmse_weighted
    S_phase <- S_phase_weighted
    S_iav <- S_iav_weighted
    S_overall <- (w.bias * S_bias + w.rmse * S_rmse + w.phase * S_phase + w.iav * S_iav + w.dist * S_dist)/(w.bias + w.rmse + 
        w.phase + w.iav + w.dist)
    scores <- data.frame(variable.name, ref.id, S_bias, S_rmse, S_phase, S_iav, S_dist, S_overall)
    scores_weighted <- scores
    # 
    scores <- rbind(scores_not.weighted, scores_weighted)
    rownames(scores) <- c("not.weighted", "weighted")
    if (outputDir != FALSE) {
        utils::write.table(scores, paste(outputDir, "/", "scorevalues", "_", meta.data.com, sep = ""))
    }
    # get all score values in case you want to compare this run against another run using a significance test
    dist.score <- rep(S_dist, length(bias.score))
    all.score.values <- data.frame(bias.score, rmse.score, phase.score, iav.score, dist.score)
    colnames(all.score.values) <- c("bias.score", "rmse.score", "phase.score", "iav.score", "dist.score")
    all.score.values[is.na(all.score.values)] <- NA  # converts all NaN to NA
    if (outputDir != FALSE) {
        utils::write.table(all.score.values, paste(outputDir, "/", "allscorevalues", "-", variable.name, "-", ref.id, sep = ""))
    }
    # selected score inputs
    scoreinputs <- data.frame(long.name, variable.name, ref.id, variable.unit, mod.mean.scalar, ref.mean.scalar, bias.scalar, 
        ref.sd.scalar, epsilon_bias.scalar, S_bias_not.weighted, rmse.scalar, crmse.scalar, ref.sd.scalar, epsilon_rmse.scalar, 
        S_rmse_not.weighted, mod.max.month.scalar, ref.max.month.scalar, phase.scalar, S_phase_not.weighted, mod.iav.scalar, 
        ref.iav.scalar, epsilon_iav.scalar, S_iav_not.weighted, mod.sigma.scalar, ref.sigma.scalar, sigma, R, S_dist)
    if (outputDir != FALSE) {
        utils::write.table(scoreinputs, paste(outputDir, "/", "scoreinputs", "_", meta.data.com, sep = ""))
    }
    # min.max.int
    mmi.bias <- intFun.min.max.int.bias(bias)
    mmi.bias.score <- c(0, 1, 0.1)
    mmi.crmse <- intFun.min.max.int(crmse)
    mmi.rmse.score <- c(0, 1, 0.1)
    mmi.phase <- c(0, 6, 1)
    mmi.phase.score <- c(0, 1, 0.1)
    mmi.iav.score <- c(0, 1, 0.1)
    # min.max.int.mod.ref
    mmi.mean <- intFun.min.max.int.mod.ref(mod.mean, ref.mean)
    mmi.max.month <- c(1, 12, 1)
    mmi.iav <- intFun.min.max.int.mod.ref(mod.iav, ref.iav)
    # add metadata:
    
    # 1. figure title (e.g. Mean_nee_ModID_123_from_1982-01_to_2008-12)
    
    # 2. min, max, interval used in legend (e.g. 0, 1, 0.1)
    
    # 3. legend bar text (e.g. 'score (-)')
    attr(mod.mean, "metadata") <- list(paste("Mean", meta.data.mod, sep = "_"), mmi.mean, variable.unit)
    attr(ref.mean, "metadata") <- list(paste("Mean", meta.data.ref, sep = "_"), mmi.mean, variable.unit)
    attr(bias, "metadata") <- list(paste("Bias", meta.data.com, sep = "_"), mmi.bias, variable.unit)
    attr(bias.score, "metadata") <- list(paste("Bias_score", meta.data.com, sep = "_"), mmi.bias.score, "score (-)")
    attr(crmse, "metadata") <- list(paste("CRMSE", meta.data.com, sep = "_"), mmi.crmse, variable.unit)
    attr(rmse.score, "metadata") <- list(paste("RMSE_score", meta.data.com, sep = "_"), mmi.rmse.score, "score (-)")
    attr(mod.max.month, "metadata") <- list(paste("Month_with_max", meta.data.mod, sep = "_"), mmi.max.month, "month")
    attr(ref.max.month, "metadata") <- list(paste("Month_with_max", meta.data.ref, sep = "_"), mmi.max.month, "month")
    attr(phase, "metadata") <- list(paste("Diff_in_max_month", meta.data.com, sep = "_"), mmi.phase, "month")
    attr(phase.score, "metadata") <- list(paste("Seasonality_score", meta.data.com, sep = "_"), mmi.phase.score, "score (-)")
    attr(mod.iav, "metadata") <- list(paste("IAV", meta.data.mod, sep = "_"), mmi.iav, variable.unit)
    attr(ref.iav, "metadata") <- list(paste("IAV", meta.data.ref, sep = "_"), mmi.iav, variable.unit)
    attr(iav.score, "metadata") <- list(paste("IAV_score", meta.data.com, sep = "_"), mmi.iav.score, "score (-)")
    # write data to file
    stat.metric <- data.frame(lon, lat, mod.mean, ref.mean, bias, bias.score, crmse, rmse.score, mod.max.month, ref.max.month, 
        phase, phase.score, mod.iav, ref.iav, iav.score)
    colnames(stat.metric) <- c("lon", "lat", "mod.mean", "ref.mean", "bias", "bias.score", "crmse", "rmse.score", "mod.max.month", 
        "ref.max.month", "phase", "phase.score", "mod.iav", "ref.iav", "iav.score")
    my.filename <- paste(variable.name, "FLUXNET", sep = "_")
    if (outputDir != FALSE) {
        utils::write.table(stat.metric, paste(outputDir, my.filename, sep = "/"))
    }
    
    land <- intFun.coast(my.xlim, my.ylim, my.projection, shp.filename)  # reproject coastline
    # loop through all layers of the raster stack
    lon2 <- stat.metric$lon
    lat2 <- stat.metric$lat
    for (i in 3:ncol(stat.metric)) {
        data <- stat.metric[i]
        cname <- base::names(data)
        data <- data.frame(lon2, lat2, data)
        data <- stats::na.omit(data)
        lon <- data$lon2
        lat <- data$lat2
        data <- data.frame(data[, cname])
        colnames(data) <- cname
        values <- data
        colnames(values) <- "values"
        data <- intFun.site.points(lon, lat, values)  # make spatial points
        x <- base::get(cname)  # get the data
        my.attributes <- attributes(x)
        meta <- my.attributes$metadata
        id <- unlist(meta[1])
        my.title <- gsub("_", " ", id)
        min.max.int <- unlist(meta[2])
        legend.bar.text <- latex2exp::TeX(meta[[3]])
        # for legend
        min <- min.max.int[1]
        max <- min.max.int[2]
        interval <- min.max.int[3]
        my.breaks <- round(seq(min, max, interval), 3)  # breaks of the colors
        my.labels <- round(seq(min, max, interval), 3)  # locations where to set the labels
        my.col <- viridis::viridis(n = length(my.breaks) - 1, direction = -1)
        my.col.bias <- scico::scico(n = length(my.breaks) - 1, palette = "vik")
        my.col.phase <- grDevices::rainbow(n = length(my.breaks) - 1)
        if (i == 5) 
            {
                my.col <- my.col.bias
            }  # divergent color scheme for bias plots
        if (i == 9) 
            {
                my.col <- my.col.phase
            }  # circular color scheme for phase
        if (i == 10) 
            {
                my.col <- my.col.phase
            }  # circular color scheme for phase
        if (i == 11) 
            {
                my.col <- my.col.phase
            }  # circular color scheme for phase
        
        my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
        my.legend.args <- list(text = legend.bar.text, side = 2, font = 1, line = 1, cex = 1)
        
        colors <- cut(values$values, breaks = my.breaks, labels = my.col, include.lowest = TRUE)
        colors <- toString(colors)
        colors <- unlist(strsplit(colors, split = ", "))
        
        # plot
        oldpar <- graphics::par(mfrow = c(1, 2))
        on.exit(graphics::par(oldpar))
        my.plotname <- paste(my.filename, cname, sep = "_")
        my.plotname <- gsub("_", "-", my.plotname)
        my.plotname <- gsub(".", "-", my.plotname, fixed = TRUE)
        if (outputDir != FALSE) {
            grDevices::pdf(paste(outputDir, "/", my.plotname, ".pdf", sep = ""), width = plot.width, height = plot.height)
        }
        graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(3, 3, 3, 4), lwd = 1, cex = 1)
        # create a dummy raster layer
        if (irregular == TRUE) {
            my.xlim <- c(raster::extent(land)[1], raster::extent(land)[2])
            my.ylim <- c(raster::extent(land)[3], raster::extent(land)[4])
        }
        # mod.mean
        dummy <- stats::runif(360 * 180, min = min, max = max)
        dummy <- matrix(dummy, nrow = 180)
        dummy <- raster::raster(dummy)
        raster::extent(dummy) <- myExtent
        raster::plot(dummy, col = NA, xlim = my.xlim, ylim = my.ylim, main = paste(long.name, my.title, sep = "\n"), axes = FALSE, 
            legend = FALSE)
        raster::plot(land, col = "grey", border = NA, add = TRUE)
        graphics::points(data, pch = 16, col = colors)
        if (irregular == FALSE) {
            graphics::axis(1, labels = TRUE, tcl = 0.3)
            graphics::axis(2, labels = TRUE, tcl = 0.3, las = 2)
        }
        plot(dummy, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args, 
            legend.width = 1.5, legend.shrink = 1, font = 1)
        if (outputDir != FALSE) {
            grDevices::dev.off()
        }
    }
    
}

