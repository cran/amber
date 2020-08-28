################################################################################
#' Scores for gridded reference data with a varying time dimension
#' @description This function compares model output against remote-sensing
#' based reference data that vary in time. The performance of a model is
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
#' @param mod.id A string that identifies the source of the reference data set, e.g. 'CLASSIC'
#' @param ref.id A string that identifies the source of the reference data set, e.g. 'MODIS'
#' @param unit.conv.mod A number that is used as a factor to convert the unit of the model data, e.g. 86400
#' @param unit.conv.ref A number that is used as a factor to convert the unit of the reference data, e.g. 86400
#' @param variable.unit A string that gives the final units using LaTeX notation, e.g. 'gC m$^{-2}$ day$^{-1}$'
#' @param score.weights R object that gives the weights of each score (\eqn{S_{bias}}, \eqn{S_{rmse}}, \eqn{S_{phase}}, \eqn{S_{iav}}, \eqn{S_{dist}})
#' that are used for computing the overall score, e.g. c(1,2,1,1,1)
#' @param outlier.factor A number that is used to define outliers, e.g. 10.
#'  Plotting raster objects that contain extreme outliers lead to figures where
#'  most grid cells are presented by a single color since the color legend covers
#'  the entire range of values. To avoid this, the user may define outliers that
#'  will be masked out and marked with a red dot. Outliers are all values that
#'  exceed the interquartile range multiplied by the outlier factor defined here.
#' @param irregular Logical. If TRUE the data is converted from an irregular to a regular grid. Default is FALSE.
#' @param my.projection A string that defines the projection of the irregular grid
#' @param numCores An integer that defines the number of cores, e.g. 2
#' @param timeInt A string that gives the time interval of the model data, e.g. 'month' or 'year'
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @param myLevel A number that determines what level of the output netCDF file to use.
#' This is relevant for files with multiple levels, which applies to soil data.
#' By default, myLevel is set to 1.
#' @param variable.name A string with the variable name, e.g. 'GPP'. If FALSE, the variable name stored in the NetCDF file will be used instead. Default is FALSE.
#' @param phaseMinMax A string (either 'phaseMax' or 'phaseMin') that determines
#' whether to assess the seasonal peak as a maximum or a minimum. The latter may be appropriate for variables
#' that tend to be negative, such as net longwave radiation or net ecosystem exchange.
#'
#' @return (1) A list that contains three elements. The first element is a
#' raster stack with model data
#' (mean, \eqn{mod.mean}; standard deviation; interannual-variability, \eqn{mod.iav}; monthly mean climatology; month of annual cycle maximum, \eqn{mod.max.month}),
#' the reference data (mean, \eqn{ref.mean}; standard deviation; interannual-variability, \eqn{ref.iav}; monthly mean climatology; month of annual cycle maximum, \eqn{ref.max.month}),
#' statistical metrics
#' (bias, \eqn{bias}; centralized root mean square error, \eqn{crmse}; time difference of the annual cycle maximum, \eqn{phase}),
#' and scores (bias score, \eqn{bias.score}; root mean square error score, \eqn{rmse.score}; inter-annual variability score \eqn{iav.score}; annual cycle score (\eqn{phase.score}).
#' The second and third element of the list are spatial
#' point data frames that give the model and reference outliers, respectively.
#' Most of the content of the list can be plotted using \link{plotGrid}. The only exception is monthly mean climatology.
#'
#' (2) NetCDF files for each of the statistical variables listed above.
#'
#' (3) Three text files: (i) score values and (ii) score inputs averaged across
#' the entire grid, and (iii) score values for each individual grid cell.
#'
#' @examples
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
#' # (1) Global plots on a regular grid
#' long.name <- 'Gross primary productivity'
#' nc.mod <- system.file('extdata/modelRegular', 'gpp_monthly.nc', package = 'amber')
#' nc.ref <- system.file('extdata/referenceRegular', 'gpp_GBAF_128x64.nc', package = 'amber')
#' mod.id <- 'CLASSIC' # define a model experiment ID
#' ref.id <- 'GBAF' # give reference dataset a name
#' unit.conv.mod <- 86400*1000 # optional unit conversion for model data
#' unit.conv.ref <- 86400*1000 # optional unit conversion for reference data
#' variable.unit <- 'gC m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
#'
#' plot.me <- scores.grid.time(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit)
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
#' my.projection <- '+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.'
#'
#' plot.me <- scores.grid.time(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit, score.weights = c(1, 2, 1, 1, 1), outlier.factor = 1000,
#' irregular = TRUE, my.projection = my.projection, numCores = 2, timeInt = 'month',
#' outputDir = FALSE, myLevel = 1, variable.name = FALSE)
#'
#' # Plot results:
#' irregular <- TRUE # data is on an irregular grid
#' my.projection <- '+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.'
#' # shp.filename <- system.file('extdata/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp',
#' #  package = 'amber')
#' shp.filename <- system.file("extdata/ne_110m_land/ne_110m_land.shp", package = "amber")
#' my.xlim <- c(-171, 0) # longitude range that you wish to plot
#' my.ylim <- c(32, 78) # latitude range that you wish to plot
#' plot.width <- 7 # plot width
#' plot.height <- 3.8 # plot height
#'
#' plotGrid(long.name, plot.me, irregular, my.projection,
#' shp.filename, my.xlim, my.ylim, plot.width, plot.height)
#' } #donttest
#'
#'
#' @export
scores.grid.time <- function(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod, unit.conv.ref, variable.unit, score.weights = c(1,
    2, 1, 1, 1), outlier.factor = 1000, irregular = FALSE, my.projection = "+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.",
    numCores = 2, timeInt = "month", outputDir = FALSE, myLevel = 1, variable.name = FALSE, phaseMinMax = "phaseMax") {

    #---------------------------------------------------------------------------

    # (I) Data preparation

    #---------------------------------------------------------------------------

    # (1) Reproject data from an irregular to a regular grid if 'irregular=TRUE'

    #---------------------------------------------------------------------------
    if (irregular == TRUE) {

        regular <- "+proj=longlat +ellps=WGS84"
        rotated <- my.projection
        # get variable name
        if (variable.name == FALSE) {
            nc <- ncdf4::nc_open(nc.mod)
            variable.name <- base::names(nc[["var"]])
            variable.name <- variable.name[1]
            variable.name <- ifelse(variable.name == "burntFractionAll", "burnt", variable.name)  # rename burntFractionAll to shorter name
            variable.name <- toupper(variable.name)  # make variable name upper-case
            ncdf4::nc_close(nc)
        }
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
            start.date <- base::as.Date(origin) + time[1] + 30  # the result is consistent with cdo sinfo
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
            suppressWarnings(lonLat <- sp::spTransform(lonLat, sp::CRS(rotated)))
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
    } else {

        #-----------------------------------------------------------------------

        # (2) Process data if 'irregular=FALSE'

        #-----------------------------------------------------------------------
        if (variable.name == FALSE) {
            nc <- ncdf4::nc_open(nc.mod)
            variable.name <- names(nc[["var"]])
            ncdf4::nc_close(nc)
            variable.name <- variable.name[length(variable.name)]  # take the last variable (relevant for CanESM5)
            variable.name <- ifelse(variable.name == "burntFractionAll", "burnt", variable.name)  # rename burntFractionAll to shorter name
            variable.name <- toupper(variable.name)  # make variable name upper-case
        }
        mod <- raster::brick(nc.mod, level = myLevel)
        suppressWarnings(mod <- raster::rotate(mod))

        dates.mod <- raster::getZ(mod)
        dates.mod <- format(as.Date(dates.mod), "%Y-%m")  # only year and month
        start.date.mod <- min(dates.mod)
        end.date.mod <- max(dates.mod)

        # reference data
        ref <- raster::brick(nc.ref)
        suppressWarnings(ref <- raster::rotate(ref))
        dates.ref <- raster::getZ(ref)
        dates.ref <- format(as.Date(dates.ref), "%Y-%m")  # only year and month
        start.date.ref <- min(dates.ref)
        end.date.ref <- max(dates.ref)
    }

    #---------------------------------------------------------------------------

    # 1.3 Remaining part applies to both regular and irregular gridded data

    #---------------------------------------------------------------------------

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
    mod.outlier.points <- intFun.grid.outliers.points(mod.mean, outlier.neg, outlier.pos)

    # reference data
    ref.mean <- raster::mean(ref, na.rm = TRUE)  # time mean
    ref.outlier_range <- intFun.grid.define.outlier(ref.mean, outlier.factor)  # define outlier range
    outlier.neg <- ref.outlier_range[1]
    outlier.pos <- ref.outlier_range[2]
    ref.mask_outliers <- intFun.grid.outliers.na(ref.mean, outlier.neg, outlier.pos)
    ref.mask_outliers <- ref.mask_outliers - ref.mask_outliers + 1
    ref <- ref * ref.mask_outliers
    names(ref) <- ref.names
    ref.outlier.points <- intFun.grid.outliers.points(ref.mean, outlier.neg, outlier.pos)

    #---------------------------------------------------------------------------

    # II Statistical analysis

    #---------------------------------------------------------------------------

    # (1) Bias

    #---------------------------------------------------------------------------
    # create a mask to excludes all grid cells that the model and reference data do not have in common.  This mask varies in
    # time.
    mask <- (mod * ref)
    mask <- mask - mask + 1
    mod <- mod * mask
    names(mod) <- mod.names  # this adds the corresponding dates
    ref <- ref * mask
    names(ref) <- ref.names  # this adds the corresponding dates
    # now mod and ref are based on the same grid cells

    mod.mean <- raster::mean(mod, na.rm = TRUE)  # time mean
    ref.mean <- raster::mean(ref, na.rm = TRUE)  # time mean
    weights <- ref.mean  # weights used for spatial integral
    bias <- mod.mean - ref.mean  # time mean
    mod.sd <- raster::calc(mod, fun = sd, na.rm = TRUE)  # standard deviation of model data
    ref.sd <- raster::calc(ref, fun = sd, na.rm = TRUE)  # standard deviation of reference data
    epsilon_bias <- abs(bias)/ref.sd
    epsilon_bias[epsilon_bias == Inf] <- NA  # relative error
    bias.score <- exp(-epsilon_bias)  # bias score as a function of space
    S_bias_not.weighted <- mean(raster::getValues(bias.score), na.rm = TRUE)  # scalar score (not weighted)
    # calculate the weighted scalar score
    a <- bias.score * weights  # this is a raster
    b <- raster::cellStats(a, "sum")  # this is a scalar, the sum up all values
    S_bias_weighted <- b/raster::cellStats(weights, "sum")  # scalar score (weighted)
    # Compute global mean values of score input(s) The mask ensures that mod and ref are based on same grid cells.
    mask <- (mod.mean * ref.mean)
    mask <- mask - mask + 1
    mod.mean.scalar <- mean(raster::getValues(mask * mod.mean), na.rm = TRUE)  # global mean value
    ref.mean.scalar <- mean(raster::getValues(mask * ref.mean), na.rm = TRUE)  # global mean value
    bias.scalar <- mean(raster::getValues(bias), na.rm = TRUE)  # global mean value
    bias.scalar.rel <- (mod.mean.scalar - ref.mean.scalar)/abs(ref.mean.scalar) * 100
    ref.sd.scalar <- mean(raster::getValues(ref.sd), na.rm = TRUE)  # global mean value
    # epsilon_bias.scalar <- stats::median(raster::getValues(epsilon_bias), na.rm = TRUE) # global mean value
    epsilon_bias.scalar <- stats::median(raster::getValues(epsilon_bias), na.rm = TRUE)  # global median value
    # Wilcox significance test
    data <- raster::stack(mod, ref)
    no.data <- raster::mean(data)
    no.data <- no.data - no.data + 1
    data <- raster::calc(data, intFun.grid.na)  # Convert NA to zero, because no NA allowed
    pvalue <- raster::calc(data, intFun.grid.wilcox)
    # set statistically significant differences to NA
    bias.significance <- raster::calc(pvalue, intFun.grid.significance)  #
    # A value of one means that differences are not statistically significant
    bias.significance <- bias.significance - bias.significance + 1
    bias.significance <- bias.significance * no.data  # this excludes grid cells with no data

    #---------------------------------------------------------------------------

    # (2) root mean square error (rmse)

    #---------------------------------------------------------------------------
    rmse <- intFun.rmse(mod, ref)  # rmse
    mod.anom <- mod - mod.mean  # anomaly
    ref.anom <- ref - ref.mean  # anomaly
    crmse <- intFun.crmse(mod.anom, ref.anom)  # centralized rmse
    epsilon_rmse <- crmse/ref.sd
    epsilon_rmse[epsilon_rmse == Inf] <- NA  # relative error
    rmse.score <- exp(-epsilon_rmse)  # rmse score as a function of space
    S_rmse_not.weighted <- mean(raster::getValues(rmse.score), na.rm = TRUE)  # scalar score (not weighted)
    # calculate the weighted scalar score
    a <- rmse.score * weights  # this is a raster
    b <- raster::cellStats(a, "sum")  # this is a scalar, the sum of all values
    S_rmse_weighted <- b/raster::cellStats(weights, "sum")  # scalar score (weighted)
    # compute global mean values of score input(s)
    rmse.scalar <- mean(raster::getValues(rmse), na.rm = TRUE)  # global mean value
    crmse.scalar <- mean(raster::getValues(crmse), na.rm = TRUE)  # global mean value
    # epsilon_rmse.scalar <- stats::median(raster::getValues(epsilon_rmse), na.rm = TRUE) # global mean value
    epsilon_rmse.scalar <- stats::median(raster::getValues(epsilon_rmse), na.rm = TRUE)  # global median value

    #---------------------------------------------------------------------------

    # (3) phase shift

    #---------------------------------------------------------------------------

    index <- format(as.Date(names(ref), format = "X%Y.%m.%d"), format = "%m")
    index <- as.numeric(index)
    mod.clim.mly <- raster::stackApply(mod, index, fun = raster::mean)
    ref.clim.mly <- raster::stackApply(ref, index, fun = raster::mean)

    mod.cycle <- mod.clim.mly  # does not necessarily start in Jan or end in Dec
    ref.cycle <- ref.clim.mly  # does not necessarily start in Jan or end in Dec

    # Ensure that the order is correct
    JanToDec <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6", "index_7", "index_8", "index_9", "index_10",
        "index_11", "index_12")
    mod.clim.mly <- mod.clim.mly[[JanToDec]]
    ref.clim.mly <- ref.clim.mly[[JanToDec]]

    bias.clim.mly <- mod.clim.mly - ref.clim.mly

    # find month of seasonal peak

    # In most cases, we are interested in the timing of the seasonal maximum value

    # In some cases, however, the seasonal peak is a minimum, e.g. NEE = RECO - GPP

    #
    if (phaseMinMax == "phaseMax") {
        mod.max.month <- raster::which.max(mod.clim.mly)
        ref.max.month <- raster::which.max(ref.clim.mly)
    }

    if (phaseMinMax == "phaseMin") {
        mod.max.month <- raster::which.min(mod.clim.mly)
        ref.max.month <- raster::which.min(ref.clim.mly)
    }

    # get shortest time distance between these months
    abs.diff <- abs(mod.max.month - ref.max.month)  # absolute difference from 0 to 12 months
    phase <- raster::calc(abs.diff, intFun.theta)  # shortest distance from 0 to 6 months (theta)
    phase.score <- 0.5 * (1 + cos(2 * pi * phase/12))  # score from 0 (6 months) to 1 (0 months)
    S_phase_not.weighted <- mean(raster::getValues(phase.score), na.rm = TRUE)  # scalar score (not weighted)
    # calculate the weighted scalar score
    a <- phase.score * weights  # this is a raster
    b <- raster::cellStats(a, "sum")  # this is a scalar, the sum up all values
    S_phase_weighted <- b/raster::cellStats(weights, "sum")  # scalar score (weighted)
    # compute global mean values of score input(s)
    mod.max.month.scalar <- mean(raster::getValues(mod.max.month), na.rm = TRUE)
    ref.max.month.scalar <- mean(raster::getValues(ref.max.month), na.rm = TRUE)
    phase.scalar <- mean(raster::getValues(phase), na.rm = TRUE)  # global mean value

    #---------------------------------------------------------------------------

    # (4) interannual variability

    #---------------------------------------------------------------------------

    years <- floor(raster::nlayers(mod)/12)  # total number of years
    months <- years * 12  # number of months considering complete years only
    mod.fullyear <- raster::subset(mod, 1:months)
    ref.fullyear <- raster::subset(ref, 1:months)
    c.mod <- raster::calc(mod.cycle, fun = function(x) {
        rep(x, years)
    })  # climatological cycle for all months (mod)
    c.ref <- raster::calc(ref.cycle, fun = function(x) {
        rep(x, years)
    })  # climatological cycle for all months (ref)
    mod.iav <- sqrt(raster::mean((mod.fullyear - c.mod)^2, na.rm = TRUE))  # interannual variability  (mod)
    ref.iav <- sqrt(raster::mean((ref.fullyear - c.ref)^2, na.rm = TRUE))  # interannual variability  (ref)
    # set values close to zero to NA
    ref.iav.na <- ref.iav
    ref.iav.na[ref.iav.na < 10^(-5)] <- NA
    epsilon_iav <- abs((mod.iav - ref.iav))/ref.iav.na
    epsilon_iav[epsilon_iav == Inf] <- NA  # I changed Eq. 26 so that epsilon_iav >= 0
    iav.score <- exp(-epsilon_iav)  # iav score as a function of space
    S_iav_not.weighted <- mean(raster::getValues(iav.score), na.rm = TRUE)  # scalar score (not weighted)
    # calculate the weighted scalar score
    a <- iav.score * weights  # this is a raster
    b <- raster::cellStats(a, "sum")  # this is a scalar, the sum up all values
    S_iav_weighted <- b/raster::cellStats(weights, "sum")  # scalar score (weighted)
    # compute global mean values of score input(s)
    mod.iav.scalar <- mean(raster::getValues(mod.iav), na.rm = TRUE)  # global mean value
    ref.iav.scalar <- mean(raster::getValues(ref.iav), na.rm = TRUE)  # global mean value
    # epsilon_iav.scalar <- stats::median(raster::getValues(epsilon_iav), na.rm = TRUE) # global mean value
    epsilon_iav.scalar <- stats::median(raster::getValues(epsilon_iav), na.rm = TRUE)  # global mean value

    #---------------------------------------------------------------------------

    # (5) dist

    #---------------------------------------------------------------------------

    mod.sigma.scalar <- stats::sd(raster::getValues(mod.mean), na.rm = TRUE)  # standard deviation of period mean data
    ref.sigma.scalar <- stats::sd(raster::getValues(ref.mean), na.rm = TRUE)  # standard deviation of period mean data
    sigma <- mod.sigma.scalar/ref.sigma.scalar
    y <- raster::getValues(mod.mean)
    x <- raster::getValues(ref.mean)
    reg <- stats::lm(y ~ x)
    R <- sqrt(summary(reg)$r.squared)
    S_dist <- 2 * (1 + R)/(sigma + 1/sigma)^2  # weighting does not apply

    #---------------------------------------------------------------------------

    # scores

    #---------------------------------------------------------------------------

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
    # get all score values in case you want to compare two runs having all values will enable you to conduct a significance test
    bias.score.values <- raster::getValues(bias.score)
    rmse.score.values <- raster::getValues(rmse.score)
    phase.score.values <- raster::getValues(phase.score)
    iav.score.values <- raster::getValues(iav.score)
    dist.score.values <- raster::getValues(mask * S_dist)
    all.score.values <- data.frame(bias.score.values, rmse.score.values, phase.score.values, iav.score.values, dist.score.values)
    colnames(all.score.values) <- c("bias.score", "rmse.score", "phase.score", "iav.score", "dist.score")
    all.score.values[is.na(all.score.values)] <- NA  # converts all NaN to NA
    if (outputDir != FALSE) {
        utils::write.table(all.score.values, paste(outputDir, "/", "allscorevalues", "-", variable.name, "-", ref.id, sep = ""))
    }
    # selected score inputs
    scoreinputs <- data.frame(long.name, variable.name, ref.id, variable.unit, mod.mean.scalar, ref.mean.scalar, bias.scalar,
        bias.scalar.rel, ref.sd.scalar, epsilon_bias.scalar, S_bias_not.weighted, rmse.scalar, crmse.scalar, ref.sd.scalar, epsilon_rmse.scalar,
        S_rmse_not.weighted, mod.max.month.scalar, ref.max.month.scalar, phase.scalar, S_phase_not.weighted, mod.iav.scalar,
        ref.iav.scalar, epsilon_iav.scalar, S_iav_not.weighted, mod.sigma.scalar, ref.sigma.scalar, sigma, R, S_dist)
    if (outputDir != FALSE) {
        utils::write.table(scoreinputs, paste(outputDir, "/", "scoreinputs", "_", meta.data.com, sep = ""))
    }
    # function that returns the min, max, and interval used in legend
    mmi.bias <- intFun.min.max.int.bias(bias)
    mmi.bias.score <- c(0, 1, 0.1)
    mmi.crmse <- intFun.min.max.int(crmse)
    mmi.rmse.score <- c(0, 1, 0.1)
    mmi.phase <- c(0, 6, 1)
    mmi.phase.score <- c(0, 1, 0.1)
    mmi.iav.score <- c(0, 1, 0.1)
    # intFun.min.max.int.mod.ref
    mmi.mean <- intFun.min.max.int.mod.ref(mod.mean, ref.mean)
    mmi.sd <- intFun.min.max.int.mod.ref(mod.sd, ref.sd)
    mmi.clim.mly <- intFun.min.max.int.mod.ref(mod.clim.mly, ref.clim.mly)
    mmi.bias.clim.mly <- intFun.min.max.int.bias(bias.clim.mly)
    mmi.max.month <- c(1, 12, 1)
    mmi.iav <- intFun.min.max.int.mod.ref(mod.iav, ref.iav)
    # add metadata: 1. filename (e.g. nee_mod.mean.nc), 2. figure title (e.g.  Mean_nee_ModID_123_from_1982-01_to_2008-12), 3.
    # min, max, interval used in legend (e.g. 0, 1, 0.1), 4.  legend bar text (e.g. 'score (-)')
    raster::metadata(mod.mean) <- list(paste(variable.name, ref.id, "mod_mean", sep = "_"), paste("Mean", meta.data.mod, sep = "_"),
        mmi.mean, variable.unit)
    raster::metadata(ref.mean) <- list(paste(variable.name, ref.id, "ref_mean", sep = "_"), paste("Mean", meta.data.ref, sep = "_"),
        mmi.mean, variable.unit)
    raster::metadata(bias) <- list(paste(variable.name, ref.id, "bias", sep = "_"), paste("Bias", meta.data.com, sep = "_"),
        mmi.bias, variable.unit)
    raster::metadata(bias.significance) <- list(paste(variable.name, ref.id, "bias_significance", sep = "_"), paste("Bias_significance",
        meta.data.com, sep = "_"))
    raster::metadata(bias.score) <- list(paste(variable.name, ref.id, "bias_score", sep = "_"), paste("Bias_score", meta.data.com,
        sep = "_"), mmi.bias.score, "score (-)")
    raster::metadata(crmse) <- list(paste(variable.name, ref.id, "crmse", sep = "_"), paste("CRMSE", meta.data.com, sep = "_"),
        mmi.crmse, variable.unit)
    raster::metadata(mod.sd) <- list(paste(variable.name, ref.id, "mod_sd", sep = "_"), paste("SD", meta.data.mod, sep = "_"),
        mmi.sd, variable.unit)
    raster::metadata(ref.sd) <- list(paste(variable.name, ref.id, "ref_sd", sep = "_"), paste("SD", meta.data.ref, sep = "_"),
        mmi.sd, variable.unit)
    raster::metadata(rmse.score) <- list(paste(variable.name, ref.id, "rmse_score", sep = "_"), paste("RMSE_score", meta.data.com,
        sep = "_"), mmi.rmse.score, "score (-)")
    raster::metadata(mod.clim.mly) <- list(paste(variable.name, ref.id, "mod_clim_mly", sep = "_"), paste("Monthly", meta.data.mod,
        sep = "_"), mmi.clim.mly, variable.unit)
    raster::metadata(ref.clim.mly) <- list(paste(variable.name, ref.id, "ref_clim_mly", sep = "_"), paste("Monthly", meta.data.ref,
        sep = "_"), mmi.clim.mly, variable.unit)
    raster::metadata(bias.clim.mly) <- list(paste(variable.name, ref.id, "bias_clim_mly", sep = "_"), paste("Monthly", meta.data.com,
        sep = "_"), mmi.bias.clim.mly, variable.unit)
    raster::metadata(mod.max.month) <- list(paste(variable.name, ref.id, "mod_max_month", sep = "_"), paste("Month_with_max",
        meta.data.mod, sep = "_"), mmi.max.month, "month")
    raster::metadata(ref.max.month) <- list(paste(variable.name, ref.id, "ref_max_month", sep = "_"), paste("Month_with_max",
        meta.data.ref, sep = "_"), mmi.max.month, "month")
    raster::metadata(phase) <- list(paste(variable.name, ref.id, "phase", sep = "_"), paste("Diff_in_max_month", meta.data.com,
        sep = "_"), mmi.phase, "month")
    raster::metadata(phase.score) <- list(paste(variable.name, ref.id, "phase_score", sep = "_"), paste("Seasonality_score",
        meta.data.com, sep = "_"), mmi.phase.score, "score (-)")
    raster::metadata(mod.iav) <- list(paste(variable.name, ref.id, "mod_iav", sep = "_"), paste("IAV", meta.data.mod, sep = "_"),
        mmi.iav, variable.unit)
    raster::metadata(ref.iav) <- list(paste(variable.name, ref.id, "ref_iav", sep = "_"), paste("IAV", meta.data.ref, sep = "_"),
        mmi.iav, variable.unit)
    raster::metadata(iav.score) <- list(paste(variable.name, ref.id, "iav_score", sep = "_"), paste("IAV_score", meta.data.com,
        sep = "_"), mmi.iav.score, "score (-)")
    # write data to netcdf
    stat.metric <- raster::stack(mod.mean, ref.mean, bias, bias.score, crmse, rmse.score, mod.max.month, ref.max.month, phase,
        phase.score, mod.iav, ref.iav, iav.score, mod.sd, ref.sd)

    # create a NetCDF file from stat.metric

    for (i in 1:raster::nlayers(stat.metric)) {
        data <- raster::subset(stat.metric, i:i)
        my.filename <- unlist(raster::metadata(data)[1])
        my.filename <- gsub("_", "-", my.filename)
        my.filename <- gsub(".", "-", my.filename, fixed = TRUE)
        my.longname <- unlist(raster::metadata(data)[2])
        if (outputDir != FALSE) {
            raster::writeRaster(data, filename = paste(outputDir, my.filename, sep = "/"), format = "CDF", varname = variable.name,
                longname = my.longname, varunit = variable.unit, overwrite = TRUE)
        }
    }

    # create NetCDF files for monthly mean climatologies
    mly.clim <- list(bias.clim.mly, mod.clim.mly, ref.clim.mly)
    for (i in 1:length(mly.clim)) {
        data <- mly.clim[[i]]
        my.filename <- unlist(raster::metadata(data)[1])
        my.filename <- gsub("_", "-", my.filename)
        my.filename <- gsub(".", "-", my.filename, fixed = TRUE)
        my.longname <- unlist(raster::metadata(data)[2])
        if (outputDir != FALSE) {
            raster::writeRaster(data, filename = paste(outputDir, my.filename, sep = "/"), format = "CDF", varname = variable.name,
                longname = my.longname, varunit = variable.unit, overwrite = TRUE)
        }
    }

    return(list(stat.metric, mod.outlier.points, ref.outlier.points, bias.significance, bias.clim.mly, mod.clim.mly, ref.clim.mly))
}
if (getRversion() >= "2.15.1") utils::globalVariables(c("getValues", "sd"))
