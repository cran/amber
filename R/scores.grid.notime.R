################################################################################ 
#' Scores for gridded reference data that do not have a varying time dimension
#' @description This function compares model output against remote-sensing
#' based reference data that do not vary in time. The performance of a model is
#' expressed through a score that ranges from zero to one, where increasing values
#' imply better performance.
#' Contrary to the function \link{scores.grid.time}, only one score is computed
#' (spatial distribution score, \eqn{S_{dist}}) since the reference data does
#' not vary with time:
#'
#' \eqn{\sigma=\sigma_{\overline{v_{mod}}}/\sigma_{\overline{v_{ref}}}}
#'
#' \eqn{S_{dist}=2(1+R)/(\sigma+\frac{1}{\sigma})^{2}}
#'
#' where \eqn{\sigma_{\overline{v_{mod}}}} and \eqn{\sigma_{\overline{v_{ref}}}} are the
#' standard deviation of the time mean values from the model and reference data,
#' and \eqn{R} is the spatial correlation coefficient of
#' \eqn{\overline{v_{ref}}(\lambda, \phi)} and
#' \eqn{\overline{v_{mod}}(\lambda, \phi)}.
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
#' @param outlier.factor A number that is used to define outliers, e.g. 10.
#'  Plotting raster objects that contain extreme outliers lead to figures where
#'  most grid cells are presented by a single color since the color legend covers
#'  the entire range of values. To avoid this, the user may define outliers that
#'  will be masked out and marked with a red dot. Outliers are all values that
#'  exceed the interquartile range multiplied by the outlier factor defined here.
#' @param irregular Logical. If TRUE the data are converted from an irregular to a regular grid. Default is FALSE.
#' @param my.projection A string that defines the projection of the irregular grid
#' @param numCores An integer that defines the number of cores, e.g. 2
#' @param period An R obect that gives the period over which to average the model data, e.g. c('1980-01', '2017-12')
#' @param timeInt A string that gives the time interval of the model data, e.g. 'month' or 'year'
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#'
#' @return (1) A list that contains three elements. The first element is a
#' raster stack with model data
#' (mean, \eqn{mod.mean}),
#' reference data
#' (mean, \eqn{ref.mean}),
#' and the corresponding bias
#' (bias, \eqn{bias}). The second and third element of the list are spatial
#' point data frames that give the model and reference outliers, respectively.
#' The content of the list can be plotted using \link{plotGrid}.
#'
#' (2) NetCDF files for each of the statistical variables listed above.
# 
#' (3) Three text files: (i) score values and (ii) score inputs averaged across
#' the entire grid, and (iii) score values for each individual grid cell.
#'
#' @examples
#'
#' library(amber)
#' library(doParallel)
#' library(foreach)
#' library(ncdf4)
#' library(parallel)
#' library(raster)
#'
#' # (1) Global plots on a regular grid
#' long.name <- 'Soil Carbon'
#' nc.mod <- system.file('extdata/modelRegular', 'cSoil_monthly.nc', package = 'amber')
#' nc.ref <- system.file('extdata/referenceRegular', 'soilc_HWSD_128x64.nc', package = 'amber')
#' mod.id <- 'CLASSIC' # define a model experiment ID
#' ref.id <- 'HWSD' # give reference dataset a name
#' unit.conv.mod <- 1 # optional unit conversion for model data
#' unit.conv.ref <- 1 # optional unit conversion for reference data
#' variable.unit <- 'kgC m$^{-2}$' # unit after conversion (LaTeX notation)
#'
#' # Short version using default settings:
#' plot.me <- scores.grid.notime(long.name, nc.mod, nc.ref, mod.id, ref.id,
#' unit.conv.mod, unit.conv.ref, variable.unit)
#' plotGrid(long.name, plot.me)
#'
#' \donttest{
#' # Additional parameters:
#' score.weights <- c(1,2,1,1,1) # score weights of S_bias, S_rmse, S_phase, S_iav, S_dist
#' outlier.factor <- 1000
#' irregular <- FALSE
#' my.projection <- '+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.'
#' numCores <- 2
#' period <- c('1980-01', '2017-12') # period over which to average the model data
#'
#' plot.me <- scores.grid.notime(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit, score.weights, outlier.factor, irregular, my.projection,
#' numCores, period)
#' plotGrid(long.name, plot.me)
#'
#'
#' # (2) Regional plots on a rotated grid
#' long.name <- 'Soil Carbon'
#' nc.mod <- system.file('extdata/modelRotated', 'cSoil_monthly.nc', package = 'amber')
#' nc.ref <- system.file('extdata/referenceRotated', 'soilc_HWSD_rotated.nc', package = 'amber')
#' mod.id <- 'CLASSIC' # define a model experiment ID
#' ref.id <- 'HWSD' # give reference dataset a name
#' unit.conv.mod <- 1 # optional unit conversion for model data
#' unit.conv.ref <- 1 # optional unit conversion for reference data
#' variable.unit <- 'kgC m$^{-2}$' # unit after conversion (LaTeX notation)
#' score.weights <- c(1,2,1,1,1) # score weights of S_bias, S_rmse, S_phase, S_iav, S_dist
#' outlier.factor <- 1000
#' irregular <- TRUE
#' my.projection <- '+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.'
#' numCores <- 2
#' period <- c('1980-01', '2017-12') # period over which to average the model data
#'
#' plot.me <- scores.grid.notime(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit, score.weights, outlier.factor, irregular, my.projection,
#' numCores, period)
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
scores.grid.notime <- function(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod, unit.conv.ref, variable.unit, score.weights = c(1, 
    2, 1, 1, 1), outlier.factor = 1000, irregular = FALSE, my.projection = "+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.", 
    numCores = 2, period = c("1980-01", "2017-12"), timeInt = "month", outputDir = FALSE) {
    
    #---------------------------------------------------------------------------
    
    # (I) Data preparation
    
    #---------------------------------------------------------------------------
    
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
                if (nTime > 1) {
                  myValues <- data[, , i]
                } else {
                  myValues <- data
                }  # first case: model (with time), second case: ref (no time)
                myValues <- base::apply(base::t(myValues), 2, rev)  # rotate values
                r <- raster::setValues(r, myValues)  # assign values
                raster::extent(r) <- myExtent  # extent using the rotated projection
                r <- r * 1  # this is necessary for the base::do.call function below
            }
            myStack <- base::do.call(raster::stack, myRaster)
            parallel::stopCluster(cl)
            myStack <- raster::setZ(myStack, dates, name = "time")
            names(myStack) <- dates
            assign(c("mod", "ref")[id], myStack)
        }
        # model data are averaged over a period since the reference data set has no time dimension
        start.date <- period[1]
        end.date <- period[2]
        
        mod <- mod[[which(format(as.Date(raster::getZ(mod)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(mod)), 
            "%Y-%m") <= end.date)]]
        mod <- raster::mean(mod, na.rm = TRUE)
        
        nc <- ncdf4::nc_open(nc.mod)
        variable.name <- names(nc[["var"]])
        ncdf4::nc_close(nc)
        variable.name <- variable.name[1]
        variable.name <- ifelse(variable.name == "burntFractionAll", "burnt", variable.name)  # rename burntFractionAll to shorter name
        variable.name <- toupper(variable.name)  # make variable name upper-case
        
        # unit conversion if appropriate
        mod <- mod * unit.conv.mod
        ref <- ref * unit.conv.ref
        
        # get layer names
        mod.names <- names(mod)
        ref.names <- names(ref)
    } else {
        
        #-----------------------------------------------------------------------
        
        # (2) Process data if 'irregular=FALSE'
        
        #-----------------------------------------------------------------------
        
        # model data are averaged over a period since the reference data set has no time dimension
        start.date <- period[1]
        end.date <- period[2]
        
        mod <- raster::brick(nc.mod)
        mod <- raster::rotate(mod)
        dates.mod <- raster::getZ(mod)
        
        mod <- raster::setZ(mod, dates.mod, name = "date")
        mod <- mod[[which(format(as.Date(raster::getZ(mod)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(mod)), 
            "%Y-%m") <= end.date)]]
        mod <- raster::mean(mod, na.rm = TRUE)
        
        nc <- ncdf4::nc_open(nc.mod)
        variable.name <- names(nc[["var"]])
        ncdf4::nc_close(nc)
        variable.name <- variable.name[length(variable.name)]  # take the last variable (relevant for CanESM5)
        variable.name <- ifelse(variable.name == "burntFractionAll", "burnt", variable.name)  # rename burntFractionAll to shorter name
        variable.name <- toupper(variable.name)  # make variable name upper-case
        
        # reference data
        ref <- raster::raster(nc.ref)  # this reference data has no time dimension
        ref <- raster::rotate(ref)
        # unit conversion if appropriate
        mod <- mod * unit.conv.mod
        ref <- ref * unit.conv.ref
        
        # get layer names
        mod.names <- names(mod)
        ref.names <- names(ref)
    }
    
    #---------------------------------------------------------------------------
    
    # 1.3 Remaining part applies to both regular and irregular gridded data
    
    #---------------------------------------------------------------------------
    
    # Make a string that summarizes metadata. This will be added to each netcdf file (longname).
    
    # The string can then be accessed like this: names(raster(file.nc))
    
    meta.data.mod <- paste(variable.name, mod.id, "from", start.date, "to", end.date, sep = "_")
    meta.data.ref <- paste(variable.name, ref.id, "from", start.date, "to", end.date, sep = "_")
    meta.data.com <- paste(variable.name, mod.id, "vs", ref.id, "from", start.date, "to", end.date, sep = "_")
    # outliers all extreme outliers are set to NA in the grid
    mod.mean <- mod  # time mean (mod is already a time-mean)
    mod.outlier_range <- intFun.grid.define.outlier(mod.mean, outlier.factor)  # define outlier range
    outlier.neg <- mod.outlier_range[1]
    outlier.pos <- mod.outlier_range[2]
    mod.mask_outliers <- intFun.grid.outliers.na(mod.mean, outlier.neg, outlier.pos)
    mod.mask_outliers <- mod.mask_outliers - mod.mask_outliers + 1
    mod <- mod * mod.mask_outliers
    names(mod) <- mod.names
    mod.outlier.points <- intFun.grid.outliers.points(mod.mean, outlier.neg, outlier.pos)
    ## reference data
    ref.mean <- ref  # time mean (ref is already a time-mean due to preprocessing)
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
    mod.mean <- mod  # data is already time mean
    ref.mean <- ref  # data is already time mean
    bias <- mod.mean - ref.mean  # time mean
    S_bias_not.weighted <- NA
    S_bias_weighted <- NA
    # Compute global mean values of score input(s) The mask ensures that mod and ref are based on same grid cells.
    mask <- (mod.mean * ref.mean)
    mask <- mask - mask + 1
    mod.mean.scalar <- mean(raster::getValues(mask * mod.mean), na.rm = TRUE)  # global mean value
    ref.mean.scalar <- mean(raster::getValues(mask * ref.mean), na.rm = TRUE)  # global mean value
    bias.scalar <- mean(raster::getValues(bias), na.rm = TRUE)  # global mean value
    ref.sd.scalar <- NA  # global mean value
    epsilon_bias.scalar <- NA  # global mean value
    
    #---------------------------------------------------------------------------
    
    # (5) dist
    
    #---------------------------------------------------------------------------
    mod.sigma.scalar <- sd(raster::getValues(mod.mean), na.rm = TRUE)  # standard deviation of period mean data
    ref.sigma.scalar <- sd(raster::getValues(ref.mean), na.rm = TRUE)  # standard deviation of period mean data
    sigma <- mod.sigma.scalar/ref.sigma.scalar
    y <- raster::getValues(mod.mean)
    x <- raster::getValues(ref.mean)
    reg <- stats::lm(y ~ x)
    R <- sqrt(summary(reg)$r.squared)
    S_dist <- 2 * (1 + R)/(sigma + 1/sigma)^2  # weighting does not apply
    
    #---------------------------------------------------------------------------
    
    # Scores
    
    #---------------------------------------------------------------------------
    
    S_bias <- S_bias_not.weighted
    S_rmse <- NA
    S_phase <- NA
    S_iav <- NA
    # weight importance of statisitcal metrics and compute overall score
    S_overall <- S_dist
    scores <- data.frame(variable.name, ref.id, S_bias, S_rmse, S_phase, S_iav, S_dist, S_overall)
    scores_not.weighted <- scores
    # weighted (except for S_dist)
    S_bias <- S_bias_weighted
    S_rmse <- NA
    S_phase <- NA
    S_iav <- NA
    S_overall <- S_dist
    scores <- data.frame(variable.name, ref.id, S_bias, S_rmse, S_phase, S_iav, S_dist, S_overall)
    scores_weighted <- scores
    # 
    scores <- rbind(scores_not.weighted, scores_weighted)
    rownames(scores) <- c("not.weighted", "weighted")
    if (outputDir != FALSE) {
        utils::write.table(scores, paste(outputDir, "/", "scorevalues", "_", meta.data.com, sep = ""))
    }
    # get all score values in case you want to compare this run against another run having all values will enable you to
    # conduct a significance test
    dist.score <- raster::getValues(mask * S_dist)
    all.score.values <- data.frame(NA, NA, NA, NA, dist.score)
    colnames(all.score.values) <- c("bias.score", "rmse.score", "phase.score", "iav.score", "dist.score")
    all.score.values[is.na(all.score.values)] <- NA  # converts all NaN to NA
    if (outputDir != FALSE) {
        utils::write.table(all.score.values, paste(outputDir, "/", "allscorevalues", variable.name, "-", ref.id, sep = ""))
    }
    # selected score inputs
    rmse.scalar <- NA
    crmse.scalar <- NA
    phase.scalar <- NA
    mod.iav.scalar <- NA
    ref.iav.scalar <- NA
    epsilon_rmse.scalar <- NA
    S_rmse_not.weighted <- NA
    mod.max.month.scalar <- NA
    ref.max.month.scalar <- NA
    S_phase_not.weighted <- NA
    epsilon_iav.scalar <- NA
    S_iav_not.weighted <- NA
    
    scoreinputs <- data.frame(long.name, variable.name, ref.id, variable.unit, mod.mean.scalar, ref.mean.scalar, bias.scalar, 
        ref.sd.scalar, epsilon_bias.scalar, S_bias_not.weighted, rmse.scalar, crmse.scalar, ref.sd.scalar, epsilon_rmse.scalar, 
        S_rmse_not.weighted, mod.max.month.scalar, ref.max.month.scalar, phase.scalar, S_phase_not.weighted, mod.iav.scalar, 
        ref.iav.scalar, epsilon_iav.scalar, S_iav_not.weighted, mod.sigma.scalar, ref.sigma.scalar, sigma, R, S_dist)
    if (outputDir != FALSE) {
        utils::write.table(scoreinputs, paste(outputDir, "/", "scoreinputs", "_", meta.data.com, sep = ""))
    }
    # function that returns the min, max, and interval used in legend function min.max.int only requires one input apply the
    # functions min.max.int and min.max.int.mod.ref min.max.int
    mmi.bias <- intFun.min.max.int.bias(bias)
    mmi.mean <- intFun.min.max.int.mod.ref(mod.mean, ref.mean)
    # add metadata: 1. filename (e.g. nee_mod.mean.nc), 2. figure title (e.g.  Mean_nee_ModID_123_from_1982-01_to_2008-12),
    # 3.  min, max, interval used in legend (e.g. 0, 1, 0.1), 4.  legend bar text (e.g. 'score (-)')
    raster::metadata(mod.mean) <- list(paste(variable.name, ref.id, "mod_mean", sep = "_"), paste("Mean", meta.data.mod, 
        sep = "_"), mmi.mean, variable.unit)
    raster::metadata(ref.mean) <- list(paste(variable.name, ref.id, "ref_mean", sep = "_"), paste("Mean", meta.data.ref, 
        sep = "_"), mmi.mean, variable.unit)
    raster::metadata(bias) <- list(paste(variable.name, ref.id, "bias", sep = "_"), paste("Bias", meta.data.com, sep = "_"), 
        mmi.bias, variable.unit)
    
    # create a netcdf file from stat.metric
    
    stat.metric <- raster::stack(mod.mean, ref.mean, bias)
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
    return(list(stat.metric, mod.outlier.points, ref.outlier.points))
}
