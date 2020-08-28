################################################################################
#' Scores for site-level reference data that do not vary in time
#' @description This function compares model output against
#' site-level measurements such as carbon stocks. The performance of a model is
#' expressed through scores that range from zero to one, where increasing values
#' imply better performance.
#' Contrary to the function \link{scores.grid.time}, only two scores are computed
#' (bias score \eqn{S_{bias}} and spatial distribution score, \eqn{S_{dist}}) since the reference data do
#' not vary with time. Contrary to \link{scores.grid.time}, the bias is relative to the absolute reference mean
#' value rather than the reference standard deviation. Again, this is because the reference data do
#' not vary with time:
#'
#' \eqn{(i) \ bias(\lambda, \phi)=\overline{v_{mod}}(\lambda, \phi)-\overline{v_{ref}}(\lambda, \phi)}
#'
#' \eqn{(ii) \ \varepsilon_{bias}=|bias(\lambda, \phi)|/|\overline{v_{ref}}(\lambda, \phi)|}
#'
#' \eqn{(iii) \ s_{bias}(\lambda, \phi)=e^{-\varepsilon_{bias}(\lambda, \phi)}}
#'
#' \eqn{(iv) \ S_{bias}=\overline{\overline{s_{bias}}}}
#'
#' @param long.name A string that gives the full name of the variable, e.g. 'Gross primary productivity'
#' @param nc.mod A string that gives the path and name of the netcdf file that contains the model output, e.g. '/home/model_gpp.nc'
#' @param ref.csv A string that gives the path and name of the csv file that contains the reference data output, e.g. '/home/reference_biomass.csv'.
#' The columns must be in the following order: Plot ID, longitude, latitude, data values.
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
#' @param period An R obect that gives the period over which to average the model data, e.g. c('1980-01', '2017-12')
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @param variable.name A string with the variable name, e.g. 'GPP'. If FALSE, the variable name stored in the NetCDF file will be used instead. Default is FALSE.
#' @param meanPerGridCell Logical. If TRUE, then values from different sites that
#' are located in the same grid cell are averaged. Default is set to TRUE.
#' @param myCex A number that determines the size of the dots in the Figure. Default is set to 0.7.
#' @param subcaption A string that defines the subcaption of the figure, e.g. '(a)'.
#' @return (1) Figures in PDF format that show maps of the model mean, reference mean, and bias.
#' (2) Four text files: (i) score values and (ii) score inputs for each individual
#' site, and (iii) score values and (iv) score inputs averaged across sites.
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
#' long.name <- 'soil carbon'
#' nc.mod <- system.file('extdata/modelRegular', 'cSoil_monthly.nc', package = 'amber')
#' ref.csv <- system.file('extdata/siteLevelRefData', 'siteLevelDataNoTime.csv', package = 'amber')
#' mod.id <- 'CLASSIC' # define a model experiment ID
#' ref.id <- 'ABC' # give reference dataset a name
#' unit.conv.mod <- 1 # optional unit conversion for model data
#' unit.conv.ref <- 1 # optional unit conversion for reference data
#' variable.unit <- 'kgC m$^{-2}$' # unit after conversion (LaTeX notation)
#'
#' # Short version using default settings:
#' scores.site.notime(long.name, nc.mod, ref.csv, mod.id, ref.id,
#' unit.conv.mod, unit.conv.ref, variable.unit)
#'
#' # To zoom into a particular region:
#' scores.site.notime(long.name, nc.mod, ref.csv, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit, score.weights = c(1, 2, 1, 1, 1),
#' my.xlim = c(-150, -60), my.ylim = c(20, 60), plot.width = 6, plot.height = 3.8)
#'
#' # (2) Regional plots on a rotated grid
#' nc.mod <- system.file('extdata/modelRotated', 'cSoil_monthly.nc', package = 'amber')
#' ref.csv <- system.file('extdata/siteLevelRefData', 'siteLevelDataNoTime.csv', package = 'amber')
#' mod.id <- 'CLASSIC' # define a model experiment ID
#' ref.id <- 'ABC' # give reference dataset a name
#' unit.conv.mod <- 1 # optional unit conversion for model data
#' unit.conv.ref <- 1 # optional unit conversion for reference data
#' variable.unit <- 'kgC m$^{-2}$' # unit after conversion (LaTeX notation)
#' rotate.me <- FALSE
#' irregular <- TRUE
#' my.projection <-'+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.'
#' # shp.filename <- system.file('extdata/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp',
#' #  package = 'amber')
#' shp.filename <- system.file("extdata/ne_110m_land/ne_110m_land.shp", package = "amber")
#' my.xlim <- c(-171, 0) # longitude range that you wish to plot
#' my.ylim <- c(32, 78) # latitude range that you wish to plot
#' plot.width <- 7
#' plot.height <- 3.8
#' numCores <- 2
#'
#' scores.site.notime(long.name, nc.mod, ref.csv, mod.id, ref.id,
#' unit.conv.mod, unit.conv.ref, variable.unit, rotate.me = TRUE, irregular = TRUE,
#' my.projection = my.projection, shp.filename = shp.filename,
#' my.xlim = my.xlim,  my.ylim = my.ylim, plot.width = plot.width, plot.height = plot.height)
#' } #donttest
#'
#'
#' @export
scores.site.notime <- function(long.name, nc.mod, ref.csv, mod.id, ref.id, unit.conv.mod, unit.conv.ref, variable.unit, score.weights = c(1,
    2, 1, 1, 1), rotate.me = TRUE, irregular = FALSE, my.projection = "+proj=longlat +ellps=WGS84", shp.filename = system.file("extdata/ne_110m_land/ne_110m_land.shp",
    package = "amber"), my.xlim = c(-180, 180), my.ylim = c(-60, 85), plot.width = 8, plot.height = 3.8, numCores = 2, period = c("1980-01",
    "2017-12"), outputDir = FALSE, variable.name = FALSE, meanPerGridCell = TRUE, myCex = 0.5, subcaption = "") {

    #---------------------------------------------------------------------------

    # (I) Data preparation

    #---------------------------------------------------------------------------

    # (1) Reproject data from an irregular to a regular grid if 'irregular=TRUE'

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
        mod <- myStack

        # model data are averaged over a period since the reference data set has no time dimension
        start.date <- period[1]
        end.date <- period[2]

        dates.mod <- raster::getZ(mod)

        mod <- raster::setZ(mod, dates.mod, name = "date")
        mod <- mod[[which(format(as.Date(raster::getZ(mod)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(mod)), "%Y-%m") <=
            end.date)]]

        mod <- mod * unit.conv.mod
        mod.mean <- raster::mean(mod, na.rm = TRUE)
        # LAI measurements reflect maximum values. Therefore, use max rather than mean values:
        if (variable.name == "LAI") {
            index <- format(as.Date(names(mod), format = "X%Y.%m.%d"), format = "%m")
            index <- as.numeric(index)
            mod.clim.mly <- raster::stackApply(mod, index, fun = raster::mean)
            mod.mean <- max(mod.clim.mly, na.rm = TRUE)
        }

        #-----------------------------------------------------------------------

        # Process reference data

        #-----------------------------------------------------------------------

        ref <- utils::read.csv(ref.csv)
        ref[ref == -9999] <- NA
        lon <- ref[[2]]
        lat <- ref[[3]]
        ref <- ref[[4]]
        ref <- ref * unit.conv.ref

        # reproject coordinates to rotated grid
        lon <- base::matrix(lon, ncol = 1)
        lat <- base::matrix(lat, ncol = 1)
        lonLat <- base::data.frame(lon, lat)
        sp::coordinates(lonLat) <- ~lon + lat
        raster::projection(lonLat) <- regular
        suppressWarnings(lonLat <- sp::spTransform(lonLat, sp::CRS(rotated)))
        lon <- sp::coordinates(lonLat)[, 1]
        lat <- sp::coordinates(lonLat)[, 2]

        ref <- intFun.site.points(lon, lat, ref)

        # compute statistics for each site

        mod <- raster::extract(mod.mean, ref, method = "bilinear")

        mod <- data.frame(mod)

        # rownames(mod) <- site.id
        ref <- as.data.frame(ref)
        colnames(ref) <- c("lon", "lat", "ref")
        ref <- ref[-c(1, 2)]  # drop lon lat column

    } else {

        #-----------------------------------------------------------------------

        # (2) Process data if 'irregular=FALSE'

        #-----------------------------------------------------------------------
        myExtent <- c(-180, 180, -90, 90)  # this can be done in a smarter way

        # model data are averaged over a period since the reference data set has no time dimension
        start.date <- period[1]
        end.date <- period[2]

        mod <- raster::brick(nc.mod)

        if (rotate.me == TRUE) {
          suppressWarnings(mod <- raster::rotate(mod))
        }
        dates.mod <- raster::getZ(mod)

        mod <- raster::setZ(mod, dates.mod, name = "date")
        mod <- mod[[which(format(as.Date(raster::getZ(mod)), "%Y-%m") >= start.date & format(as.Date(raster::getZ(mod)), "%Y-%m") <=
            end.date)]]

        mod <- mod * unit.conv.mod
        mod.mean <- raster::mean(mod, na.rm = TRUE)

        # model data
        if (variable.name == FALSE) {
            nc <- ncdf4::nc_open(nc.mod)
            variable.name <- names(nc[["var"]])
            ncdf4::nc_close(nc)
            variable.name <- variable.name[length(variable.name)]  # take the last variable (relevant for CanESM5)
            variable.name <- ifelse(variable.name == "burntFractionAll", "burnt", variable.name)  # rename burntFractionAll to shorter name
            variable.name <- toupper(variable.name)  # make variable name upper-case
        }

        # LAI measurements reflect maximum values. Therefore, use max rather than mean values:
        if (variable.name == "LAI") {
            index <- format(as.Date(names(mod), format = "X%Y.%m.%d"), format = "%m")
            index <- as.numeric(index)
            mod.clim.mly <- raster::stackApply(mod, index, fun = raster::mean)
            mod.mean <- max(mod.clim.mly, na.rm = TRUE)
        }

        # Process reference data

        ref <- utils::read.csv(ref.csv)
        ref[ref == -9999] <- NA
        lon <- ref[[2]]
        lat <- ref[[3]]
        ref <- ref[[4]]
        ref <- ref * unit.conv.ref
        ref <- intFun.site.points(lon, lat, ref)

        if (meanPerGridCell == TRUE) {

            # get cell number and coordinates for each grid cell
            myGrid <- mod.mean
            myGrid <- raster::setValues(myGrid, 1)  # this ensures that all grid cells are considered, including ocean grid cells
            modLonLat <- raster::rasterToPoints(myGrid, spatial = TRUE)
            modcn <- raster::extract(myGrid, modLonLat, cellnumbers = TRUE)
            lonLatCN <- data.frame(sp::coordinates(modLonLat), modcn[, 1])
            colnames(lonLatCN) <- c("lon", "lat", "Group.1")

            # get cell number for each site
            cn <- raster::extract(mod, ref, cellnumbers = TRUE)
            cn <- data.frame(cn[, 1])
            colnames(cn) <- "Group.1"

            # add cell numbers to reference data
            ref.cn <- data.frame(cn, ref)

            # compute mean value per grid cell (mpg)
            index <- list(ref.cn$Group.1)
            mpg <- stats::aggregate(ref$data, by = index, FUN = mean, na.rm = TRUE)

            # add the grid cell coordinates that correspond to each site
            mpg <- merge(mpg, lonLatCN, by = "Group.1")
            lon <- mpg$lon
            lat <- mpg$lat

            mpg <- subset(mpg, select = -c(Group.1, lon, lat))
            ref <- intFun.site.points(lon, lat, mpg)
        }

        # compute statistics for each site

        mod <- raster::extract(mod.mean, ref, method = "simple")
        mod <- data.frame(mod)

        # rownames(mod) <- site.id
        ref <- as.data.frame(ref)
        colnames(ref) <- c("lon", "lat", "ref")
        ref <- ref[-c(1, 2)]  # drop lon lat column

    }

    # check whether number of columns and rows match between mod and ref
    stopifnot(ncol(ref) == ncol(mod))
    stopifnot(nrow(ref) == nrow(mod))
    # make a string that summarizes metadata
    meta.data.mod <- paste(variable.name, mod.id, "from", start.date, "to", end.date, sep = "_")
    meta.data.ref <- paste(variable.name, ref.id, "from", start.date, "to", end.date, sep = "_")
    meta.data.com <- paste(variable.name, mod.id, "vs", ref.id, "from", start.date, "to", end.date, sep = "_")
    #---------------------------------------------------------------------------

    # (II) Statistical analysis

    #---------------------------------------------------------------------------

    # The time period of the reference data varies among sites.

    # Some sites are located in the model ocean.

    # Exclude data that do not have a corresponding data set in mod or ref.

    mask.mod <- (mod - mod + 1)
    mask.ref <- (ref - ref + 1)
    mask <- mask.mod * mask.ref
    mod <- mod * mask
    ref <- ref * mask

    #---------------------------------------------------------------------------

    # (1) bias

    #---------------------------------------------------------------------------
    mod.mean <- mod$mod  # already time mean
    ref.mean <- ref$ref  # already time mean
    weights <- ref.mean  # weights used for spatial integral
    bias <- mod.mean - ref.mean  # time mean
    epsilon_bias <- abs(bias)/abs(ref.mean)
    epsilon_bias[epsilon_bias == Inf] <- NA  # relative error
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
    bias.scalar.rel <- (mod.mean.scalar - ref.mean.scalar)/abs(ref.mean.scalar) * 100
    ref.sd.scalar <- NA
    epsilon_bias.scalar <- stats::median(epsilon_bias, na.rm = TRUE)

    #---------------------------------------------------------------------------

    # (5) dist

    #---------------------------------------------------------------------------

    mod.sigma.scalar <- sd(mod.mean, na.rm = TRUE)  # standard deviation of period mean data
    ref.sigma.scalar <- sd(ref.mean, na.rm = TRUE)  # standard deviation of period mean data
    sigma <- mod.sigma.scalar/ref.sigma.scalar
    y <- mod.mean
    x <- ref.mean
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
    S_overall <- (S_bias + S_dist)/2
    scores <- data.frame(variable.name, ref.id, S_bias, S_rmse, S_phase, S_iav, S_dist, S_overall)
    scores_not.weighted <- scores
    # weighted (except for S_dist)
    S_bias <- S_bias_weighted
    S_rmse <- NA
    S_phase <- NA
    S_iav <- NA
    S_overall <- (S_bias + S_dist)/2

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
    all.score.values <- data.frame(bias.score, NA, NA, NA, dist.score)
    colnames(all.score.values) <- c("bias.score", "rmse.score", "phase.score", "iav.score", "dist.score")
    all.score.values[is.na(all.score.values)] <- NA  # converts all NaN to NA
    if (outputDir != FALSE) {
        utils::write.table(all.score.values, paste(outputDir, "/", "allscorevalues", "-", variable.name, "-", ref.id, sep = ""))
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
        bias.scalar.rel, ref.sd.scalar, epsilon_bias.scalar, S_bias_not.weighted, rmse.scalar, crmse.scalar, ref.sd.scalar, epsilon_rmse.scalar,
        S_rmse_not.weighted, mod.max.month.scalar, ref.max.month.scalar, phase.scalar, S_phase_not.weighted, mod.iav.scalar,
        ref.iav.scalar, epsilon_iav.scalar, S_iav_not.weighted, mod.sigma.scalar, ref.sigma.scalar, sigma, R, S_dist)
    if (outputDir != FALSE) {
        utils::write.table(scoreinputs, paste(outputDir, "/", "scoreinputs", "_", meta.data.com, sep = ""))
    }
    # min.max.int
    mmi.bias <- intFun.min.max.int.bias(bias)
    mmi.bias.score <- c(0, 1, 0.1)
    mmi.mean <- intFun.min.max.int.mod.ref(mod.mean, ref.mean)

    #---------------------------------------------------------------------------

    # Add metadata:

    # 1. figure title (e.g. Mean_nee_ModID_123_from_1982-01_to_2008-12)

    # 2. min, max, interval used in legend (e.g. 0, 1, 0.1)

    # 3. legend bar text (e.g. 'score (-)')
    attr(mod.mean, "metadata") <- list(paste("Mean", meta.data.mod, sep = "_"), mmi.mean, variable.unit)
    attr(ref.mean, "metadata") <- list(paste("Mean", meta.data.ref, sep = "_"), mmi.mean, variable.unit)
    attr(bias, "metadata") <- list(paste("Bias", meta.data.com, sep = "_"), mmi.bias, variable.unit)
    attr(bias.score, "metadata") <- list(paste("Bias_score", meta.data.com, sep = "_"), mmi.bias.score, "score (-)")

    # write data to file
    stat.metric <- data.frame(lon, lat, mod.mean, ref.mean, bias, bias.score)
    colnames(stat.metric) <- c("lon", "lat", "mod.mean", "ref.mean", "bias", "bias.score")
    my.filename <- paste(variable.name, ref.id, sep = "_")
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
        my.breaks.bias <- round(seq(min - interval, max + interval, interval), 3)  # location of color breaks
        my.labels <- round(seq(min, max, interval), 3)  # locations where to set the labels
        my.col <- viridis::viridis(n = length(my.breaks) - 1, direction = -1)
        my.col.bias <- scico::scico(n = length(my.breaks) - 1, palette = "vik")
        my.col.bias <- c("magenta", my.col.bias, "black")  # add colors for outliers
        my.col.phase <- grDevices::rainbow(n = length(my.breaks) - 1)
        if (i == 5)
            {
                my.col <- my.col.bias
                my.breaks <- my.breaks.bias
                # Change values of outliers for the purpose of the legend
                data <- data[[1]]
                data[data > max] <- max + interval/2
                data[data < min] <- min - interval/2
                values <- data.frame(data)
                colnames(values) <- "values"
                # data <- intFun.site.points(lon, lat, values) # make spatial points

                data <- data.frame(lon, lat, values)
                bias.pos <- subset(data, data$values >= 0)
                bias.neg <- subset(data, data$values < 0)

                bias.pos <- intFun.site.points(bias.pos$lon, bias.pos$lat, bias.pos$values)
                bias.neg <- intFun.site.points(bias.neg$lon, bias.neg$lat, bias.neg$values)

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

        if (i == 5) {
            colors.pos <- cut(bias.pos$data, breaks = my.breaks, labels = my.col, include.lowest = TRUE)
            colors.pos <- toString(colors.pos)
            colors.pos <- unlist(strsplit(colors.pos, split = ", "))

            colors.neg <- cut(bias.neg$data, breaks = my.breaks, labels = my.col, include.lowest = TRUE)
            colors.neg <- toString(colors.neg)
            colors.neg <- unlist(strsplit(colors.neg, split = ", "))
        }

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
        raster::plot(dummy, col = NA, xlim = my.xlim, ylim = my.ylim, main = paste(paste(subcaption, long.name, sep = " "), my.title,
            sep = "\n"), axes = FALSE, legend = FALSE)
        raster::plot(land, col = "grey85", border = NA, add = TRUE)

        # Show positive and negative biases with different symbols
        if (i == 5) {
            graphics::points(bias.pos, pch = 24, bg = colors.pos, cex = myCex, lwd = 0.1)
            graphics::points(bias.neg, pch = 25, bg = colors.neg, cex = myCex, lwd = 0.1)
        }

        # all other metrics use same symbol
        if (i != 5) {
            graphics::points(data, pch = 21, bg = colors, cex = myCex * 1.2, lwd = 0.1)
        }

        # add values to plot my.pointLabel <- toString(round(data$values, 1)) my.pointLabel <- unlist(strsplit(my.pointLabel, ', '))
        # maptools::pointLabel(data$lon, data$lat, labels = my.pointLabel, cex = 0.5, col = 'black')

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
utils::globalVariables(c("%dopar%", "Group.1", "variable.name"))

