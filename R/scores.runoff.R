################################################################################
#' Scores for runoff
#' @description This function compares modelled runoff and measured streamflow
#' on an annual basis for 50 river basins. Modelled runoff is converted to
#' streamflow by spatially integrating annual runoff across each river basins.
#' The performance of a model is expressed through scores that range from
#' zero to one, where increasing values imply better performance. These scores
#' are computed in five steps:
#' \eqn{(i)} computation of a statistical metric,
#' \eqn{(ii)} nondimensionalization,
#' \eqn{(iii)} conversion to unit interval,
#' \eqn{(iv)} spatial integration, and
#' \eqn{(v)} averaging scores computed from different statistical metrics.
#' The latter includes the bias, root-mean-square error,
#' inter-annual variability, and spatial distribution. The corresponding equations
#' are documented in \code{\link{amber-package}}. Contrary to the function
#' \link{scores.grid.time}, no phase score is computed since data is evaluated on
#' an annual basis. Assessing annual rather than monthly values causes the standard
#' deviation \eqn{sigma_{ref}(\lambda, \phi)} to be very small. For this reason,
#' the relative errors for runoff are computed with repect to the mean rather than
#' the standard deviation:
#'
#' \eqn{(ii) \ \varepsilon_{bias}=|bias(\lambda, \phi)|/\overline{v_{ref}}(\lambda, \phi)}
#' \eqn{(ii) \ \varepsilon_{rmse}(\lambda, \phi)=crmse(\lambda, \phi)/\overline{v_{ref}}(\lambda, \phi).}
#'
#' @param long.name A string that gives the full name of the variable, e.g. 'Gross primary productivity'
#' @param nc.mod A string that gives the path and name of the netcdf file that contains the model output, e.g. '/home/model_gpp.nc'
#' @param nc.ref A string that gives the path and name of the netcdf file that contains the reference data output, e.g. '/home/reference_gpp.nc'
#' @param nc.basins A string that gives the path and name of the netcdf file that contains the river basins, e.g. '/home/basins.nc'
#' @param mod.id A string that identifies the source of the reference data set, e.g. 'CanESM2'
#' @param ref.id A string that identifies the source of the reference data set, e.g. 'MODIS'
#' @param unit.conv.mod A number that is used as a factor to convert the unit of the model data, e.g. 86400
#' @param unit.conv.ref A number that is used as a factor to convert the unit of the reference data, e.g. 86400
#' @param variable.unit A string that gives the final units using LaTeX notation, e.g. 'gC m$^{-2}$ day$^{-1}$'
#' @param score.weights R object that gives the weights of each score (\eqn{S_{bias}}, \eqn{S_{rmse}}, \eqn{S_{phase}}, \eqn{S_{iav}}, \eqn{S_{dist}})
#' that are used for computing the overall score, e.g. c(1,2,1,1,1)
#' @param shp.filename A string that gives the coastline shapefile
#' @param my.xlim An R object that gives the longitude range that you wish to
#' plot, e.g. c(-180, 180)
#' @param my.ylim An R object that gives the longitude range that you wish to
#' plot, e.g. c(-90, 90)
#' @param plot.width Number that gives the plot width, e.g. 8
#' @param plot.height Number that gives the plot height, e.g. 4
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory
#'
#' @return (1) Figures in PDF format that show global maps of
#' the model data
#' (mean, \eqn{mod.mean}; interannual-variability, \eqn{mod.iav}),
#' the reference data
#' (mean, \eqn{ref.mean}; interannual-variability, \eqn{ref.iav}),
#' statistical metrics
#' (bias, \eqn{bias}; centralized root mean square error, \eqn{crmse}),
#' and scores
#' (bias score, \eqn{bias.score}; root mean square error score, \eqn{rmse.score}; inter-annual variability score \eqn{iav.score}).
#'
#' (2) Three text files: (i) score values and (ii) score inputs averaged across
#' the entire grid, and (iii) score values for each individual river basin.

#' @examples
#'
#' \donttest{
#' library(amber)
#' library(ncdf4)
#' library(raster)
#' library(latex2exp)
#'
#' long.name <- 'Streamflow'
#' nc.mod <- system.file('extdata/modelRegular', 'mrro_monthly.nc', package = 'amber')
#' nc.ref <- system.file('extdata/referenceRegular', 'runoff.nc', package = 'amber')
#' nc.basins <- system.file('extdata/referenceRegular', 'basins.nc', package = 'amber')
#' mod.id <- 'CLASSIC' # model name
#' ref.id <- 'GRDC' # give reference dataset a name
#' unit.conv.mod <- 86400 # optional unit conversion for model data
#' unit.conv.ref <- 86400 # optional unit conversion for reference data
#' variable.unit <- 'kg m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
#' score.weights <- c(1,2,1,1,1) # define score weights
#'
#' scores.runoff(long.name, nc.mod, nc.ref, nc.basins, mod.id, ref.id, unit.conv.mod,
#' unit.conv.ref, variable.unit, score.weights)
#' }
#' @export
scores.runoff <- function(long.name, nc.mod, nc.ref, nc.basins, mod.id, ref.id, unit.conv.mod, unit.conv.ref, variable.unit,
    score.weights = c(1, 2, 1, 1, 1), shp.filename = system.file("extdata/ne_110m_land/ne_110m_land.shp", package = "amber"),
    my.xlim = c(-180, 180), my.ylim = c(-60, 85), plot.width = 8, plot.height = 3.8, outputDir = FALSE) {

    # Main steps:

    # (1) model data: get variable name, dates, and data

    # (2) reference data: get variable name, dates, data, and river basin names (site-level)

    # (3) get basin polygon, basin_index, and basin area

    # (4) subset common time period

    # (5) convert reference data from monthly to yearly means

    # (6) compute sum of modeled runoff for each basin

    # (7) convert model data from monthly to yearly means

    # (8) statistical analysis

    # comment: each of the 50 basins has a basin_index, which ranges from 0 to 49

    # All data must be sorted such that the basin_index increases from 0 to 49 (basins <-
    # basins[order(basins$basin_index),]).

    # projection
    regular <- "+proj=longlat +ellps=WGS84"

    # (1) model data: get variable name, dates, and data
    nc <- ncdf4::nc_open(nc.mod)
    variable.name <- names(nc[["var"]])
    ncdf4::nc_close(nc)
    variable.name <- variable.name[length(variable.name)]  # take the last variable (relevant for CanESM5)
    variable.name <- ifelse(variable.name == "burntFractionAll", "burnt", variable.name)  # rename burntFractionAll to shorter name
    variable.name <- toupper(variable.name)  # make variable name upper-case
    mod <- raster::brick(nc.mod)  # unit: kg m-2 s-1
    area.gridcell <- raster::area(raster::subset(mod, 1:1)) * 1000 * 1000  # area of grid cell in m^2
    dummy <- raster::subset(mod, 1:1)  # used later for legend
    dates.mod <- raster::getZ(mod)
    dates.mod <- format(as.Date(dates.mod), "%Y-%m")  # only year and month
    start.date.mod <- min(dates.mod)  # format: 'YYYY-MM'
    end.date.mod <- max(dates.mod)  # format: 'YYYY-MM'

    # (2) reference data: get variable name, dates, data, and river basin names

    nc <- ncdf4::nc_open(nc.ref)
    variable.name.ref <- names(nc[["var"]])
    variable.name.ref <- variable.name.ref[4]
    # get reference date
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
    data <- data.frame(data * unit.conv.ref)  # unit: kg m-2 day-1
    colnames(data) <- c(dates.ref)
    # get river basin names
    nc_atts <- ncdf4::ncatt_get(nc, 0)
    basin.names <- nc_atts$site_name
    basin.names <- gsub("\\)", "", as.character(basin.names))  # omit parentheses
    basin.names <- gsub("\\(", "-", as.character(basin.names))  # omit parentheses
    basin.names <- gsub(" -", "-", as.character(basin.names))  # omit white space
    basin.names <- gsub("Rio Grande-Bra", "Rio Grande", as.character(basin.names))  # this is Rio Grande in N. America, not S. America (Brazil)
    basin.names <- gsub("Colorado-AR", "Colorado", as.character(basin.names))  # this is Colorado in N. America, not S. America (Argentina)
    basin.names <- strsplit(basin.names, ",")  # split
    basin.names <- data.frame(basin.names)
    colnames(basin.names) <- "basin.names"
    data <- data.frame(basin.names, data)
    data <- intFun.site.points(lon, lat, data)
    raster::projection(data) <- regular
    ref <- data
    ncdf4::nc_close(nc)

    # (3) get basin polygon, basin_index, and basin area

    basins <- raster::raster(nc.basins)
    basins <- raster::rasterToPolygons(basins, fun = NULL, n = 4, na.rm = TRUE, digits = 12, dissolve = TRUE)
    basins <- basins[order(basins$basin_index), ]  # sort by basin_index
    raster::projection(basins) <- regular
    # overlay polygon and points to get basin_id
    basin_index <- data.frame(sp::over(ref, basins))  # basin index ranges from 0 to 49)
    ref <- as.data.frame(ref)
    ref <- data.frame(basin_index, ref)
    ref <- ref[order(ref$basin_index), ]
    basin.area <- raster::area(basins)  # get area of each basin in m^2
    basin.area <- data.frame(basin_index, basin.area)
    basin.area <- basin.area[order(basin.area$basin_index), ]
    basin.area <- t(data.frame(basin.area))
    colnames(basin.area) <- c()

    # (4) subset common time period

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
    basin_index <- ref$basin_index
    basin.name <- ref$basin.name
    lon <- ref$lon
    lat <- ref$lat
    ref.info <- data.frame(basin_index, basin.name, lon, lat)
    # subset common time period for reference data
    ref <- data.frame(ref[4 + column.id])  # this subsets the common dates
    colnames(ref) <- c(dates)

    # (5) convert reference data from monthly to yearly means

    ref <- t(ref)
    year <- rownames(ref)
    year <- as.Date(year, c("%Y-%M"))
    year <- format(year, "%Y")
    year <- as.numeric(year)
    ref <- data.frame(year, ref)
    rownames(ref) <- c()
    ref <- apply(ref[, -1], 2, function(x) tapply(x, ref$year, mean, na.rm = FALSE))
    ref.info <- t(ref.info)
    basin_index <- ref.info[1, ]
    colnames(ref) <- basin_index

    # (6) compute sum of modeled runoff for each basin

    mod <- mod * unit.conv.mod  # unit: kg m-2 day-1
    mod <- mod * area.gridcell  # unit: kg gridcell-1 day-1
    mod <- raster::rotate(mod)
    basins <- raster::raster(nc.basins)
    basins <- raster::rasterToPolygons(basins, fun = NULL, n = 4, na.rm = TRUE, digits = 12, dissolve = TRUE)
    basins <- basins[order(basins$basin_index), ]
    basin.area <- raster::area(basins)  # get area of each basin in m^2
    basin.area <- t(data.frame(basin.area))
    colnames(basin.area) <- c()
    mod <- raster::extract(mod, basins)
    # this gives a list of 50 elements, where one element is one basin the columns of each element are months the rows of
    # each element give the number of grid cells next, sum up the rows for each column
    mod <- data.frame(lapply(mod, colSums, na.rm = TRUE))  # unit: kg basin-1 day-1
    # ncol(mod); nrow(mod)
    colnames(mod) <- basin_index
    rownames(mod) <- paste("X", common.dates$dates, sep = "")
    colnames(basin.area) <- basin.names$basin.names
    rownames(basin.area) <- c()
    mod <- apply(mod, 1, function(x) x/basin.area)  # divide by river basin area (unit: kg m-2 day-1)
    mod <- t(mod)
    colnames(mod) <- basin_index

    # (7) convert model data from monthly to yearly means

    year <- rownames(mod)
    year <- gsub("X", "", year)
    year <- as.Date(year, c("%Y-%M"))
    year <- format(year, "%Y")
    year <- as.numeric(year)
    mod <- data.frame(year, mod)
    rownames(mod) <- c()
    mod <- apply(mod[, -1], 2, function(x) tapply(x, mod$year, mean, na.rm = TRUE))
    colnames(mod) <- basin_index
    # check whether number of columns and rows match between mod and ref
    stopifnot(ncol(ref) == ncol(mod))
    stopifnot(nrow(ref) == nrow(mod))
    # head(mod) head(ref) make a string that summarizes metadata
    meta.data.mod <- paste(variable.name, mod.id, "from", start.date, "to", end.date, sep = "_")
    meta.data.ref <- paste(variable.name, ref.id, "from", start.date, "to", end.date, sep = "_")
    meta.data.com <- paste(variable.name, mod.id, "vs", ref.id, "from", start.date, "to", end.date, sep = "_")

    # (8) statistical analysis

    # The time period of the reference data varies among sites.

    # Some sites are located in the model ocean.

    # Exclude data in one data set that is not available in the other data set.

    mask.mod <- (mod - mod + 1)
    mask.ref <- (ref - ref + 1)
    mask <- mask.mod * mask.ref
    mod <- mod * mask
    ref <- ref * mask
    mod <- data.frame(mod)
    ref <- data.frame(ref)
    # get number of years of data coverage of each station
    n.years <- apply(mask, 2, sum, na.rm = TRUE)
    n.years <- intFun.site.points(lon, lat, n.years)
    # (1) bias
    mod.mean <- apply(mod, 2, mean, na.rm = TRUE)  # time mean
    ref.mean <- apply(ref, 2, mean, na.rm = TRUE)  # time mean
    weights <- ref.mean  # weights used for spatial integral
    bias <- mod.mean - ref.mean  # time mean
    ref.sd <- apply(ref, 2, sd, na.rm = TRUE)  # standard deviation of reference data
    # epsilon_bias <- abs(bias)/ref.sd  # relative error
    epsilon_bias <- abs(bias)/ref.mean  # relative error
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
    ref.sd.scalar <- mean(ref.sd, na.rm = TRUE)  # global mean value
    epsilon_bias.scalar <- mean(epsilon_bias, na.rm = TRUE)  # global mean value

    # (2) root mean square error (rmse)

    rmse <- mapply(intFun.rmse, mod, ref)
    mod.anom <- data.frame(apply(mod, 2, intFun.anom))  # compute anomalies
    ref.anom <- data.frame(apply(ref, 2, intFun.anom))  # compute anomalies
    crmse <- mapply(intFun.crmse, mod.anom, ref.anom)
    # epsilon_rmse <- crmse/ref.sd  # relative error
    epsilon_rmse <- crmse/ref.mean  # relative error
    rmse.score <- exp(-epsilon_rmse)  # rmse score as a function of space
    S_rmse_not.weighted <- mean(rmse.score, na.rm = TRUE)  # scalar score (not weighted)
    # calculate the weighted scalar score
    a <- rmse.score * weights  # this is a raster
    b <- sum(a, na.rm = TRUE)  # this is a scalar, the sum up all values
    S_rmse_weighted <- b/sum(weights, na.rm = TRUE)  # scalar score (weighted)
    # compute global mean values of score input(s)
    rmse.scalar <- mean(rmse, na.rm = TRUE)  # global mean value
    crmse.scalar <- mean(crmse, na.rm = TRUE)  # global mean value
    epsilon_rmse.scalar <- mean(epsilon_rmse, na.rm = TRUE)  # global mean value

    # (3) phase shift does not apply since the analysis is based on annual data

    mod.max.month <- rep(NA, ncol(mod))
    ref.max.month <- rep(NA, ncol(mod))
    phase <- rep(NA, ncol(mod))
    phase.score <- rep(NA, ncol(mod))
    S_phase_not.weighted <- NA
    S_phase_weighted <- NA
    phase.scalar <- NA

    # (4) interannual variability based on annual data

    mod.iav <- apply(mod.anom, 2, intFun.iav)
    ref.iav <- apply(ref.anom, 2, intFun.iav)
    # set values close to zero to NA
    ref.iav.na <- ref.iav
    ref.iav.na[ref.iav.na < 10^(-5)] <- NA
    #
    epsilon_iav <- abs((mod.iav - ref.iav))/ref.iav.na  # I changed Eq. 26 so that epsilon_iav >= 0
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
    # for plotting purposes: mod.iav.points <- intFun.site.points(lon, lat, mod.iav) ref.iav.points <-
    # intFun.site.points(lon, lat, ref.iav) iav.score.points <- intFun.site.points(lon, lat, iav.score)

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
    S_overall <- (w.bias * S_bias + w.rmse * S_rmse + w.iav * S_iav + w.dist * S_dist)/(w.bias + w.rmse + w.iav + w.dist)
    scores <- data.frame(variable.name, ref.id, S_bias, S_rmse, S_phase, S_iav, S_dist, S_overall)
    scores_not.weighted <- scores
    # weighted (except for S_dist)
    S_bias <- S_bias_weighted
    S_rmse <- S_rmse_weighted
    S_phase <- S_phase_weighted
    S_iav <- S_iav_weighted
    S_overall <- (w.bias * S_bias + w.rmse * S_rmse + w.iav * S_iav + w.dist * S_dist)/(w.bias + w.rmse + w.iav + w.dist)
    scores <- data.frame(variable.name, ref.id, S_bias, S_rmse, S_phase, S_iav, S_dist, S_overall)
    scores_weighted <- scores
    #
    scores <- rbind(scores_not.weighted, scores_weighted)
    rownames(scores) <- c("not.weighted", "weighted")
    if (outputDir != FALSE) {
        utils::write.table(scores, paste(outputDir, "/", "scorevalues", "_", meta.data.com, sep = ""))
    }

    # Get all score values in case you want to compare this run against another run.

    # Having all values will enable you to conduct a significance test.

    dist.score <- rep(S_dist, length(bias.score))
    all.score.values <- data.frame(bias.score, rmse.score, phase.score, iav.score, dist.score)
    colnames(all.score.values) <- c("bias.score", "rmse.score", "phase.score", "iav.score", "dist.score")
    all.score.values[is.na(all.score.values)] <- NA  # converts all NaN to NA
    if (outputDir != FALSE) {
        utils::write.table(all.score.values, paste(outputDir, "/", "allscorevalues", "-", variable.name, "-", ref.id, sep = ""))
    }
    # Selected score inputs
    mod.max.month.scalar <- NA
    ref.max.month.scalar <- NA
    phase.scalar <- NA
    S_phase_not.weighted <- NA

    scoreinputs <- data.frame(long.name, variable.name, ref.id, variable.unit, mod.mean.scalar, ref.mean.scalar, bias.scalar,
        ref.sd.scalar, epsilon_bias.scalar, S_bias_not.weighted, rmse.scalar, crmse.scalar, ref.sd.scalar, epsilon_rmse.scalar,
        S_rmse_not.weighted, mod.max.month.scalar, ref.max.month.scalar, phase.scalar, S_phase_not.weighted, mod.iav.scalar,
        ref.iav.scalar, epsilon_iav.scalar, S_iav_not.weighted, mod.sigma.scalar, ref.sigma.scalar, sigma, R, S_dist)

    if (outputDir != FALSE) {
        utils::write.table(scoreinputs, paste(outputDir, "/", "scoreinputs", "_", meta.data.com, sep = ""))
    }

    # Apply the functions min.max.int and min.max.int.mod.ref min.max.int

    mmi.bias <- intFun.min.max.int.bias(bias)
    mmi.bias.score <- c(0, 1, 0.1)
    mmi.crmse <- intFun.min.max.int(crmse)
    mmi.rmse.score <- c(0, 1, 0.1)
    # mmi.phase <- c(0,6,1) mmi.phase.score <- intFun.min.max.int(phase.score)
    mmi.iav.score <- c(0, 1, 0.1)
    # min.max.int.mod.ref
    mmi.mean <- intFun.min.max.int.mod.ref(mod.mean, ref.mean)
    # mmi.max.month <- c(1,12,1)
    mmi.iav <- intFun.min.max.int.mod.ref(mod.iav, ref.iav)
    # Add metadata:

    # 1. figure title (e.g. Mean_nee_ModID_123_from_1982-01_to_2008-12)

    # 2. min, max, interval used in legend (e.g. 0, 1, 0.1)

    # 3. legend bar text (e.g. 'score (-)')

    attr(mod.mean, "metadata") <- list(paste("Mean", meta.data.mod, sep = "_"), mmi.mean, variable.unit)
    attr(ref.mean, "metadata") <- list(paste("Mean", meta.data.ref, sep = "_"), mmi.mean, variable.unit)
    attr(bias, "metadata") <- list(paste("Bias", meta.data.com, sep = "_"), mmi.bias, variable.unit)
    attr(bias.score, "metadata") <- list(paste("Bias_score", meta.data.com, sep = "_"), mmi.bias.score, "score (-)")
    attr(crmse, "metadata") <- list(paste("CRMSE", meta.data.com, sep = "_"), mmi.crmse, variable.unit)
    attr(rmse.score, "metadata") <- list(paste("RMSE_score", meta.data.com, sep = "_"), mmi.rmse.score, "score (-)")
    # attr(mod.max.month,'metadata') <- list(paste('Month_with_max', meta.data.mod, sep='_'), mmi.max.month, 'month')
    # attr(ref.max.month,'metadata') <- list(paste('Month_with_max', meta.data.ref, sep='_'), mmi.max.month, 'month')
    # attr(phase,'metadata') <- list(paste('Diff_in_max_month', meta.data.com, sep='_'), mmi.phase, 'month')
    # attr(phase.score,'metadata') <- list(paste('Seasonality_score', meta.data.com, sep='_'), mmi.phase.score, 'score (-)')
    attr(mod.iav, "metadata") <- list(paste("IAV", meta.data.mod, sep = "_"), mmi.iav, variable.unit)
    attr(ref.iav, "metadata") <- list(paste("IAV", meta.data.ref, sep = "_"), mmi.iav, variable.unit)
    attr(iav.score, "metadata") <- list(paste("IAV_score", meta.data.com, sep = "_"), mmi.iav.score, "score (-)")

    # write data to file

    stat.metric <- data.frame(basin_index, basin.names, lon, lat, mod.mean, ref.mean, bias, bias.score, crmse, rmse.score,
        mod.max.month, ref.max.month, phase, phase.score, mod.iav, ref.iav, iav.score)
    colnames(stat.metric) <- c("index", "basin", "lon", "lat", "mod.mean", "ref.mean", "bias", "bias.score", "crmse", "rmse.score",
        "mod.max.month", "ref.max.month", "phase", "phase.score", "mod.iav", "ref.iav", "iav.score")

    my.filename <- paste(variable.name, "runoff", sep = "_")
    if (outputDir != FALSE) {
        utils::write.table(stat.metric, paste(outputDir, my.filename, sep = "/"))
    }

    # plot data read in land polygon used in plots
    land <- raster::shapefile(shp.filename)
    regular <- "+proj=longlat +ellps=WGS84"
    raster::projection(land) <- regular

    # Drop columns related to phase, as phase is not evaluated for annual runoff

    stat.metric <- subset(stat.metric, select = -c(mod.max.month, ref.max.month, phase, phase.score))
    stat.metric <- stat.metric[order(stat.metric$index), ]

    for (i in 5:ncol(stat.metric)) {
        data <- stat.metric[i]
        basins <- basins[order(basins$basin_index), ]  # sort basins by basin_index
        values <- data
        colnames(values) <- "values"
        data <- data.frame(stat.metric$lon, stat.metric$lat, stat.metric$index, values$values)
        colnames(data) <- c("lon", "lat", "basin_index", "values")
        data <- data[order(data$basin_index), ]  # sort data by basin_index
        data <- as.data.frame(data)
        data.index <- data  # used later for plotting basin_index
        data <- intFun.site.points(data$lon, data$lat, data$values)  # make spatial points
        data <- as.data.frame(data)
        cname <- names(stat.metric[i])
        x <- get(cname)  # get the data
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

        if (i == 7)
            {
                my.col <- my.col.bias
            }  # divergent color scheme for bias plots
        color.id <- seq(1, length(my.col), 1)
        color.id <- data.frame(color.id, my.col)
        my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
        my.legend.args <- list(text = legend.bar.text, side = 2, font = 1, line = 1, cex = 1)
        fun.colors <- function(x, breaks) {
            viridis::viridis(length(breaks) - 1, direction = -1)[cut(x, breaks)]
        }
        if (i == 7) {
            fun.colors <- function(x, breaks) {
                scico::scico(n = length(my.breaks) - 1, palette = "vik")[cut(x, breaks)]
            }
        }
        colors <- fun.colors(data$data, my.breaks)
        data <- intFun.site.points(data$lon, data$lat, data)
        # create a dummy raster file used for plotting the legend
        dummy <- raster::rotate(dummy)
        my.extent <- raster::extent(dummy)
        random.values <- stats::runif(ncol(dummy) * nrow(dummy), min, max)
        dummy <- matrix(random.values, ncol = ncol(dummy))
        dummy <- raster::raster(dummy)
        raster::extent(dummy) <- my.extent
        # plot
        my.plotname <- paste(my.filename, cname, sep = "-")
        my.plotname <- gsub("_", "-", my.plotname)
        my.plotname <- gsub(".", "-", my.plotname, fixed = TRUE)

        # plot
        oldpar <- graphics::par(mfrow = c(1, 2))
        on.exit(graphics::par(oldpar))
        if (outputDir != FALSE) {
            grDevices::pdf(paste(outputDir, "/", my.plotname, ".pdf", sep = ""), width = plot.width, height = plot.height)
        }
        graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(3, 3, 3, 4), lwd = 1, cex = 1)
        plot(dummy, col = NA, legend = FALSE, ylim = c(-50, 85), main = paste(long.name, my.title, sep = "\n"), axes = FALSE)
        plot(land, col = "grey", border = NA, add = TRUE)
        plot(basins, border = "white", lwd = 0.1, col = colors, add = TRUE)
        values <- data.frame(data)
        graphics::points(data, pch = 16, cex = 0.5, col = "red")
        # add values to plot
        my.pointLabel <- toString(round(data$data, 1))
        my.pointLabel <- unlist(strsplit(my.pointLabel, ", "))
        # pointLabel(data$lon, data$lat, labels = my.pointLabel, cex=0.8, col='black')
        axis(1, labels = TRUE, tcl = 0.3)
        axis(2, labels = TRUE, tcl = 0.3, las = 2)
        graphics::box()
        plot(dummy, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args,
            legend.width = 1.5, legend.shrink = 1, font = 1)
        if (outputDir != FALSE) {
            dev.off()
        }
    }

    # Additional plots scatter plot of annual mean streamflow
    x <- ref.mean
    y <- mod.mean
    my.lm <- stats::lm(y ~ x)
    R2 <- summary(my.lm)$r.squared
    R2 <- round(R2, 2)

    my.xlab <- expression(paste("Reference data ", "(kg m"^{
        -2
    } ~ "day"^{
        -1
    } ~ ")", sep = ""))
    my.ylab <- expression(paste("Model data ", "(kg m"^{
        -2
    } ~ "day"^{
        -1
    } ~ ")", sep = ""))
    my.legend <- substitute(paste("R"^{
        2
    } ~ "= ", R2.value, ", ", italic("n"), " = 50"), list(R2.value = R2))

    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", "MRRO-scatterplot.pdf", sep = ""), width = 4, height = 4)
    }
    graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(5, 5, 2, 2), lwd = 1, cex.main = 1, tck = 0.03)
    plot(x, y, pch = 16, xlim = c(min = mmi.mean[1], max = mmi.mean[2]), ylim = c(min = mmi.mean[1], max = mmi.mean[2]),
        xlab = my.xlab, ylab = my.ylab, main = "Annual mean streamflow of 50 river basins", las = 1)
    graphics::legend("topleft", legend = my.legend, bty = "n")
    graphics::abline(0, 1)
    graphics::box()
    # Ticks
    axis(3, labels = FALSE, tcl = 0.3)
    axis(4, labels = FALSE, tcl = 0.3)
    if (outputDir != FALSE) {
        dev.off()
    }
}
