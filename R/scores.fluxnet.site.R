################################################################################
#' Scores for FLUXNET reference data when model run at FLUXNET site
#' @description This function compares model output in CSV format against
#' FLUXNET measurements in CSV format. Use this function when running your model
#' at each FLUXNET site individually. The performance of a model is
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
#' @param mod.csv A string that gives the name of the model output in csv format, e.g. 'gpp_monthly.csv'
#' @param mod.csv.path A string that gives the path to the model output for site-level runs
#' @param ref.csv A string that gives the path and name of the csv file that contains the reference data output, e.g. '/home/reference_gpp.csv'. The columns of this file should contain latitude, longitude, date, variable of interest, and site name.
#' @param mod.id A string that identifies the source of the reference data set, e.g. 'CanESM2'
#' @param ref.id A string that identifies the source of the reference data set, e.g. 'MODIS'
#' @param unit.conv.mod A number that is used as a factor to convert the unit of the model data, e.g. 86400
#' @param unit.conv.ref A number that is used as a factor to convert the unit of the reference data, e.g. 86400
#' @param variable.unit A string that gives the final units using LaTeX notation, e.g. 'gC m$^{-2}$ day$^{-1}$'
#' @param sites A vector of strings that give the fluxnet site names, e.g.  c('AU-Tum','BR-Sa1','CA-Qfo')
#' @param score.weights R object that gives the weights of each score (\eqn{S_{bias}}, \eqn{S_{rmse}}, \eqn{S_{phase}}, \eqn{S_{iav}}, \eqn{S_{dist}})
#' that are used for computing the overall score, e.g. c(1,2,1,1,1)
#' @param my.xlim An R object that gives the longitude range that you wish to
#' plot, e.g. c(-180, 180)
#' @param my.ylim An R object that gives the longitude range that you wish to
#' plot, e.g. c(-90, 90)
#' @param plot.width Number that gives the plot width, e.g. 8
#' @param plot.height Number that gives the plot height, e.g. 4
#' @param numCores An integer that defines the number of cores, e.g. 2
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @param phaseMinMax A string (either 'phaseMax' or 'phaseMin') that determines
#' whether to assess the seasonal peak as a maximum or a minimum. The latter may be appropriate for variables
#' that tend to be negative, such as net longwave radiation or net ecosystem exchange.
#' @param myCex A number that determines the size of the dots in the Figure. Default is set to 0.7.
#' @return (1) Figures in PDF format that show maps of
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
#' long.name <- 'Gross primary productivity'
#' mod.csv <- 'gpp_monthly.csv'
#' mod.csv.path <- system.file('extdata/siteLevelRun', package = 'amber')
#' ref.csv <- system.file('extdata/referenceRegular', 'gpp_monthly_fluxnet.csv', package = 'amber')
#' mod.id <- 'CLASSIC-Sitelevel' # define a model experiment ID
#' ref.id <- 'FLUXNET' # give reference dataset a name
#' unit.conv.mod <- 86400*1000 # optional unit conversion for model data
#' unit.conv.ref <- 1 # optional unit conversion for reference data
#' variable.unit <- 'gC m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)
#'
#' sites <- c('AU-Tum','CA-TPD', 'US-WCr')
#'
#' # Short version using default settings:
#' scores.fluxnet.site(long.name, mod.csv, mod.csv.path, ref.csv, mod.id, ref.id,
#' unit.conv.mod, unit.conv.ref, variable.unit, sites)
#'
#' # Additional parameters:
#' score.weights <- c(1,2,1,1,1) # score weights of S_bias, S_rmse, S_phase, S_iav, S_dist
#' my.xlim <- c(-180, 180)
#' my.ylim <- c(-60, 85)
#' plot.width <- 8
#' plot.height <- 3.8
#' numCores <- 2
#'
#' scores.fluxnet.site(long.name, mod.csv, mod.csv.path, ref.csv, mod.id, ref.id,
#' unit.conv.mod, unit.conv.ref, variable.unit, sites, score.weights,
#' my.xlim, my.ylim, plot.height, numCores)
#'
#' # To zoom into a particular region:
#' scores.fluxnet.site(long.name, mod.csv, mod.csv.path, ref.csv, mod.id, ref.id,
#' unit.conv.mod, unit.conv.ref, variable.unit, sites,
#' my.xlim = c(-150, -60), my.ylim = c(20, 60), plot.width = 6, plot.height = 3.8)
#' } #donttest
#' @export
scores.fluxnet.site <- function(long.name, mod.csv, mod.csv.path, ref.csv, mod.id, ref.id, unit.conv.mod, unit.conv.ref, variable.unit,
    sites, score.weights = c(1, 2, 1, 1, 1), my.xlim = c(-180, 180), my.ylim = c(-60, 85), plot.width = 8, plot.height = 3.8,
    numCores = 2, outputDir = FALSE, phaseMinMax = "phaseMax", myCex = 0.7) {

    #---------------------------------------------------------------------------

    # (I) Data preparation

    #---------------------------------------------------------------------------

    # Comments:

    # Model data consist of one CSV file per site.

    # Reference data consist of a single CSV file.

    # The raw data have the row names lat, lon, time, data, siteID, ...

    # The first section of this script produces the data frames 'mod' and 'ref'.

    # The corresponding column names are Fluxnet site names.

    # The corresponding row names are dates.

    #---------------------------------------------------------------------------

    # Model data: merge individual model output files

    cl <- parallel::makePSOCKcluster(numCores)
    doParallel::registerDoParallel(cl)
    nTime <- length(sites)
    allSites <- foreach::foreach(i = 1:nTime) %dopar% {
        sitename <- sites[i]
        data.path <- file.path(mod.csv.path, sitename, "csv", mod.csv)
        data <- utils::read.csv(data.path)
        sitename <- rep(sitename, nrow(data))
        data <- data.frame(data, sitename)
    }
    mod <- base::do.call(rbind, allSites)
    parallel::stopCluster(cl)
    mod.sites <- data.frame(unique(mod$sitename))
    colnames(mod.sites) <- "sitename"

    variable.name <- colnames(mod[4])
    variable.name <- ifelse(variable.name == "burntFractionAll", "burnt", variable.name)  # rename burntFractionAll to shorter name
    variable.name <- toupper(variable.name)  # make variable name upper-case

    #---------------------------------------------------------------------------

    # Reference data: read in file

    ref <- utils::read.csv(ref.csv)
    ref.sites <- data.frame(unique(ref$sitename))
    colnames(ref.sites) <- "sitename"

    #---------------------------------------------------------------------------

    # Merge mod and ref data:

    # get common sites

    commonSites <- merge(mod.sites, ref.sites, by = "sitename")

    # subset common sites
    mod <- merge(mod, commonSites, by = "sitename")
    ref <- merge(ref, commonSites, by = "sitename")

    # merge by sites and dates
    modRef <- merge(mod, ref, by = c("sitename", "time"))

    # mod <- data.frame(modRef$sitename, modRef$lat.x, modRef$lon.x, modRef$time, modRef[, 5])
    mod <- data.frame(modRef$sitename, modRef$lat.y, modRef$lon.y, modRef$time, modRef[, 5])
    colnames(mod) <- c("sitename", "lat", "lon", "time", "variable")

    ref <- data.frame(modRef$sitename, modRef$lat.y, modRef$lon.y, modRef$time, modRef[, 8])
    colnames(ref) <- c("sitename", "lat", "lon", "time", "variable")

    modRef <- list(mod, ref)

    # At this stage, mod and ref columns are: sitename, lat, lon, time, variable

    #---------------------------------------------------------------------------

    for (x in 1:2) {
        data <- modRef[[x]]
        site.id <- unique(data$sitename)
        nSites <- length(site.id)
        allSites <- split(data, data$sitename)

        # Restructure data
        cl <- parallel::makePSOCKcluster(numCores)
        doParallel::registerDoParallel(cl)
        eachSite <- foreach::foreach(i = 1:nSites) %dopar% {
            singleSite <- as.data.frame(allSites[i])  # single site
            date <- t(singleSite[4])
            date <- as.Date(date)
            date <- format(as.Date(date), "%Y-%m")  # only year and month
            values <- t(singleSite[5])
            values[values == -9999] <- NA
            ID <- as.character(singleSite[1, 1])
            lat <- singleSite[1, 2]
            lon <- singleSite[1, 3]
            singleSite <- data.frame(ID, lon, lat, values)
            colnames(singleSite) <- c("siteID", "lon", "lat", date)
            assign(paste("site", i, sep = "_"), singleSite)
        }
        parallel::stopCluster(cl)

        # Make the merge function to use all.x = TRUE and all.y = TRUE
        fun.merge <- function(x, y) {
            myMerge <- merge(x, y, all.x = TRUE, all.y = TRUE)
            return(myMerge)
        }

        # merge all individual sites and sort row by date
        data <- Reduce(fun.merge, eachSite)
        siteID <- as.character(data$siteID)
        lon <- data$lon
        lat <- data$lat
        data <- data[4:ncol(data)]  # data only
        data <- data[, order(names(data))]
        data <- data.frame(siteID, lon, lat, data * c(unit.conv.mod, unit.conv.ref)[x])
        # make siteID to rownames
        myRownames <- data[, 1]
        data <- data[, -1]  # drop siteID column
        rownames(data) <- myRownames
        assign(paste(c("mod", "ref")[x], sep = ""), data)
    }

    # At this stage, mod and ref columns are lon, lat, dates (X1996.01, X1996.02, etc.).

    #---------------------------------------------------------------------------

    # Check whether mod and ref have the same number of rows and columns.

    stopifnot(ncol(ref) == ncol(mod))
    stopifnot(nrow(ref) == nrow(mod))

    #---------------------------------------------------------------------------

    # Transpose data, where rows are dates and columns are sites

    # this is necessary for using mapply below

    lon <- ref[1]
    lat <- ref[2]

    mod <- data.frame(t(mod))
    ref <- data.frame(t(ref))

    dates <- data.frame(rownames(mod))

    mod <- mod[-(1:2), ]  # drop lon and lat
    ref <- ref[-(1:2), ]  # drop lon and lat
    dates <- dates[-(1:2), ]  # drop lon and lat

    stopifnot(colnames(ref) == colnames(mod))
    stopifnot(rownames(ref) == rownames(mod))

    # At this stage, mod and ref comlumn and row names are Fluxnet site names and dates, respectively.

    #---------------------------------------------------------------------------

    # Exclude data that are not covered by both data sets to make mod and ref comparable.

    mask.mod <- (mod - mod + 1)
    mask.ref <- (ref - ref + 1)
    mask <- mask.mod * mask.ref
    mod <- mod * mask
    ref <- ref * mask

    mod <- data.frame(mod)  # necessary in case I am looking at a single site only
    ref <- data.frame(ref)  # necessary in case I am looking at a single site only
    #---------------------------------------------------------------------------

    # Get all dates, months, starting date, and ending date.

    # dates <- rownames(mod)
    dates <- gsub("X", "", dates)
    dates <- gsub("\\.", "-", dates)
    dates <- paste(dates, "15", sep = "-")
    dates <- as.Date(dates)

    month <- format(as.Date(dates), format = "%m")
    month <- as.numeric(month)

    start.date <- min(dates)
    end.date <- max(dates)
    start.date <- format(as.Date(start.date), "%Y-%m")
    end.date <- format(as.Date(end.date), "%Y-%m")

    #---------------------------------------------------------------------------

    # make a string that summarizes metadata

    meta.data.mod <- paste(variable.name, mod.id, "from", start.date, "to", end.date, sep = "_")
    meta.data.ref <- paste(variable.name, ref.id, "from", start.date, "to", end.date, sep = "_")
    meta.data.com <- paste(variable.name, mod.id, "vs", ref.id, "from", start.date, "to", end.date, sep = "_")

    #---------------------------------------------------------------------------

    # (II) Statistical analysis

    #---------------------------------------------------------------------------

    # (1) bias

    #---------------------------------------------------------------------------

    mod.mean <- apply(mod, 2, mean, na.rm = TRUE)  # time mean
    ref.mean <- apply(ref, 2, mean, na.rm = TRUE)  # time mean
    weights <- ref.mean  # weights used for spatial integral
    bias <- mod.mean - ref.mean  # time mean
    ref.sd <- apply(ref, 2, sd, na.rm = TRUE)  # standard deviation of reference data
    epsilon_bias <- abs(bias)/ref.sd
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
    ref.sd.scalar <- mean(ref.sd, na.rm = TRUE)
    epsilon_bias.scalar <- stats::median(epsilon_bias, na.rm = TRUE)

    #---------------------------------------------------------------------------

    # (2) root mean square error (rmse)

    #---------------------------------------------------------------------------

    rmse <- mapply(intFun.rmse, mod, ref)  # compute rmse
    mod.anom <- data.frame(apply(mod, 2, intFun.anom))  # compute anomalies
    ref.anom <- data.frame(apply(ref, 2, intFun.anom))  # compute anomalies
    crmse <- mapply(intFun.crmse, mod.anom, ref.anom)
    epsilon_rmse <- crmse/ref.sd
    epsilon_rmse[epsilon_rmse == Inf] <- NA  # relative error
    rmse.score <- exp(-epsilon_rmse)  # rmse score as a function of space
    S_rmse_not.weighted <- mean(rmse.score, na.rm = TRUE)  # scalar score (not weighted)
    # calculate the weighted scalar score
    a <- rmse.score * weights  # this is a raster
    b <- sum(a, na.rm = TRUE)  # this is a scalar, the sum up all values
    S_rmse_weighted <- b/sum(weights, na.rm = TRUE)  # scalar score (weighted)
    # compute global mean values of score input(s)
    rmse.scalar <- mean(rmse, na.rm = TRUE)  # global mean value
    crmse.scalar <- mean(crmse, na.rm = TRUE)  # global mean value
    epsilon_rmse.scalar <- stats::median(epsilon_rmse, na.rm = TRUE)

    #---------------------------------------------------------------------------

    # (3) phase shift

    #---------------------------------------------------------------------------

    mod <- data.frame(month, mod)  # add month
    ref <- data.frame(month, ref)  # add month
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

    # find month of seasonal peak

    # In most cases, we are interested in the timing of the seasonal maximum value

    # In some cases, however, the seasonal peak is a minimum, e.g. NEE = RECO - GPP

    if (phaseMinMax == "phaseMax") {
        mod.max.month <- apply(mod.clim.mly, 2, which.max)
        ref.max.month <- apply(ref.clim.mly, 2, which.max)
    }

    if (phaseMinMax == "phaseMin") {
        mod.max.month <- apply(mod.clim.mly, 2, which.min)
        ref.max.month <- apply(ref.clim.mly, 2, which.min)
    }

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

    #---------------------------------------------------------------------------

    # (4) interannual variability

    #---------------------------------------------------------------------------

    # This approach assumes that all data start in Jan.

    # All months after the last Dec will be dropped if data does not end in Dec.


    mod.anom <- intFun.anom.mly(mod, mod.clim.mly.mm)
    ref.anom <- intFun.anom.mly(ref, ref.clim.mly.mm)
    mod.iav <- apply(mod.anom, 2, intFun.iav)
    ref.iav <- apply(ref.anom, 2, intFun.iav)

    # set values close to zero to NA
    ref.iav.na <- ref.iav
    ref.iav.na[ref.iav.na < 10^(-5)] <- NA

    epsilon_iav <- abs((mod.iav - ref.iav))/ref.iav.na
    epsilon_iav[epsilon_iav == Inf] <- NA  # I changed Eq. 26 so that epsilon_iav > =  0
    iav.score <- exp(-epsilon_iav)  # iav score as a function of space
    S_iav_not.weighted <- mean(iav.score, na.rm = TRUE)  # scalar score (not weighted)
    # calculate the weighted scalar score
    a <- iav.score * weights  # this is a raster
    b <- sum(a, na.rm = TRUE)  # this is a scalar, the sum up all values
    S_iav_weighted <- b/sum(weights, na.rm = TRUE)  # scalar score (weighted)
    # compute global mean values of score input(s)
    mod.iav.scalar <- mean(mod.iav, na.rm = TRUE)  # global mean value
    ref.iav.scalar <- mean(ref.iav, na.rm = TRUE)  # global mean value
    epsilon_iav.scalar <- stats::median(epsilon_iav, na.rm = TRUE)  # global mean value

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

    scores <- rbind(scores_not.weighted, scores_weighted)
    rownames(scores) <- c("not.weighted", "weighted")
    if (outputDir != FALSE) {
        utils::write.table(scores, paste(outputDir, "/", "scorevalues", "_", meta.data.com, sep = ""))
    }

    #---------------------------------------------------------------------------

    # Get all score values to compare two runs using a significance test

    #---------------------------------------------------------------------------

    dist.score <- rep(S_dist, length(bias.score))
    all.score.values <- data.frame(bias.score, rmse.score, phase.score, iav.score, dist.score)
    colnames(all.score.values) <- c("bias.score", "rmse.score", "phase.score", "iav.score", "dist.score")
    all.score.values[is.na(all.score.values)] <- NA  # converts all NaN to NA
    if (outputDir != FALSE) {
        utils::write.table(all.score.values, paste(outputDir, "/", "allscorevalues", "-", variable.name, "-", ref.id, sep = ""))
    }

    #---------------------------------------------------------------------------

    # selected score inputs
    scoreinputs <- data.frame(long.name, variable.name, ref.id, variable.unit, mod.mean.scalar, ref.mean.scalar, bias.scalar,
        bias.scalar.rel, ref.sd.scalar, epsilon_bias.scalar, S_bias_not.weighted, rmse.scalar, crmse.scalar, ref.sd.scalar, epsilon_rmse.scalar,
        S_rmse_not.weighted, mod.max.month.scalar, ref.max.month.scalar, phase.scalar, S_phase_not.weighted, mod.iav.scalar,
        ref.iav.scalar, epsilon_iav.scalar, S_iav_not.weighted, mod.sigma.scalar, ref.sigma.scalar, sigma, R, S_dist)
    if (outputDir != FALSE) {
        utils::write.table(scoreinputs, paste(outputDir, "/", "scoreinputs", "_", meta.data.com, sep = ""))
    }

    #---------------------------------------------------------------------------
    # Only make plots if I assess > 1 site
    if (length(bias) > 1) {
        # For legend: minimum, maximum, interval.

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

        #---------------------------------------------------------------------------

        # Metadata:

        # 1. figure title (e.g. Mean_nee_ModID_123_from_1982-01_to_2008-12)

        # 2. min, max, interval used in legend (e.g. 0, 1, 0.1)

        # 3. legend bar text (e.g. 'score (-)')

        #---------------------------------------------------------------------------

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

        #---------------------------------------------------------------------------

        # coastline

        shp.filename = system.file("extdata/ne_110m_land/ne_110m_land.shp", package = "amber")
        land <- raster::shapefile(shp.filename)

        #---------------------------------------------------------------------------

        # Plot data

        #---------------------------------------------------------------------------

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
            # create a dummy raster layer mod.mean
            dummy <- stats::runif(360 * 180, min = min, max = max)
            dummy <- matrix(dummy, nrow = 180)
            dummy <- raster::raster(dummy)
            raster::extent(dummy) <- c(-180, 180, -90, 90)
            raster::plot(dummy, col = NA, xlim = my.xlim, ylim = my.ylim, main = paste(long.name, my.title, sep = "\n"), axes = FALSE,
                legend = FALSE)
            raster::plot(land, col = "grey", border = NA, add = TRUE)
            graphics::points(data, pch = 16, col = colors, cex = myCex)
            graphics::axis(1, labels = TRUE, tcl = 0.3)
            graphics::axis(2, labels = TRUE, tcl = 0.3, las = 2)
            plot(dummy, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args,
                legend.width = 1.5, legend.shrink = 1, font = 1)
            if (outputDir != FALSE) {
                grDevices::dev.off()
            }
        }
    }
}

