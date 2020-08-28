#-------------------------------------------------------------------------------
# This file contains all internal functions, which consist of a few lines only.  They are used by other functions in the
# package that do not start with 'intFun.'
#-------------------------------------------------------------------------------
# intFun.grid.define.outlier
#' Upper and lower threshold values that define outliers
#' @description Plotting raster objects that contain extreme outliers lead to
#' plots where most grid cells are presented by a single colour since the color
#' legend covers the entire range of values. To avoid this, the user may define
#' upper and lower threshold values, which can then be set to NA and marked as
#' outliers. The threshold values are defined as the interquartile range
#' multiplied by a user-specified factor.
#' @param x A raster object
#' @param y A number
#' @return An upper and a lower threshold value that enclose a range beyond
#' which values are considered as outliers.
#' @examples

#' library(raster)
#' # make some data
#' data <- runif(100,-10,10)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' # compute upper and lower threshold values for outliers
#' threshold.values <- intFun.grid.define.outlier(data, 2)

#' @keywords internal
#' @export
intFun.grid.define.outlier <- function(x, y) {
    my.values <- raster::getValues(x)
    q75 <- stats::quantile(my.values, probs = 0.75, names = FALSE, na.rm = TRUE)
    q25 <- stats::quantile(my.values, probs = 0.25, names = FALSE, na.rm = TRUE)
    iqr <- q75 - q25
    outlier.pos <- q75 + abs(y * iqr)
    outlier.neg <- q25 - abs(y * iqr)
    y <- c(outlier.neg, outlier.pos)
    return(y)
}

#-------------------------------------------------------------------------------

# intFun.grid.outliers.na
#' Set all values outside a range to NA
#' @description Plotting raster objects that contain extreme outliers lead to
#' plots where most grid cells are presented by a single color since the color
#' legend covers the entire range of values. To avoid this, the user may define
#' upper and lower threshold values, which can then be set to NA and marked as
#' outliers.
#' @param x A raster object
#' @param outlier.neg A number that presents the lower bound of the interquartile range
#' @param outlier.pos A number that presents the upper bound of the interquartile range
#' @return A raster object where all outliers are set to NA
#' @examples
#'
#' library(raster)
#' # create a raster object with outliers
#' data <- c(runif(98,-10,10), -1000, 1000)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' # compute upper and lower threshold values for outliers
#' threshold.values <- intFun.grid.define.outlier(data, 2)
#' # Set all values outside the outlier range to NA
#' outlier.neg <- threshold.values[1]
#' outlier.pos <- threshold.values[2]
#' data.no.outliers <- intFun.grid.outliers.na(data, outlier.neg, outlier.pos)
#' plot(data.no.outliers)
#'
#' @keywords internal
#' @export
intFun.grid.outliers.na <- function(x, outlier.neg, outlier.pos) {
    x[x < outlier.neg] <- NA
    x[x > outlier.pos] <- NA
    return(x)
}

#-------------------------------------------------------------------------------

# intFun.grid.outliers.points
#' Convert outliers to spatial points
#' @description Plotting raster objects that contain extreme outliers lead to
#' plots where most grid cells are presented by a single color since the color
#' legend covers the entire range of values. To avoid this, the user may define
#' upper and lower threshold values, which can then be set to NA and marked as
#' outliers. This function converts outliers to spatial points.
#' @param x A raster object
#' @param outlier.neg A number
#' @param outlier.pos A number
#' @return A raster object where all outliers are set to NA
#' @examples
#'
#' library(raster)
#' # create a raster object with outliers
#' data <- c(runif(98,-10,10), -1000, 1000)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' # compute upper and lower threshold values for outliers
#' threshold.values <- intFun.grid.define.outlier(data, 2)
#' # Set all values outside the outlier range to NA
#' outlier.points <- intFun.grid.outliers.points(data, threshold.values[1], threshold.values[2])
#' plot(data); points(outlier.points)
#'
#' @keywords internal
#' @export
intFun.grid.outliers.points <- function(x, outlier.neg, outlier.pos) {
    x[x > outlier.neg & x < outlier.pos] <- NA
    outlier.points <- raster::rasterToPoints(x)
    return(outlier.points)
}

#-------------------------------------------------------------------------------

# intFun.grid.na
#' Replace NA with zero in raster object
#' @description This function sets all NA's in a raster object to zero.
#' @param x A raster object
#' @return A raster object where all NA's have been replaced with zero
#' @examples
#'
#' library(raster)
#' # create a raster object with NA
#' data <- c(runif(99,-10,10),NA)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' # replace NA with zero in raster object
#' data <- intFun.grid.na(data)
#' plot(data); text(data)
#'
#' @keywords internal
#' @export
intFun.grid.na <- function(x) {
    x[is.na(x)] = 0
    return(x)
}  # relaces zero with NA

#-------------------------------------------------------------------------------

# intFun.grid.phase
#' Time distance in months
#' @description Consider two rasters that show the month of the maximum value of
#' a variable during the climatological mean annual cycle. Calculating the
#' absolute difference between both rasters may yield values that range from
#' 0 to 11 (type a). A more meaningful range that considers the circular nature
#' of the annual cycle would range from 0 to 6 (type b). This function converts
#' a raster of type (a) to type (b).
#' @param x A raster object of type (a)
#' @return A raster object of type (b)
#' @examples
#'
#' library(raster)
#' # make some data
#' data <- runif(100,0,11)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' # convert raster of type (a) to type (b)
#' data <- calc(data,intFun.grid.phase)
#' plot(data)
#' plot(data); text(data)
#'
#' @keywords internal
#' @export
intFun.grid.phase <- function(x) {
    ifelse(x <= 6, x, (12 - x))
}

#-------------------------------------------------------------------------------

# intFun.grid.significance
#' Set significant differences to NA
#' @description This function sets all values that are equal or less than 0.05
#' to NA. Applied to a raster that shows the p-values of a significance
#' test, the function sets all statistically significant differences to NA
#' at the 5 percent level.
#' @param x A raster object
#' @return A raster object
#' @examples
#'
#' library(raster)
#' # make some data
#' data <- runif(100,0,1)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' # repeat raster stack 10 times
#' data <- intFun.grid.significance(data)
#' plot(data); text(data, digits=2)
#'
#' @keywords internal
#' @export
intFun.grid.significance <- function(x) {
    x[x <= 0.05] <- NA
    return(x)
}

#-------------------------------------------------------------------------------

# intFun.grid.wilcox
#' Wilcox significance test
#' @description This function conducts the Wilcox two-sided significance test for
#' two data sets that are stored in a single raster stack. The resulting raster
#' object shows the corresponding p-value.
#' @param x A raster stack object
#' @return A raster object
#' @examples
#'
#' library(raster)
#' # create a raster stack
#' for(i in 1:100)
#' {
#'  data <- raster::raster(matrix(runif(100,0,1), ncol=10))
#'  assign(paste('layer', i , sep='_'), data)
#' }
#' my.list <- lapply(ls(pattern='layer_'), get)
#' data <- do.call(stack, my.list)
#' # define a, b, and c
#' a <- nlayers(data)/2
#' b <- a+1
#' c <- nlayers(data)
#' a;b;c
#' # conduct two-sided Wilcox significance test
#' pvalue <- calc(data, intFun.grid.wilcox)
#' plot(pvalue); text(data, digits=2)
#'
#' @keywords internal
#' @export
intFun.grid.wilcox <- function(x) {
    a <- length(x)/2
    b <- a + 1
    c <- length(x)
    y <- stats::wilcox.test(x[1:a], x[b:c], alternative = c("two.sided"))$p.value
    return(y)
}

#-------------------------------------------------------------------------------

# intFun.site.points
#' Make spatial points
#' @description This function creates a spatial points data frame for a regular gaussian grid.
#' @param lon A vector (longitude)
#' @param lat A vector (latitude)
#' @param data A vector (data values)
#' @return A spatial points data frame
#' @examples
#'
#' library(raster)
#' # make some data
#' lon <- runif(100,-180,180)
#' lat <- runif(100, -90, 90)
#' data <- runif(100, 0, 10)
#' # convert data to spatial points data frame
#' data <- intFun.site.points(lon, lat, data)
#' plot(data)
#'
#' @keywords internal
#' @export
intFun.site.points <- function(lon, lat, data) {
    regular <- "+proj=longlat +ellps=WGS84"
    points <- data.frame(lon, lat, data)
    sp::coordinates(points) <- ~lon + lat
    raster::projection(points) <- regular
    return(points)
}

#-------------------------------------------------------------------------------

# intFun.rmse
#' Root mean square error (RMSE)
#' @description This function computes the root mean square error (RMSE), which is defined as:
#'
#' \eqn{$rmse(\lambda, \phi)=\sqrt{\frac{1}{t_{f}-t_{0}}\int_{t_{0}}^{t_{f}}(v_{mod}(t,\lambda, \phi)-v_{ref}(t,\lambda, \phi))^{2}dt}$}
#'
#' where \eqn{\lambda} is the longitude, \eqn{\phi} is the latitude, \eqn{t}
#' is the time, \eqn{t_0} is the initial time step, \eqn{t_f} is the final time
#' time step, \eqn{v_{mod}} is a modelled variable and \eqn{v_{ref}} is the
#' corresponding reference variable.
#'
#' @param mod An R object (model output data)
#' @param ref An R object (reference data)
#' @return An R object that gives the root mean square error when comparing
#' \code{mod} against \code{ref}.
#' @examples
#'
#' library(raster)
#' # create two raster stacks
#' for(i in 1:100)
#' {
#'  mod <- raster::raster(matrix(runif(100,-1,1), ncol=10))
#'  ref <- raster::raster(matrix(runif(100,-2,2), ncol=10))
#'  assign(paste('mod', i , sep='_'), mod)
#'  assign(paste('ref', i , sep='_'), ref)
#' }
#' my.list.mod <- lapply(ls(pattern='mod_'), get)
#' my.list.ref <- lapply(ls(pattern='ref_'), get)
#' mod <- do.call(stack, my.list.mod)
#' ref <- do.call(stack, my.list.ref)
#' # compute RMSE
#' rmse <- intFun.rmse(mod,ref)
#' plot(rmse); text(rmse, digits=2)
#'
#' @keywords internal
#' @export
intFun.rmse <- function(mod, ref) {
    if (intFun.isRaster(ref) == TRUE)
        sqrt(raster::mean((mod - ref)^2, na.rm = TRUE)) else sqrt(mean((mod - ref)^2, na.rm = TRUE))
}

#-------------------------------------------------------------------------------

# intFun.anom
#' Anomalies
#' @description This function computes anomalies.
#' @param x An R object
#' @return An R object that gives the anomalies
#' @examples
#'
#' # make some data
#' month <- seq(1,12,1)
#' month <- rep(month,10)
#' mod <- runif(length(month), 0,100)
#' mod <- data.frame(month, mod)
#' # compute anomalies
#' mod.anom <- data.frame(apply(mod[2], 2, intFun.anom))
#'
#' @keywords internal
#' @export
intFun.anom <- function(x) {
    if (intFun.isRaster(x) == TRUE)
        x - raster::mean(x, na.rm = TRUE) else x - mean(x, na.rm = TRUE)
}

#-------------------------------------------------------------------------------

# intFun.crmse
#' Centralized root mean square error (CRMSE)
#' @description This function computes the centralized root mean square error,
#' which is defined as:
#'
#' \eqn{$crmse(\lambda, \phi) = \sqrt{\frac{1}{t_{f}-t_{0}}\int_{t_{0}}^{t_{f}}[(v_{mod}(t,\lambda, \phi)-\overline{v_{mod}}(\lambda, \phi))-(v_{ref}(t,\lambda, \phi)-\overline{v_{ref}}(\lambda, \phi))]^{2}dt}$}
#'
#' where \eqn{\lambda} is the longitude, \eqn{\phi} is the latitude, \eqn{t}
#' is the time, \eqn{t_0} is the initial time step, \eqn{t_f} is the final time
#' time step, \eqn{v_{mod}} is a modelled variable, \eqn{v_{ref}} is the
#' corresponding reference variable, \eqn{\overline{v_{mod}}} is the time-mean
#' modelled variable, and \eqn{\overline{v_{ref}}} is the time-mean reference
#' variable.
#' @param mod.anom An R object (e.g. monthly anomalies from model output)
#' @param ref.anom An R object (e.g. monthly anomalies from reference data)
#' @return An R object that shows the centralized root mean square error
#' @examples
#'
#' library(raster)
#' # create two raster stacks
#' for(i in 1:100)
#' {
#'  mod <- raster::raster(matrix(runif(100,0,10), ncol=10))
#'  ref <- raster::raster(matrix(runif(100,0,10), ncol=10))
#'  assign(paste('mod', i , sep='_'), mod)
#'  assign(paste('ref', i , sep='_'), ref)
#' }
#' my.list.mod <- lapply(ls(pattern='mod_'), get)
#' my.list.ref <- lapply(ls(pattern='ref_'), get)
#' mod <- do.call(stack, my.list.mod)
#' ref <- do.call(stack, my.list.ref)
#' # compute anomalies
#' mod.anom <- intFun.anom(mod)
#' ref.anom <- intFun.anom(ref)
#' # compute CRMSE
#' crmse <- intFun.crmse(mod.anom, ref.anom)
#' plot(crmse); text(crmse, digits=2)
#'
#' @keywords internal
#' @export
intFun.crmse <- function(mod.anom, ref.anom) {
    if (intFun.isRaster(ref.anom) == TRUE)
        sqrt(raster::mean((mod.anom - ref.anom)^2, na.rm = TRUE)) else sqrt(mean((mod.anom - ref.anom)^2, na.rm = TRUE))
}

#-------------------------------------------------------------------------------

# intFun.theta
#' Time distance in months
#' @description Consider two rasters that show the month of the maximum value of
#' a variable during the climatological mean annual cycle. Calculating the
#' absolute difference between both rasters may yield values that range from
#' 0 to 11 (type a). A more meaningful range that considers the circular nature
#' of the annual cycle would range from 0 to 6 (type b). This function converts
#' a raster of type (a) to type (b).
#' @param x A raster object of type (a)
#' @return A raster object of type (b)
#' @examples
#'
#' library(raster)
#' # create a raster object with NA
#' data <- runif(100,0,11)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' # replace NA with zero in raster object
#' data <- calc(data,intFun.theta)
#' plot(data); text(data)
#'
#' @keywords internal
#' @export
intFun.theta <- function(x) {
    ifelse(x <= 6, x, (12 - x))
}

#-------------------------------------------------------------------------------

# intFun.anom.mly
#' Monthly anomalies
#' @description This function computes monthly anomalies.
#' @param mly An R object with a monthly time series
#' @param clim An R object with climatological monthly mean data
#' @return A raster object with monthly anomalies
#' @examples
#'
#' # make some data
#' month <- seq(1,12,1)
#' month <- rep(month,10)
#' mod <- runif(length(month), 0,10)
#' mod <- data.frame(month, mod)
#' # make an index
#' index <- list(mod$month)
#' clim.mly <- apply(mod, 2, function(x) {tapply(x, index, mean, na.rm=TRUE)})
#' mod.anom <- intFun.anom.mly(mod, clim.mly)
#'
#' @keywords internal
#' @export
intFun.anom.mly <- function(mly, clim) {
    data <- merge(mly, clim, by = "month")
    n <- (ncol(data) - 1)/2 + 1
    mod <- data[2:n]
    clim <- data[(n + 1):ncol(data)]
    anom <- mod - clim
    return(anom)
}

#-------------------------------------------------------------------------------

# intFun.iav
#' Inter-annual variability
#' @description This function computes the inter-annual variability. All data must start in January. All months after the last Dec will be dropped if data does not end in December.
#' @param anom R object with monthly anomalies
#' @return R object with inter-annual variability
#' @examples
#'
#' # make some data
#' month <- seq(1,12,1)
#' month <- rep(month,10)
#' mod <- runif(length(month), 0,100)
#' mod <- data.frame(month, mod)
#' # compute climatological mean monthly values
#' index <- list(mod$month)
#' mod.clim.mly <- apply(mod, 2, function(x) {tapply(x, index, mean, na.rm=TRUE)})
# compute anomalies
#' mod.anom <- intFun.anom.mly(mod, mod.clim.mly)
# compute monthly inter-annual variability
#' mod.iav <- intFun.iav(mod.anom)
#'
#' @keywords internal
#' @export
intFun.iav <- function(anom) {
    sqrt(mean((anom)^2, na.rm = TRUE))
}

#-------------------------------------------------------------------------------

# intFun.min.max.int
#' Range and interval for color bar legend
#' @description This function returns the minimum, maximum, and interval value
#' that can be used to define a color bar legend.
#' @param x An R object
#' @return Minimum, maximum, and interval value for color bar legend
#' @examples
#'
#' library(raster)
#' library(classInt)
#' # create a raster object
#' data <- runif(100,-23,864)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' # Get min, max, and interval for color bar legend
#' mmi <- intFun.min.max.int(data)
#'
#' @keywords internal
#' @export
intFun.min.max.int <- function(x) {
    if (intFun.isRaster(x) == TRUE) {
        values <- raster::getValues(x)
    } else {
        values = x
    }
    values <- stats::na.omit(values)
    suppressWarnings(classes <- classInt::classIntervals(values, style = "pretty"))
    classes <- data.frame(classes[2])
    min <- min(classes[])
    max <- max(classes[])
    int <- classes[2, ] - classes[1, ]
    min <- round(min, 3)
    max <- round(max, 3)
    int <- round(int, 3)
    y <- c(min, max, int)
    return(y)
}

#-------------------------------------------------------------------------------

# intFun.min.max.int.ext
#' Range and interval for color bar legend
#' @description This function returns the minimum, maximum, and interval value
#' that can be used to define a color bar legend. The minimum and maximum values are computed for the first and 99th percentile.
#' @param x An R object
#' @return Minimum, maximum, and interval value for color bar legend
#' @examples
#'
#' library(raster)
#' library(classInt)
#' # create a raster object
#' data <- runif(100,-23,864)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' # Get min, max, and interval for color bar legend
#' mmi <- intFun.min.max.int.ext(data)
#'
#' @keywords internal
#' @export
intFun.min.max.int.ext <- function(x) {
    if (intFun.isRaster(x) == TRUE) {
        values <- raster::getValues(x)
    } else {
        values = x
    }
    values <- stats::na.omit(values)

    # exclude extremes

    qUpper <- stats::quantile(values, probs = 0.99)
    qLower <- stats::quantile(values, probs = 0.01)
    values[values > qUpper] <- NA
    values[values < qLower] <- NA
    values <- stats::na.omit(values)

    suppressWarnings(classes <- classInt::classIntervals(values, style = "pretty"))
    classes <- data.frame(classes[2])
    min <- min(classes[])
    max <- max(classes[])
    int <- classes[2, ] - classes[1, ]
    min <- round(min, 3)
    max <- round(max, 3)
    int <- round(int, 3)
    y <- c(min, max, int)
    return(y)
}

#-------------------------------------------------------------------------------

# intFun.min.max.int.mod.ref
#' Range and interval for color bar legend for two raster objects
#' @description This function returns the minimum, maximum, and interval value
#' that can be used to define a color bar legend for two raster objects. This is
#' useful when plotting a variable from different data sets (e.g. model and
#' observation-based data).
#' @param mod An R object
#' @param ref An R object
#' @return Minimum, maximum, and interval value for color bar legend
#' @examples
#'
#' library(raster)
#' library(classInt)
#' # create a raster object
#' data <- runif(100,-23,864)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' mod <- data
#' ref <- data+300
#' # Get min, max, and interval for color bar legend
#' mmi <- intFun.min.max.int.mod.ref(mod, ref)
#'
#' @keywords internal
#' @export
intFun.min.max.int.mod.ref <- function(mod, ref) {
    if (intFun.isRaster(ref) == TRUE) {
        values <- c(raster::getValues(mod), raster::getValues(ref))
    } else {
        values <- unlist(c(mod, ref))
    }
    values <- stats::na.omit(values)
    suppressWarnings(classes <- classInt::classIntervals(values, style = "pretty"))
    classes <- data.frame(classes[2])
    min <- min(classes[])
    max <- max(classes[])
    int <- classes[2, ] - classes[1, ]
    min <- round(min, 3)
    max <- round(max, 3)
    int <- round(int, 3)
    y <- c(min, max, int)
    return(y)
}

#-------------------------------------------------------------------------------

# intFun.min.max.int.bias
#' Range and interval for color bar legend (bias)
#' @description This function returns the minimum, maximum, and interval value
#' that can be used to define a color bar legend for plotting biases. The
#' function ensures that the absolute values of the minimum and maximum are
#' identical. For instance, consider a bias that ranges from -10 to 300.
#' The resulting minimum and maximum value used by the color scheme are then
#' -300 and +300, respectively. The function excludes extremes identified as 1st
#' and 99th percentiles.
#' @param x An R object
#' @return Minimum, maximum, and interval value for color bar legend used for
#' plotting biases.
#' @examples
#'
#' library(raster)
#' library(classInt)
#' # create a raster object
#' data <- runif(100,-23,864)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' # Get min, max, and interval for color bar legend
#' mmi <- intFun.min.max.int.bias(data)
#'
#' @keywords internal
#' @export
intFun.min.max.int.bias <- function(x) {
    if (intFun.isRaster(x) == TRUE) {
        values <- raster::getValues(x)
    } else {
        values = x
    }
    values <- stats::na.omit(values)

    # exclude extremes
    qUpper <- stats::quantile(values, probs = 0.99)
    qLower <- stats::quantile(values, probs = 0.01)
    values[values > qUpper] <- NA
    values[values < qLower] <- NA
    values <- stats::na.omit(values)

    max.abs <- max(abs(values))
    values <- c(-max.abs, max.abs, values)
    suppressWarnings(classes <- classInt::classIntervals(values, style = "pretty"))
    classes <- data.frame(classes[2])
    min <- min(classes[])
    max <- max(classes[])
    int <- classes[2, ] - classes[1, ]
    min <- round(min, 3)
    max <- round(max, 3)
    int <- round(int, 3)
    y <- c(min, max, int)
    return(y)
}

#-------------------------------------------------------------------------------

# intFun.min.max.int.diff
#' Range and interval for color bar legend (difference)
#' @description This function returns the minimum, maximum, and interval value
#' that can be used to define a color bar legend for plotting differences. The
#' function ensures that the absolute values of the minimum and maximum are
#' identical. For instance, consider a bias that ranges from -10 to 300.
#' The resulting minimum and maximum value used by the color scheme are then
#' -300 and 300, respectively. The function is similar to intFun.min.max.int.bias,
#' but does not exclude extremes.
#' @param x An R object
#' @return Minimum, maximum, and interval value for color bar legend used for
#' plotting biases.
#' @examples
#'
#' library(raster)
#' library(classInt)
#' # create a raster object
#' data <- runif(100,-23,864)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' # Get min, max, and interval for color bar legend
#' mmi <- intFun.min.max.int.diff(data)
#'
#' @keywords internal
#' @export
intFun.min.max.int.diff <- function(x) {
    if (intFun.isRaster(x) == TRUE) {
        values <- raster::getValues(x)
    } else {
        values = x
    }
    values <- stats::na.omit(values)
    max.abs <- max(abs(values))
    values <- c(-max.abs, max.abs, values)
    suppressWarnings(classes <- classInt::classIntervals(values, style = "pretty"))
    classes <- data.frame(classes[2])
    min <- min(classes[])
    max <- max(classes[])
    int <- classes[2, ] - classes[1, ]
    min <- round(min, 3)
    max <- round(max, 3)
    int <- round(int, 3)
    y <- c(min, max, int)
    return(y)
}

#-------------------------------------------------------------------------------

# intFun.min.max.int.raw
#' Range and interval for color bar legend
#' @description This function returns the minimum, maximum, and interval value
#' that can be used to define a color bar legend. Contrary to \link{intFun.min.max.int},
#' this function does not use any rounding to find pretty values for min, max, int,
#' which is an advantage for plotting data with extremely small numerical values.
#' @param x An R object
#' @return Minimum, maximum, and interval value for color bar legend
#' @examples
#'
#' library(raster)
#' library(classInt)
#' # create a raster object
#' data <- runif(100,0,10^(-5))
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' # Get min, max, and interval for color bar legend
#' mmi <- intFun.min.max.int.raw(data)
#'
#' @keywords internal
#' @export
intFun.min.max.int.raw <- function(x) {
    if (intFun.isRaster(x) == TRUE) {
        values <- raster::getValues(x)
    } else {
        values = x
    }
    values <- stats::na.omit(values)
    suppressWarnings(classes <- classInt::classIntervals(values, style = "pretty"))
    classes <- data.frame(classes[2])
    min <- min(classes[])
    max <- max(classes[])
    int <- classes[2, ] - classes[1, ]
    y <- c(min, max, int)
    return(y)
}

#-------------------------------------------------------------------------------

# intFun.bin
#' Make bins
#' @description This function converts values into bins.
#' @param x A raster object
#' @param bin A number
#' @return A raster object with binned values
#' @examples
#'
#' library(raster)
#' data <- raster::raster(matrix(runif(100,0,100), ncol=10))
#' my.bin <- 25
#' # convert values to bins
#' bins <- intFun.bin(data, my.bin)
#' # plot data
#' plot(data);text(data)
#' plot(bins); text(bins)
#'
#' @keywords internal
#' @export
intFun.bin <- function(x, bin) {
    bin * round(x/bin)
}

#-------------------------------------------------------------------------------

# intFun.rel.error
#' Relative error used for binned data
#' @description This function computes the relative error for binned data used
#' for assessing functional relationships:
#' \eqn{$\varepsilon_{func}^u=\sqrt{\frac{\int(f_{mod}(u)-f_{ref}(u))^2du}{\int f_{ref}(u)^2du}}$}
#'
#' where \eqn{f_{mod}(u)} and \eqn{f_{ref}(u)} are the binned time means of
#' the model and reference data, respectively.
#'
#' @param mod An R object (binned time means of model output)
#' @param ref An R object (binned time means of reference data)
#' @return An R object
#' @examples
#'
#' # make some data
#' mod <- runif(100,0,10)
#' ref <- runif(100,0,10)
#' # compute relative error
#' r <- intFun.rel.error(mod,ref)
#'
#' @keywords internal
#' @export
intFun.rel.error <- function(mod, ref) sqrt(sum((mod - ref)^2)/sum(ref^2))

#-------------------------------------------------------------------------------

# intFun.coast
#' Reproject coastline
#' @description This function reprojects the coastline for a specified domain.
#' @param my.xlim An R object with minimum and maximum longitude, e.g. c(-171, -23)
#' @param my.ylim An R object with minimum and maximum latitude, e.g. c(32, 75)
#' @param my.projection An R string that defines the desired projection, e.g. '+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.'
#' @param shp.filename An R string that gives the name of a shapefile that should be reprojected
#' @return Reprojected coastline
#' @examples
#'
#' library(rgdal)
#' library(rgeos)
#' my.xlim <- c(-171, -23)
#' my.ylim <- c(32, 75)
#' my.projection <- '+proj=ob_tran +o_proj=longlat +o_lon_p=83. +o_lat_p=42.5 +lon_0=263.'
#' shp.filename <- system.file('extdata/ne_110m_land/ne_110m_land.shp', package = 'amber')
#' land <- intFun.coast(my.xlim, my.ylim, my.projection, shp.filename)
#' raster::plot(land)
#'
#' @keywords internal
#' @export
intFun.coast <- function(my.xlim, my.ylim, my.projection = "+proj=longlat +ellps=WGS84", shp.filename) {
    # regular projection
    regular <- "+proj=longlat +ellps=WGS84"
    # make a box
    lon <- seq(min(my.xlim), max(my.xlim), 0.1)
    lat <- seq(min(my.ylim), max(my.ylim), 0.1)
    top <- cbind(lon, rep(max(lat), length(lon)))
    right <- cbind(rep(max(lon), length(lat)), rev(lat))
    bottom <- cbind(rev(lon), rep(min(lat), length(lon)))
    left <- cbind(rep(min(lon), length(lat)), lat)
    box <- rbind(top, right, bottom, left)
    box <- sp::Polygon(box)
    box <- sp::Polygons(list(box), 1)
    box <- sp::SpatialPolygons(list(box))
    suppressWarnings(raster::projection(box) <- regular)
    # coastline
    land <- raster::shapefile(shp.filename)
    suppressWarnings(raster::projection(land) <- regular)
    suppressWarnings(land <- rgeos::gBuffer(land, width = 0))  # this avoids a Ring Self-intersection error
    land <- raster::crop(land, box)
    suppressWarnings(land <- sp::spTransform(land, sp::CRS(my.projection)))
    return(land)
}

#-------------------------------------------------------------------------------

# random sampling in time

#-------------------------------------------------------------------------------

# intFun.sampleGridTime
#' Random sampling of gridded reference data along time axis
#' @description This function conducts random sampling
#' @param ref A raster object
#' @return A raster object
#' @examples
#'
#' library(raster)
#' # make some data
#' data01 <- c(seq(1,99,1), NA)
#' data02 <- c(NA, seq(101,199,1))
#' data03 <- data01 + data02
#' data01 <- raster::raster(matrix(data01, ncol=10))
#' data02 <- raster::raster(matrix(data02, ncol=10))
#' data03 <- raster::raster(matrix(data03, ncol=10))
#' data <- raster::stack(data01, data02, data03)
#'
#' randomSample <- raster::calc(data, fun = intFun.sampleGridTime)
#'
#' plot(data)
#' plot(randomSample)
#'
#' @keywords internal
#' @export
intFun.sampleGridTime <- function(ref) {
    mask.ref <- ref - ref + 1
    mySize <- length(ref)
    mySample <- sample(ref, size = mySize, replace = TRUE)
    y <- mySample
    y <- mySample * mask.ref
    return(y)
}
#-------------------------------------------------------------------------------

# intFun.sampleGridTimeAndSpace
#' random sampling across space and time (if available) for reference grid data
#' @description This function conducts random sampling with replacement for a single raster layer
#' @param ref A raster object
#' @return A raster object
#' @examples
#'
#' library(raster)
#' # make some data
#' data01 <- c(seq(1,99,1), NA)
#' data02 <- c(NA, seq(101,199,1))
#' data03 <- data01 + data02
#' data01 <- raster::raster(matrix(data01, ncol=10))
#' data02 <- raster::raster(matrix(data02, ncol=10))
#' data03 <- raster::raster(matrix(data03, ncol=10))
#' data <- raster::stack(data01, data02, data03)
#'
#' randomSample <- intFun.sampleGridTimeAndSpace(ref = data)
#' plot(data)
#' plot(randomSample)
#'
#' @keywords internal
#' @export
intFun.sampleGridTimeAndSpace <- function(ref) {
    mask.ref <- ref - ref + 1
    refValues <- raster::values(ref)
    mySize <- length(refValues)
    refValues <- matrix(refValues, ncol = 1)
    refValues <- stats::na.omit(refValues)
    mySample <- sample(refValues, size = mySize, replace = TRUE)
    y <- raster::setValues(ref, mySample)
    y <- y * mask.ref
    return(y)
}

#-------------------------------------------------------------------------------

# intFun.sampleDataFrame
#' random sampling across space and time for a data frame (FLUXNET and runoff).
#' @description This function conducts random sampling with replacement for a single raster layer
#' @param ref A raster object
#' @return A raster object
#' @examples
#'
#' # make some data
#' data <- seq(1, 29, 1)
#' data <- c(data, NA)
#' data <- data.frame(matrix(data, ncol = 10))
#'
#' randomSample <- intFun.sampleDataFrame(ref = data)
#' print(data)
#' print(randomSample)
#'
#' @keywords internal
#' @export
intFun.sampleDataFrame <- function(ref) {
    mask <- ref - ref + 1
    refValues <- unlist(ref)
    mySize <- length(refValues)
    refValues <- data.frame(refValues)
    refValues <- stats::na.omit(refValues)
    refValues <- unlist(refValues)
    refSample <- sample(refValues, size = mySize, replace = TRUE)
    refSample <- matrix(refSample, ncol = ncol(ref))
    refSample <- data.frame(refSample)
    refSample <- refSample * mask
    y <- refSample
    return(y)
}

#-------------------------------------------------------------------------------

# intFun.addBWtext
#' Add text to raster
#' @description This function adds text to raster. A threshold value determines
#' whether the text is black or white
#' @param myRaster A raster object
#' @param myDigits An integer that gives the number of desired digits
#' @param myDigits A number that determines the size of the text, e.g. 0.7
#' @return plots text of a raster object
#' @examples
#'
#' library(raster)
#' # make some data
#' data <- runif(100,-1,1)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#'
#' plot(data)
#' intFun.addBWtext(myRaster = data, myDigits = 1, myCex = 0.7)
#'
#' @keywords internal
#' @export
intFun.addBWtext <- function(myRaster, myDigits = 1, myCex = 0.7) {
    myThreshold <- max(abs(raster::values(myRaster)), na.rm = TRUE)/2
    if (min(abs(raster::values(myRaster)), na.rm = TRUE) < myThreshold) {
        raster::text(myRaster, digits = myDigits, fun = function(x) {
            abs(x) < myThreshold
        }, cex = myCex)
    }
    if (max(abs(raster::values(myRaster)), na.rm = TRUE) >= myThreshold) {
        raster::text(myRaster, digits = myDigits, fun = function(x) {
            abs(x) >= myThreshold
        }, col = "white", cex = myCex)
    }
}

#-------------------------------------------------------------------------------
# Functions copied from other packages to reduce dependencies:
#-------------------------------------------------------------------------------
# intFun.isRaster
#' Reproject coastline
#' @description This function assesses whether an R object is a raster. The original
#' code was copied from intFun.isRaster (spatial.tools_1.6.0)
#' @param x An R object such as a raster or a number
#' @examples
#'
#' x <- 1
#' intFun.isRaster(x)
#' y <- raster::raster(matrix(seq(1,10), ncol = 2))
#' intFun.isRaster(y)
#'
#' @keywords internal
#' @export
intFun.isRaster <- function(x) {
    return((class(x) == "RasterLayer" || class(x) == "RasterBrick" || class(x) == "RasterStack"))
}

#-------------------------------------------------------------------------------
