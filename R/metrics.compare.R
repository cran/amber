################################################################################
#' Compares global mean statistical meterics for multiple model runs
#' @description This function plots statistical meterics for multiple model runs.
#' This is useful for comparing the impact of changes in model settings or input data.
#' @param mod.path.list A list with paths for each model run, e.g. mod.path.list <- list(mod01.path, mod02.path, mod03.path)
#' @param plot.width Number that gives the plot width, e.g. 6
#' @param plot.height Number that gives the plot height, e.g. 5
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @param defineVariableOrder Logical. If TRUE, variables are sorted according to the parameter myVariables defined below. Default setting is FALSE.
#' @param myVariables An R object with variable names of variables that should be included in table, e.g. c('GPP', 'RECO', 'NEE')
#' @return A figure in PDF format that shows statistical metrics for multiple model runs.
#' @examples
#'
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
#' mod01.path <- system.file('extdata/SIMod01', package = 'amber')
#' mod02.path <- system.file('extdata/SIMod02', package = 'amber')
#'
#' mod.path.list <- list(mod01.path, mod02.path)
#' myVariables <- c('RNS', 'RSS', 'RLS', 'ALBS', 'HFLS', 'HFSS', 'HFG', 'GPP',
#' 'RECO', 'NEE', 'FIRE', 'AGB', 'CVEG', 'CSOIL', 'LAI', 'BURNT', 'SNW', 'MRSLL', 'MRRO')
#'
#' myVariables <- c('ALBS', 'GPP', 'LAI')
#'
#' metrics.compare(mod.path.list = mod.path.list, plot.width = 5, plot.height = 15,
#' outputDir = FALSE, defineVariableOrder = TRUE, myVariables = myVariables)
#'
#' @export
metrics.compare <- function(mod.path.list = mod.path.list, plot.width = 13, plot.height = 14, outputDir = FALSE, defineVariableOrder = TRUE,
    myVariables = myVariables) {

    nmod <- length(mod.path.list)  # number of model runs

    result <- foreach::foreach(i = 1:nmod) %do% {
        if (defineVariableOrder == FALSE) {
            # get variable and model and reference data name
            my.list <- list.files(path = mod.path.list[[i]], pattern = "scoreinputs_")
            split <- strsplit(my.list, "_")
            split <- unlist(split)
            split <- matrix(split, nrow = 9)
            split <- t(split)
            split <- as.data.frame(split)
            mod.id <- split[3]
        }

        #---------------------------------------------------------------------------
        # The user may want to define the exact order of variables, which is done here.  Otherwise, all variables are listed
        # alphabetically.
        #---------------------------------------------------------------------------
        if (defineVariableOrder == TRUE) {
            # define order
            myOrder <- seq(1, length(myVariables), 1)
            myOrder <- data.frame(myVariables, myOrder)
            colnames(myOrder) <- c("variable", "order")
            # get variable and model and reference data name
            my.list <- list.files(path = mod.path.list[[i]], pattern = "scoreinputs_")
            split <- strsplit(my.list, "_")
            split <- unlist(split)
            split <- matrix(split, nrow = 9)
            split <- t(split)
            split <- as.data.frame(split)
            split <- split[c(2, 3, 5)]  # select columns
            split <- data.frame(split, unlist(my.list))
            colnames(split) <- c("variable", "mod.id", "ref.id", "fileName")
            split <- merge(split, myOrder, by = "variable")
            mod.id <- split$mod.id
            split <- split[order(split$order, split$ref.id), ]
            my.list <- as.character(split$fileName)
        }
        #---------------------------------------------------------------------------

        # get data (statistical metrics for each variable)
        data.list <- lapply(paste(mod.path.list[[i]], my.list, sep = "/"), utils::read.table)

        # convert list to matrix
        data <- do.call("rbind", data.list)
        data <- data.frame(mod.id, data)
    }

    # The goal is to create a matrix like this: bias model 1, bias model 2, rmse model03, (...)  variable 1 variable 2 (...)

    # Check whether all models have an identical set of variables and reference data
    variable.name.mod01 <- result[[1]]$variable.name
    ref.id.mod01 <- result[[1]]$ref.id
    checkVariable <- lapply(result, function(x) identical(variable.name.mod01, x$variable.name))
    checkRef <- lapply(result, function(x) identical(ref.id.mod01, x$ref.id))
    lapply(checkVariable, function(x) stopifnot(checkVariable == TRUE))
    lapply(checkVariable, function(x) stopifnot(checkRef == TRUE))

    # extract columns of interest
    variable.name <- result[[1]]$variable.name
    variable.unit <- result[[1]]$variable.unit
    mod.id <- lapply(result, "[", , 1)
    mod.id <- unique(unlist(mod.id))
    ref.id <- result[[1]]$ref.id
    ref.mean <- result[[1]]$ref.mean.scalar
    mod.mean <- lapply(result, "[", , 6)
    bias.abs <- lapply(result, "[", , 8)
    bias.rel <- lapply(result, "[", , 9)
    epsilon_rmse <- lapply(result, "[", , 16)
    phase <- lapply(result, "[", , 20)
    epsilon_iav <- lapply(result, "[", , 24)
    sigma <- lapply(result, "[", , 28)
    R <- lapply(result, "[", , 29)

    variable.name <- data.frame(variable.name)
    variable.unit <- data.frame(variable.unit)
    ref.id <- data.frame(ref.id)
    ref.mean <- data.frame(ref.mean)
    mod.mean <- data.frame(mod.mean)
    bias.abs <- data.frame(bias.abs)
    bias.rel <- data.frame(bias.rel)
    epsilon_rmse <- data.frame(epsilon_rmse)
    phase <- data.frame(phase)
    epsilon_iav <- data.frame(epsilon_iav)
    sigma <- data.frame(sigma)
    R <- data.frame(R)
    #-------------------------------------------------------------------------------

    colnames(mod.mean) <- paste("mean", mod.id, sep = "_")
    colnames(bias.abs) <- paste("bias.abs", mod.id, sep = "_")
    colnames(bias.rel) <- paste("bias", mod.id, sep = "_")
    colnames(epsilon_rmse) <- paste("epsilon_rmse", mod.id, sep = "_")
    colnames(phase) <- paste("phase", mod.id, sep = "_")
    colnames(epsilon_iav) <- paste("epsilon_iav", mod.id, sep = "_")
    colnames(sigma) <- paste("sigma", mod.id, sep = "_")
    colnames(R) <- paste("R", mod.id, sep = "_")

    data <- data.frame(variable.name, variable.unit, ref.id, ref.mean, mod.mean, bias.abs, bias.rel, epsilon_rmse, phase, epsilon_iav,
        sigma, R)
    units <- data$variable.unit

    # Omit backslashes in units (e.g. fractional area burnt is in '\\% per month').
    fun.omit.backslash <- function(x, y) {
        y <- gsub("\\", "", x, fixed = TRUE)
        return(y)
    }
    units <- sapply(X = units, FUN = fun.omit.backslash)

    units <- latex2exp::TeX(units)
    variable.ref <- paste(data$variable.name, data$ref.id, sep = "-")
    my.lab.metrics <- latex2exp::TeX(c("mean (units)", "bias (abs)", "bias (\\%)", "$\\epsilon_{rmse} (-)$", "$\\theta$ (months)",
        "$\\epsilon_{iav}$ (-)", "$\\sigma$ (-)", "$R$ (-)"))
    data <- data[-c(1, 2, 3)]  # drop first three columns, i.e. variable name, unit, and ref.id
    my.lab.models <- rep(mod.id, (ncol(data) - 1)/length(mod.id))
    my.lab.models <- c("Reference", as.character(my.lab.models))

    #-------------------------------------------------------------------------------

    # convert to raster

    data <- as.matrix(data)
    data <- raster::raster(data)

    # raster extent
    my.width <- 2
    my.extent <- c(0, my.width * ncol(data), 0.5, nrow(data) + 0.5)
    raster::extent(data) <- my.extent

    data.all <- data

    # split raster by statistical metric for the use of different color schemes

    # ref.mean, mod.mean, bias.abs, bias.rel, epsilon_rmse, phase, epsilon_iav, sigma, R

    firstChunk <- 1 + 2 * nmod  # number of rows with anbsolute numbers (1 x ref.mean, n.mod x mod.mean, n.mod x bias.abs)

    data <- data.all
    data[raster::cellFromCol(data, c((firstChunk + 1):ncol(data.all)))] <- NA
    refModBias <- data  # 4 columns

    data <- data.all
    data[raster::cellFromCol(data, c(1:firstChunk, (3 * nmod + 2):ncol(data.all)))] <- NA
    bias.rel <- data

    data <- data.all
    data[raster::cellFromCol(data, c(1:(3 * nmod + 1), (4 * nmod + 2):ncol(data.all)))] <- NA
    epsilon_rmse <- data

    data <- data.all
    data[raster::cellFromCol(data, c(1:(4 * nmod + 1), (5 * nmod + 2):ncol(data.all)))] <- NA
    phase <- data

    data <- data.all
    data[raster::cellFromCol(data, c(1:(5 * nmod + 1), (6 * nmod + 2):ncol(data.all)))] <- NA
    epsilon_iav <- data

    data <- data.all
    data[raster::cellFromCol(data, c(1:(6 * nmod + 1), (7 * nmod + 2):ncol(data.all)))] <- NA
    sigma <- data

    data <- data.all
    data[raster::cellFromCol(data, c(1:(7 * nmod + 1)))] <- NA
    R <- data

    myList <- list(refModBias, bias.rel, epsilon_rmse, phase, epsilon_iav, sigma, R)

    # plot

    # colors and breaks (colBre) this is a nested list colBre[[i]][[1]] gives colors colBre[[i]][[2]] gives color breaks
    colBre <- foreach::foreach(i = 1:length(myList)) %do% {
        data <- myList[[i]]
        mmi <- intFun.min.max.int(data)
        if (i == 2) {
            mmi <- intFun.min.max.int.diff(data)
        }
        min <- mmi[1]
        max <- mmi[2]
        interval <- mmi[3]
        my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
        my.labels <- round(seq(min, max, interval), 3)  # location of labels
        my.col <- rev(viridis::viridis(n = length(my.labels) - 1))

        if (i == 2) {
            my.col <- scico::scico(n = length(my.labels) - 1, palette = "vik")
        }

        if (i == 1) {
            my.col <- NA
            my.breaks <- c(NA, NA)
        }

        my.col <- list(my.col, my.breaks)
    }

    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))

    if (outputDir != FALSE) {
        grDevices::pdf("metrics_compare.pdf", width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(13, 10, 5, 5), lwd = 1, cex = 1, family = "Helvetica")

    # plot data
    raster::plot(data.all, col = "grey90", legend = FALSE, main = NA, axes = FALSE)
    foreach::foreach(i = 1:length(myList)) %do% {
        data <- myList[[i]]
        my.col <- colBre[[i]][[1]]
        my.breaks <- colBre[[i]][[2]]
        graphics::par(new = TRUE)
        raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = NA, axes = FALSE)
        if (i == 1) {
            raster::text(data, digits = 2, cex = 0.7)
        }
        if (i == 2) {
            intFun.addBWtext(myRaster = data, myDigits = 0, myCex = 0.7)
        }
        if (i > 2) {
            intFun.addBWtext(myRaster = data, myDigits = 2, myCex = 0.7)
        }
    }

    graphics::axis(side = 1, at = seq(1, my.width * ncol(data.all), my.width), labels = my.lab.models, las = 2)
    graphics::axis(side = 2, at = rev(seq(1, nrow(data.all), 1)), labels = variable.ref, las = 2)
    graphics::axis(side = 4, at = rev(seq(1, nrow(data.all), 1)), labels = units, las = 2)
    # graphics::axis(side = 3, at = my.width + seq(my.width * nmod/2, my.width * ncol(data.all), my.width * nmod), labels =
    # my.lab.metrics, las = 1, tcl = 0)

    a <- my.width * nmod/2
    b <- my.width * (ncol(data.all) - 1)
    c <- my.width * nmod

    graphics::axis(side = 3, at = my.width + seq(a, b, c), labels = my.lab.metrics, las = 1, tcl = 0)

    # graphics::abline(v = my.width)
    graphics::abline(v = my.width + 1 * my.width * nmod)
    graphics::abline(v = my.width + 2 * my.width * nmod)
    graphics::abline(v = my.width + 3 * my.width * nmod)
    graphics::abline(v = my.width + 4 * my.width * nmod)
    graphics::abline(v = my.width + 5 * my.width * nmod)
    graphics::abline(v = my.width + 6 * my.width * nmod)
    graphics::abline(v = my.width + 7 * my.width * nmod)

    # add horizontal lines between the different variables
    hoLiLo <- table(variable.name)  # horizontal line location
    hoLiLo <- rev(hoLiLo)
    hoLiLo <- cumsum(hoLiLo) + 0.5
    foreach::foreach(i = 1:length(hoLiLo)) %do% {
        graphics::abline(h = hoLiLo[i])
    }

    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}
