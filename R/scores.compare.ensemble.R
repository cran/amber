################################################################################
#' Summarize scores from multiple model runs in single figure.
#' @description This function produces a figure that summarizes score values from
#' multiple model runs. The figure has four columns, which give the multi-model
#' mean scores, the total score range, the model with the lowest score and the
#' model with the highest score. The respective inputs are created by the functions
#' \link{scores.fluxnet.csv} or \link{scores.fluxnet.nc},
#' \link{scores.grid.notime},
#' \link{scores.grid.time}, and
#' \link{scores.runoff}.
#' @param mod.path.list A list with paths for each model run, e.g. mod.path.list <- list(mod01.path, mod02.path, mod03.path)
#' @param modelIDs An R object with the different model run IDs, e.g. c('CLASSIC.CRUJRAv2', 'CLASSIC.GSWP3W5E5', 'CLASSIC.CRUNCEP')
#' @param myVariables An R object with variable names of variables that should be included in table, e.g. c('GPP', 'RECO', 'NEE')
#' @param plot.width Number that gives the plot width, e.g. 6
#' @param plot.height Number that gives the plot height, e.g. 5
#' @param myMargin An R object that gives the figure margins, e.g. c(4, 13, 3, 4)
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @return A figure in PDF format that shows the ensemble scores,
#' the total score range, the model with the lowest score and the model with the highest score.
#' @examples
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
#' mod01.path <- system.file('extdata/model01', package = 'amber')
#' mod02.path <- system.file('extdata/model02', package = 'amber')
#' mod.path.list <- list(mod01.path, mod02.path)
#' modelIDs <- c('CLASSIC.CRUJRAv2', 'CLASSIC.GSWP3W5E5')
#' myVariables <- c('GPP', 'BURNT')
#' #myVariables <- c('RNS', 'RSS', 'RLS', 'ALBS', 'HFLS', 'HFSS', 'HFG', 'GPP', 'RECO',
#' #'NEE', 'FIRE', 'AGB', 'CVEG', 'CSOIL', 'LAI', 'BURNT', 'SNW', 'MRSLL', 'MRRO')
#'
#' scores.compare.ensemble(mod.path.list = mod.path.list, modelIDs = modelIDs,
#' myVariables = myVariables, plot.width = 9.3, plot.height = 10,
#' myMargin = c(12, 0, 3, 0), outputDir = FALSE)
#'
#' @export
scores.compare.ensemble <- function(mod.path.list = mod.path.list, modelIDs = modelIDs, myVariables = myVariables, plot.width = 10,
    plot.height = 10, myMargin = c(12, 0, 3, 0), outputDir = FALSE) {


    nmod <- length(mod.path.list)  # number of model runs
    eachModel <- foreach::foreach(m = 1:nmod) %do% {

        # get a list of file names
        my.list <- list.files(path = mod.path.list[[m]], pattern = "scorevalues_")

        # define order
        myOrder <- seq(1, length(myVariables), 1)
        myOrder <- data.frame(myVariables, myOrder)
        colnames(myOrder) <- c("variable", "order")

        # get file names and sort according to defined order
        split <- strsplit(my.list, "_")
        split <- unlist(split)
        split <- matrix(split, nrow = 9)
        split <- t(split)
        split <- as.data.frame(split)
        split <- split[c(2, 5)]  # select columns
        split <- data.frame(split, unlist(my.list))
        colnames(split) <- c("variable", "ref.id", "fileName")
        split <- merge(split, myOrder, by = "variable")
        split <- split[order(split$order, split$ref.id), ]
        my.list <- as.character(split$fileName)

        # create vertical axis labels
        split <- split[c(1, 2)]  # select columns
        my.row.lab <- paste(split$variable, split$ref.id, sep = "-")

        # get data by looping through all models

        my.files <- paste(mod.path.list[[m]], my.list, sep = "/")
        data <- lapply(my.files, utils::read.table)
        data <- do.call("rbind", data)
        weights <- rep(c("not.weighted", "weighted"), nrow(data)/2)
        data <- data.frame(weights, data)

        row.names(data) <- c()  # drop row names
        data.nw <- subset(data, weights == "not.weighted")
        data.nw <- subset(data.nw, select = -c(weights))
        rownames(data.nw) <- c()  # omit rownames
        data <- data.nw

        # convert to raster
        data <- data[3:ncol(data)]
        colnames(data) <- c()
        data <- as.matrix(data)
        data <- raster::raster(data)
        raster::metadata(data) <- list(my.row.lab)
        my.width <- 2
        raster::extent(data) <- c(0, 6 * my.width, 0.5, nrow(data) + 0.5)
        data
    }
    data <- do.call(raster::stack, eachModel)
    score.mean <- raster::mean(data)
    score.min <- min(data)
    score.max <- max(data)
    score.range <- score.max - score.min
    mod.best <- raster::which.max(data)
    mod.worst <- raster::which.min(data)

    # get labels
    my.row.lab <- unlist(raster::metadata(data[[1]]))
    my.col.lab <- latex2exp::TeX(c("$S_{bias}$", "$S_{rmse}$", "$S_{phase}$", "$S_{iav}$", "$S_{dist}$", "$S_{overall}$"))

    #---------------------------------------------------------------------------

    # Figure

    #---------------------------------------------------------------------------

    my.filename <- "scores_compare_ensemble"

    # par
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, ".pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(1, 4), font.main = 1, oma = c(0, 8, 0, 0), mar = myMargin, lwd = 1, cex = 1, family = "Helvetica")

    # (a) ensemble scores
    data <- score.mean
    min <- 0
    max <- 1
    interval <- 0.1
    my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
    my.labels <- round(seq(min, max, interval), 3)  # location of labels
    my.col <- viridis::plasma(n = length(my.breaks) - 1, direction = -1)
    legend.bar.text <- "Score (-)"
    my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
    my.legend.args <- list(text = legend.bar.text, side = 3, font = 1, line = 1, cex = 0.7)

    raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = NA, axes = FALSE)
    intFun.addBWtext(myRaster = data, myDigits = 2, myCex = 0.7)
    graphics::title("(a)", line = 0.5, adj = 0, font.main = 2, cex.main = 1.25)
    raster::plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args,
        legend.width = 1, legend.shrink = 1, font = 1, horizontal = TRUE)

    graphics::axis(side = 1, at = seq(1, my.width * ncol(data), my.width), labels = my.col.lab, las = 2)
    # graphics::axis(side = 3, at = 3 * my.width, labels = ('(a) Mean ensemble scores'), las = 1, tcl = 0)
    graphics::axis(side = 2, at = rev(seq(1, nrow(data), 1)), labels = my.row.lab, las = 2)

    # (b) score ranges
    data <- score.range
    mmi <- intFun.min.max.int(data)
    min <- mmi[1]
    max <- mmi[2]
    interval <- mmi[3]
    my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
    my.labels <- round(seq(min, max, interval), 3)  # location of labels
    my.col <- viridis::viridis(n = length(my.breaks) - 1, direction = -1)
    legend.bar.text <- "Max. minus Min. Score (-)"
    my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
    my.legend.args <- list(text = legend.bar.text, side = 3, font = 1, line = 1, cex = 0.7)

    raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = NA, axes = FALSE)
    intFun.addBWtext(myRaster = data, myDigits = 2, myCex = 0.7)
    graphics::title("(b)", line = 0.5, adj = 0, font.main = 2, cex.main = 1.25)
    raster::plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args,
        legend.width = 1, legend.shrink = 1, font = 1, horizontal = TRUE)
    graphics::axis(side = 1, at = seq(1, my.width * ncol(data), my.width), labels = my.col.lab, las = 2)
    # graphics::axis(side = 3, at = 3 * my.width, labels = ('(b) Max. minus min. scores'), las = 1, tcl = 0)

    # (c) best model
    data <- mod.best
    min <- 0
    max <- nmod
    interval <- 1
    my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
    my.labels <- round(seq(min + 0.5, max - 0.5, interval), 3)  # location of labels
    my.col <- viridis::viridis(n = length(my.breaks) - 1, direction = -1)
    legend.bar.text <- "Ensemble member ID with best score"
    my.axis.args <- list(at = my.labels, labels = seq(1, nmod, 1), cex.axis = 1)
    my.legend.args <- list(text = legend.bar.text, side = 3, font = 1, line = 1, cex = 0.7)

    raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = NA, axes = FALSE)
    intFun.addBWtext(myRaster = data, myDigits = 0, myCex = 0.7)
    graphics::title("(c)", line = 0.5, adj = 0, font.main = 2, cex.main = 1.25)
    raster::plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args,
        legend.width = 1, legend.shrink = 1, font = 1, horizontal = TRUE)
    graphics::axis(side = 1, at = seq(1, my.width * ncol(data), my.width), labels = my.col.lab, las = 2)
    # graphics::axis(side = 3, at = 3 * my.width, labels = ('(c) Ensemble member with best score'), las = 1, tcl = 0)

    # (d) worst model
    legend.bar.text <- "Ensemble member ID with worst score"
    my.axis.args <- list(at = my.labels, labels = seq(1, nmod, 1), cex.axis = 1)
    my.legend.args <- list(text = legend.bar.text, side = 3, font = 1, line = 1, cex = 0.7)
    data <- mod.worst
    raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = NA, axes = FALSE)
    intFun.addBWtext(myRaster = data, myDigits = 0, myCex = 0.7)
    graphics::title("(d)", line = 0.5, adj = 0, font.main = 2, cex.main = 1.25)
    raster::plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args,
        legend.width = 1, legend.shrink = 1, font = 1, horizontal = TRUE)

    graphics::axis(side = 1, at = seq(1, my.width * ncol(data), my.width), labels = my.col.lab, las = 2)
    # graphics::axis(side = 3, at = 3 * my.width, labels = ('(d) Ensemble member with worst score'), las = 1, tcl = 0)

    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}

