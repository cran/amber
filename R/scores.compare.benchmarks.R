################################################################################
#' Compares model scores against scores reference scores.
#' @description Interpreting scores is challenging as reference data are subject to uncertainty.
#' This function compares two types of scores. The first set of scores expresses model performance
#' and is based on comparing model output against reference data. The second set of scores is based
#' on a comparison of two independent reference data (e.g. remotely sensed GPP against FLUXNET). The
#' difference between both scores reflect the uncertainty of reference data and indicates how well a model could
#' perform given this uncertainty. Scores that are based on reference data only are here referred to
#' as benchmark scores.
#'
#' @param bench.path A string that gives the path where benchmarks (i.e. reference vs. reference data) are stored
#' @param model.path A string that gives the path where the output from \link{scores.tables} is stored (model)
#' @param model.id A string that gives the name of a model, e.g. 'CLASSIC'
#' @param plot.width Number that gives the plot width, e.g. 7.3
#' @param plot.height Number that gives the plot height, e.g. 6.5
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @param defineVariableOrder Logical. If TRUE, variables are sorted according to the parameter myVariables defined below. Default setting is FALSE.
#' @param myVariables An R object that defines the variables and their order in the score table, e.g. c('GPP', 'RECO', 'NEE').
#' @return A figure in PDF format that shows the (a) benchmark skill score,
#' (b) model skill score, and (c) score difference.
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
#' bench.path <- system.file('extdata/scoresBenchmarks', package = 'amber')
#' model.path <- system.file('extdata/scores', package = 'amber')
#'
#' model.id <- 'CLASSIC'
#' myVariables <- c('ALBS', 'GPP', 'LAI')
#'
#' scores.compare.benchmarks(bench.path, model.path, model.id, plot.width = 8, plot.height = 4,
#' defineVariableOrder = TRUE, myVariables = myVariables)
#'
#' @export
scores.compare.benchmarks <- function(bench.path, model.path, model.id, plot.width = 8, plot.height = 4.8, outputDir = FALSE,
    defineVariableOrder = FALSE, myVariables = myVariables) {


    if (defineVariableOrder == FALSE) {
        # get variable and reference data name
        my.list <- list.files(path = bench.path, pattern = "scorevalues_")
        split <- strsplit(my.list, "_")
        split <- unlist(split)
        split <- matrix(split, ncol = length(my.list))
        split <- t(split)
        split <- as.data.frame(split)
        split <- split[c(2, 3, 5)]  # select columns
        colnames(split) <- c("variable.name", "mod.id", "ref.id")
    }

    #---------------------------------------------------------------------------
    # The user may want to define the exact order of variables, which is done here.  Otherwise, all variables are listed
    # alphabetically.
    #---------------------------------------------------------------------------
    # Get file names
    #---------------------------------------------------------------------------
    if (defineVariableOrder == TRUE) {
        # define order
        myOrder <- seq(1, length(myVariables), 1)
        myOrder <- data.frame(myVariables, myOrder)
        colnames(myOrder) <- c("variable.name", "order")
        # get file names and sort according to defined order
        my.list <- list.files(path = bench.path, pattern = "scorevalues_")
        split <- strsplit(my.list, "_")
        split <- unlist(split)
        split <- matrix(split, ncol = length(my.list))
        split <- t(split)
        split <- as.data.frame(split)
        split <- split[c(2, 3, 5)]  # select columns
        split <- data.frame(split, unlist(my.list))
        colnames(split) <- c("variable.name", "mod.id", "ref.id", "fileName")

        split <- merge(split, myOrder, by = "variable.name")
        split <- split[order(split$order, split$mod.id), ]
        my.list <- as.character(split$fileName)
        split <- split[c(1, 2, 3)]  # select columns
    }
    #---------------------------------------------------------------------------
    # Get benchmark scores
    #---------------------------------------------------------------------------

    # Get score values
    data.list <- lapply(paste(bench.path, my.list, sep = "/"), utils::read.table)

    # use not-weighted score values only
    data <- lapply(data.list, function(x) x["not.weighted", ])

    # convert list to matrix
    data <- do.call("rbind", data)
    rownames(data) <- c()

    data <- subset(data, select = -c(1, 2))
    myOrder <- seq(1, nrow(data), 1)
    df <- data.frame(split, data, myOrder)
    bench.scores <- df

    #---------------------------------------------------------------------------
    # Get model scores
    #---------------------------------------------------------------------------

    # Get score values
    nmod <- length(model.path)  # number of model runs
    eachModel <- foreach::foreach(m = 1:nmod) %do% {


        my.list <- list.files(path = model.path[[m]], pattern = "scorevalues_")
        data.list <- lapply(paste(model.path[[m]], my.list, sep = "/"), utils::read.table)

        # use not-weighted score values only
        data <- lapply(data.list, function(x) x["not.weighted", ])

        # convert list to matrix
        data <- do.call("rbind", data)
        rownames(data) <- c()

        model.scores <- data
    }

    # compute mean ensemble score, if multiple model runs are considered
    variable.name <- eachModel[[1]][1]
    ref.id <- eachModel[[1]][2]

    myFun.getValues <- function(x) {
        scores <- x[3:ncol(x)]
        return(scores)
    }
    scores <- lapply(eachModel, myFun.getValues)
    mean.ensemble.score <- Reduce(`+`, scores)/length(scores)
    model.scores <- data.frame(variable.name, ref.id, mean.ensemble.score)

    #---------------------------------------------------------------------------
    # Subset and sort model.scores in accordance to bench.scores
    #---------------------------------------------------------------------------

    bench.model.scores <- merge(bench.scores, model.scores, by = c("variable.name", "ref.id"))
    bench.model.scores <- bench.model.scores[order(bench.model.scores$myOrder), ]

    bench.scores <- bench.model.scores[4:9]

    model.scores <- bench.model.scores[11:16]

    bms <- bench.model.scores

    my.row.lab.bench <- paste(bms$variable, bms$mod.id, bms$ref.id, sep = "-")
    my.row.lab.model <- paste(bms$variable, model.id, bms$ref.id, sep = "-")

    #---------------------------------------------------------------------------
    # comparison
    #---------------------------------------------------------------------------

    # convert to raster

    # bench means
    data <- bench.scores
    data <- as.matrix(data)
    data <- raster::raster(data)
    bench <- data

    # model data
    data <- model.scores
    data <- as.matrix(data)
    data <- raster::raster(data)
    model <- data

    # delta
    delta <- model - bench

    # prepare inputs for plot
    my.col.lab <- latex2exp::TeX(c("$S_{bias}$", "$S_{rmse}$", "$S_{phase}$", "$S_{iav}$", "$S_{dist}$", "$S_{overall}$"))

    # raster extent
    my.width <- 2
    data <- delta
    my.extent <- c(0, my.width * ncol(data), 0.5, nrow(data) + 0.5)
    raster::extent(bench) <- my.extent
    raster::extent(model) <- my.extent
    raster::extent(delta) <- my.extent

    # plot
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))
    if (outputDir != FALSE) {
        grDevices::pdf("scores_compare_benchmarks.pdf", width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(1, 3), font.main = 1, oma = c(0, 9, 0, 9), mar = c(13, 0, 3, 0), lwd = 1, cex = 1, family = "Helvetica")

    # (a) plot benchmark scores

    # legend
    min <- 0
    max <- 1
    interval <- 0.1
    my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
    my.labels <- round(seq(min, max, interval), 3)  # location of labels
    my.col <- viridis::plasma(n = length(my.breaks) - 1, direction = -1)

    legend.bar.text <- "Benchmark Scores (-)"
    my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
    my.legend.args <- list(text = legend.bar.text, side = 3, font = 1, line = 1, cex = 0.75)

    data <- bench
    my.row.lab <- my.row.lab.bench
    raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = NA, axes = FALSE)
    intFun.addBWtext(data, myDigits = 2)
    axis(side = 2, at = rev(seq(1, nrow(data), 1)), labels = my.row.lab, las = 2)
    axis(side = 1, at = seq(1, my.width * ncol(data), my.width), labels = my.col.lab, las = 2)
    graphics::title("(a)", line = 0.5, adj = 0, font.main = 2, cex.main = 1.25)
    raster::plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args,
        legend.width = 1, legend.shrink = 1, font = 1, horizontal = TRUE)

    # (b) plot delta, i.e. model - benchmark

    data <- delta
    # legend
    mmi.delta <- intFun.min.max.int(data)
    min <- mmi.delta[1]
    max <- mmi.delta[2]
    interval <- mmi.delta[3]
    my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
    my.labels <- round(seq(min, max, interval), 3)  # location of labels
    my.col <- scico::scico(n = length(my.breaks) - 1, palette = "vik")
    # my.col <- viridis::viridis(n = length(my.breaks) - 1, direction = -1)
    legend.bar.text <- paste(model.id, "minus Benchmark Scores (-)", sep = " ")
    my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
    my.legend.args <- list(text = legend.bar.text, side = 3, font = 1, line = 1, cex = 0.75)

    # For additional analysis, sort score difference table to see best/worst performance df <- data.frame(as.matrix(data))
    # colnames(df) <- c('S_bias', 'S_rmse', 'S_phase', 'S_iav', 'S_dist', 'S_overall') rownames(df) <- my.row.lab variableName <-
    # strsplit(my.row.lab, '-') variableName <- t(as.data.frame(variableName)) row.names(variableName) <- c() variableName <-
    # variableName[,1] df <- data.frame(variableName, df) drop 'NEE−FluxCom−CT2019' as NEE-FluxCom wrong in tropics df <-
    # subset(df, rownames(df) != 'NEE-FluxCom-CT2019') compute mean for each variable df <- tapply(df$S_overall, INDEX =
    # list(df$variableName), FUN = mean) sort df <- sort(df) round(df, 2)

    # plot
    raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = NA, axes = FALSE)

    intFun.addBWtext(delta, myDigits = 2)

    # axis(side=2, at=rev(seq(1,nrow(data),1)), labels = my.row.lab, las=2)
    axis(side = 1, at = seq(1, my.width * ncol(data), my.width), labels = my.col.lab, las = 2)
    graphics::title("(b)", line = 0.5, adj = 0, font.main = 2, cex.main = 1.25)
    raster::plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args,
        legend.width = 1, legend.shrink = 1, font = 1, horizontal = TRUE)

    # (c) plot model legend
    min <- 0
    max <- 1
    interval <- 0.1
    my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
    my.labels <- round(seq(min, max, interval), 3)  # location of labels
    my.col <- viridis::plasma(n = length(my.breaks) - 1, direction = -1)
    my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
    my.legend.args <- list(text = legend.bar.text, side = 3, font = 1, line = 1, cex = 0.75)

    data <- model
    my.row.lab <- my.row.lab.model
    legend.bar.text <- paste(model.id, "Scores (-)", sep = " ")
    my.legend.args <- list(text = legend.bar.text, side = 3, font = 1, line = 1, cex = 0.75)

    raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = NA, axes = FALSE)
    intFun.addBWtext(data, myDigits = 2)
    axis(side = 4, at = rev(seq(1, nrow(data), 1)), labels = my.row.lab, las = 2)
    axis(side = 1, at = seq(1, my.width * ncol(data), my.width), labels = my.col.lab, las = 2)
    graphics::title("(c)", line = 0.5, adj = 0, font.main = 2, cex.main = 1.25)
    raster::plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args,
        legend.width = 1, legend.shrink = 1, font = 1, horizontal = TRUE)


    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}

