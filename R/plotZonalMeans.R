################################################################################
#' Plot zonal mean values of model and reference data in a single figure
#' @description This function plots zonal mean values of model and reference data, computed from \link{zonalMeanStats}, for a selection of variables.
#' @param mod.path.list A List of directories where AMBER output is stored for different model runs,
#' e.g. list(mod01.path, mod02.path, mod03.path)
#' @param modelIDs An R object with the different model run IDs, e.g. c('CLASSIC.CRUJRAv2', 'CLASSIC.GSWP3W5E5', 'CLASSIC.CRUNCEP')
#' @param myVariables  An R object with the variable names of interest, e.g. c('GPP.FluxCom', 'RECO.FluxCom').
#' @param metric Specify for what statistical metric you wish to plot zonal means.
#' Current options are 'mean', 'bias', 'crmse', 'phase', 'iav', 'sd', 'bias.score', 'rmse.score', 'phase.score', and 'iav.score'.
#' @param lat.range Latitudinal range of ticks and labels on the horizontal axis, e.g. c(-50, 80).
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @param plot.width A number that gives the desired plot width, e.g. 6.
#' @param plot.height A number that gives the desired plot height, e.g. 10.
#' @param myFilename A string that gives name of the Figure, e.g. 'zonalMeans.pdf'.
#' @param legendLocation A string that specifies the location of the model ID legend, e.g. 'topright'.
#' @param plotRef Logical. If FALSE, reference data is not plotted. Default is TRUE.
#' @param subcaption A string that defines the subcaption of the figure, e.g. '(a)'.
#' @return A figure that shows zonal mean values of model and reference data.
#'
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
#' mod01.path <- system.file('extdata/SIMod01', package = 'amber')
#' mod02.path <- system.file('extdata/SIMod02', package = 'amber')
#' mod.path.list <- list(mod01.path, mod02.path)
#' modelIDs <- c('CLASSIC.CRUJRAv2', 'CLASSIC.GSWP3W5E5')
#' myVariables <- c('GPP.MODIS', 'LAI.AVHRR')
#' plotZonalMeans(mod.path.list, modelIDs, myVariables, metric = 'mean')
#'
#' @export
plotZonalMeans <- function(mod.path.list, modelIDs, myVariables, metric = "mean", lat.range = c(-50, 80), plot.width = 6, plot.height = 10,
    outputDir = FALSE, myFilename = "zonalMeans.pdf", legendLocation = "topright", plotRef = TRUE, subcaption = "") {

    myPalette <- "plasma"

    nmod <- length(mod.path.list)  # number of model runs
    myList <- foreach::foreach(i = 1:nmod) %do% {

        mod.id <- modelIDs[i]

        units <- utils::read.table(file.path(mod.path.list[[i]], "zonalMeanStatsUnits"))
        data <- utils::read.table(file.path(mod.path.list[[i]], "zonalMeanStats"))

        zone <- data$zone
        names <- colnames(data)
        names <- names[-1]  # drop column zone

        variables <- data.frame(myVariables)
        colnames(variables) <- "variables"
        id <- seq(1, nrow(variables), 1)
        variables <- data.frame(id, variables)

        # merge variable names and units
        variables <- merge(variables, units, by = "variables")
        variables <- variables[order(variables$id), ]

        names <- toString(variables[, 1])
        names <- unlist(strsplit(names, ", "))
        units <- variables[, 3]
        units <- gsub("\\", "", units, fixed = TRUE)  # NEW

        list(mod.id, data, names, units)
    }

    units <- myList[[1]][[4]]  # same for all model runs
    names <- myList[[1]][[3]]  # same for all model runs

    # plot
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))

    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", myFilename, sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(length(names), 1), font.main = 1, mar = c(0, 0, 0, 0), oma = c(5, 6, 2, 5), lwd = 1, cex = 1, xpd = NA,
        tck = 0.03, las = 1)

    # loop for each variable name (i), e.g. GPP.MODIS
    for (i in 1:length(names)) {

        # function that checks whether integer is even or odd this will be used for deciding where to add labels
        is.even <- function(x) x%%2 == 0
        # replace dot with dash in filename to avoid problems when using file in LaTeX
        variable.name <- gsub(pattern = ".", replacement = "-", names[i], fixed = TRUE)
        # myTitle <- gsub(pattern = '.', replacement = ' ', names[i], fixed = TRUE) mean my.ylab <- 'mean'
        my.unit <- latex2exp::TeX(units[i])
        if (metric == "phase") {
            my.unit <- "months"
        }

        if (metric == "score") {
            my.unit <- "(-)"
        }

        # loop for each model run (j)
        myListMetric <- foreach::foreach(j = 1:nmod) %do% {
            data <- myList[[j]][[2]]
            names <- myList[[j]][[3]]
            if (metric == "mean" | metric == "bias") {
                colname.mod <- paste(names[i], ".mod.mean", sep = "")
                colname.ref <- paste(names[i], ".ref.mean", sep = "")
            }

            if (metric == "phase") {
                colname.mod <- paste(names[i], ".mod.max.month", sep = "")
                colname.ref <- paste(names[i], ".ref.max.month", sep = "")
            }

            if (metric == "iav") {
                colname.mod <- paste(names[i], ".mod.iav", sep = "")
                colname.ref <- paste(names[i], ".ref.iav", sep = "")
            }

            if (metric == "sd") {
                colname.mod <- paste(names[i], ".mod.sd", sep = "")
                colname.ref <- paste(names[i], ".ref.sd", sep = "")
            }

            if (metric == "crmse") {
                colname.mod <- paste(names[i], ".crmse", sep = "")
                colname.ref <- paste(names[i], ".crmse", sep = "")
            }

            if (metric == "bias.score") {
                colname.mod <- paste(names[i], ".bias.score", sep = "")
                colname.ref <- paste(names[i], ".bias.score", sep = "")
            }

            if (metric == "rmse.score") {
                colname.mod <- paste(names[i], ".rmse.score", sep = "")
                colname.ref <- paste(names[i], ".rmse.score", sep = "")
            }

            if (metric == "phase.score") {
                colname.mod <- paste(names[i], ".phase.score", sep = "")
                colname.ref <- paste(names[i], ".phase.score", sep = "")
            }

            if (metric == "iav.score") {
                colname.mod <- paste(names[i], ".iav.score", sep = "")
                colname.ref <- paste(names[i], ".iav.score", sep = "")
            }

            data.mod <- data[colname.mod]
            data.ref <- data[colname.ref]
            zone <- data[1]
            data.metric <- data.frame(data[1], data.mod[1], data.ref[1])
            if (metric == "bias") {
                bias <- data.mod[1] - data.ref[1]
                data.metric <- data.frame(data[1], bias, bias)
            }
            colnames(data.metric) <- c("zone", "data.mod", "data.ref")
            data.metric
        }

        y.min <- min(unlist(lapply(myListMetric, `[`, 2:3)), na.rm = TRUE)
        y.max <- max(unlist(lapply(myListMetric, `[`, 2:3)), na.rm = TRUE)
        my.ylim <- c(y.min, y.max)
        myCol <- viridis::viridis(n = nmod, begin = 0.15, end = 0.85, option = myPalette)
        myRefCol <- "black"

        # Set color to NA for cases when there are no corresponding reference values:
        if (metric == "bias") {
            myRefCol <- NA
        }
        if (metric == "crmse") {
            myRefCol <- NA
        }
        if (metric == "bias.score") {
            myRefCol <- NA
        }
        if (metric == "rmse.score") {
            myRefCol <- NA
        }
        if (metric == "phase.score") {
            myRefCol <- NA
        }
        if (metric == "iav.score") {
            myRefCol <- NA
        }
        if (plotRef == FALSE) {
            myRefCol <- NA
        }

        plot(myListMetric[[1]]$zone, myListMetric[[1]]$data.ref, col = myRefCol, type = "l", main = NA, ylim = my.ylim, xlab = NA,
            ylab = NA, xaxt = "n", yaxt = "n", panel.first = graphics::rect(-23, y.min, 23, y.max, density = NULL, angle = 0,
                col = "grey90", border = NA))

        # Add zero line for bias plots:

        if (metric == "bias") {
            min.lat <- min(myListMetric[[1]]$zone)
            max.lat <- max(myListMetric[[1]]$zone)
            if (min.lat < 0 && max.lat > 0) {
                graphics::segments(min.lat, 0, max.lat, 0, lty = 3)
            }
        }


        # loop for each model run (j)
        for (j in 1:nmod) {
            graphics::lines(myListMetric[[j]]$zone, myListMetric[[j]]$data.mod, col = myCol[j])
        }
        graphics::legend("topleft", variable.name, bty = "n", pch = 16, col = myRefCol)

        if (i == 1) {
            graphics::legend(legendLocation, modelIDs, col = myCol[1:nmod], pch = 16, bty = "n")
            graphics::axis(3, at = seq(lat.range[1], lat.range[2], 10), labels = FALSE, tcl = 0.3)
            graphics::title(subcaption, line = 0.5, adj = 0, font.main = 2, cex.main = 1.25)
        }

        #-------------------------------------------------------------------------------
        # ticks and labels
        #-------------------------------------------------------------------------------

        graphics::axis(1, at = seq(lat.range[1], lat.range[2], 10), labels = FALSE, tcl = 0.3)

        # for subplots when i is odd:
        if (is.even(i) == FALSE) {
            graphics::axis(2, labels = TRUE, tcl = 0.3)
            graphics::axis(4, labels = FALSE, tcl = 0.3)
            graphics::mtext(my.unit, 2, line = 3.5, las = 3)
        }
        # for subplots when i is even:
        if (is.even(i) == TRUE) {
            graphics::axis(2, labels = FALSE, tcl = 0.3)
            graphics::axis(4, labels = TRUE, tcl = 0.3)
            graphics::mtext(my.unit, 4, line = 3.5, las = 3)
        }

        # for the last subplot
        if (i == max(length(names))) {
            graphics::axis(1, at = seq(lat.range[1], lat.range[2], 10), labels = TRUE, tcl = 0.3)
            graphics::mtext("Degrees Latitude", 1, line = 3, las = 1)
        }
    }
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}

utils::globalVariables("%do%")

