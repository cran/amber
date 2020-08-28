################################################################################
#' Plots that show statistical metrics for FLUXNET sites
#' @description This function plots statistical metrics for the comparison against FLUXNET data.
#' @param inputDir A string that gives the location of text files produced by \link{scores.fluxnet.csv}, e.g. '/home/project/study'.
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @param variableNames A string of six variables that should be plotted. Default is set to c('GPP', 'RECO', 'NEE', 'RNS', 'HFLS', 'HFSS').
#' @param mod.id A string that identifies the source of the reference data set, e.g. 'CLASSIC'
#' @param plot.width Number that gives the plot width, e.g. 12
#' @param plot.height Number that gives the plot height, e.g. 8
#' @return Figures that show statistical metrics produced by the function \link{scores.fluxnet.csv}. The function expects six input variables.
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
#' inputDir <- paste(system.file('extdata', package = 'amber'), 'scores', sep = '/')
#' plotFluxnetStats(inputDir, outputDir = FALSE, mod.id = 'CLASSIC.CRUJRAv2')
#'
#' @export
plotFluxnetStats <- function(inputDir, outputDir = FALSE, mod.id = "CLASSIC", variableNames = c("GPP", "RECO", "NEE", "RNS",
    "HFLS", "HFSS"), plot.width = 12, plot.height = 8) {

    for (i in 1:length(variableNames)) {
        # data
        dataFileName <- paste(variableNames[i], "FLUXNET", sep = "_")
        data <- utils::read.table(file.path(inputDir, dataFileName))
        assign(variableNames[i], data)

        # unit
        unitFilePattern <- paste("scoreinputs", variableNames[i], mod.id, "vs_FLUXNET", sep = "_")
        unitFileName <- list.files(path = inputDir, pattern = unitFilePattern)
        unitFile <- utils::read.table(file.path(inputDir, unitFileName))
        unit <- unitFile$variable.unit  # in LaTeX format
        unit <- latex2exp::TeX(unit[[1]])
        assign(paste(variableNames[i], "Unit", sep = ""), unit)
    }

    #---------------------------------------------------------------------------
    # plot
    #---------------------------------------------------------------------------
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))
    #---------------------------------------------------------------------------
    # scatter plots for all variables (mean)
    #---------------------------------------------------------------------------
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", "FLUXNET_scatterplot_mean.pdf", sep = ""), width = plot.width, height = plot.height)
    }

    graphics::par(mfrow = c(2, 3), font.main = 1, mar = c(4, 4, 1, 1), lwd = 1, cex = 1, tck = 0.03, las = 1)

    for (i in 1:length(variableNames)) {

        # get data and unit
        data <- get(variableNames[i])
        unit <- get(paste(variableNames[i], "Unit", sep = ""))

        # plot characteristics
        minLim <- min(data$ref.mean, data$mod.mean, na.rm = TRUE)
        maxLim <- max(data$ref.mean, data$mod.mean, na.rm = TRUE)
        myLim <- c(minLim, maxLim)
        my.xlab <- paste(variableNames[i], "(FLUXNET)", sep = " ")
        my.ylab <- paste(variableNames[i], " (", mod.id, ")", sep = "")
        subplot <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")

        plot(data$ref.mean, data$mod.mean, xlim = myLim, ylim = myLim, xlab = my.xlab, ylab = my.ylab)
        graphics::abline(0, 1)
        graphics::legend("topleft", subplot[i])
        graphics::legend("bottomright", unit, bty = "n")
        graphics::axis(3, labels = FALSE, tcl = 0.3)
        graphics::axis(4, labels = FALSE, tcl = 0.3)
    }
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }

    #---------------------------------------------------------------------------
    # scatter plots for all variables (iav)
    #---------------------------------------------------------------------------
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", "FLUXNET_scatterplot_iav.pdf", sep = ""), width = plot.width, height = plot.height)
    }

    graphics::par(mfrow = c(2, 3), font.main = 1, mar = c(4, 4, 1, 1), lwd = 1, cex = 1, tck = 0.03, las = 1)

    for (i in 1:length(variableNames)) {

        # get data and unit
        data <- get(variableNames[i])
        unit <- get(paste(variableNames[i], "Unit", sep = ""))

        # plot characteristics
        minLim <- min(data$ref.iav, data$mod.iav, na.rm = TRUE)
        maxLim <- max(data$ref.iav, data$mod.iav, na.rm = TRUE)
        myLim <- c(minLim, maxLim)
        my.xlab <- paste(variableNames[i], "(FLUXNET)", sep = " ")
        my.ylab <- paste(variableNames[i], " (", mod.id, ")", sep = "")
        subplot <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")

        plot(data$ref.iav, data$mod.iav, xlim = myLim, ylim = myLim, xlab = my.xlab, ylab = my.ylab)
        graphics::abline(0, 1)
        graphics::legend("topleft", subplot[i])
        graphics::legend("bottomright", unit, bty = "n")
        graphics::axis(3, labels = FALSE, tcl = 0.3)
        graphics::axis(4, labels = FALSE, tcl = 0.3)
    }
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }

    #---------------------------------------------------------------------------
    # histogram (bias)
    #---------------------------------------------------------------------------
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", "FLUXNET_histogram_bias.pdf", sep = ""), width = plot.width, height = plot.height)
    }

    graphics::par(mfrow = c(2, 3), font.main = 1, mar = c(4, 4, 1, 1), lwd = 1, cex = 1, tck = -0.03, las = 1)

    for (i in 1:length(variableNames)) {

        # get data and unit
        data <- get(variableNames[i])
        unit <- get(paste(variableNames[i], "Unit", sep = ""))

        # plot characteristics
        my.xlab <- paste(variableNames[i], "Bias (Model minus FLUXNET)", sep = " ")
        my.ylab <- "Number of Sites"
        subplot <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")

        graphics::hist(data$bias, col = "grey", border = "grey", xlab = my.xlab, ylab = my.ylab, main = NA)
        graphics::abline(v = 0)
        graphics::legend("topleft", subplot[i])
        graphics::legend("topright", unit, bty = "n")
        graphics::box()
    }
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }

    #---------------------------------------------------------------------------
    # histogram (crmse)
    #---------------------------------------------------------------------------
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", "FLUXNET_histogram_crmse.pdf", sep = ""), width = plot.width, height = plot.height)
    }

    graphics::par(mfrow = c(2, 3), font.main = 1, mar = c(4, 4, 1, 1), lwd = 1, cex = 1, tck = -0.03, las = 1)

    for (i in 1:length(variableNames)) {

        # get data and unit
        data <- get(variableNames[i])
        unit <- get(paste(variableNames[i], "Unit", sep = ""))

        # plot characteristics
        my.xlab <- paste(variableNames[i], "CRMSE", sep = " ")
        my.ylab <- "Number of Sites"
        subplot <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")

        graphics::hist(data$crmse, col = "grey", border = "grey", xlab = my.xlab, ylab = my.ylab, main = NA)
        graphics::legend("topleft", subplot[i])
        graphics::legend("topright", unit, bty = "n")
        graphics::box()
    }
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }

    #---------------------------------------------------------------------------
    # histogram (phase)
    #---------------------------------------------------------------------------
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", "FLUXNET_histogram_phase.pdf", sep = ""), width = plot.width, height = plot.height)
    }

    graphics::par(mfrow = c(2, 3), font.main = 1, mar = c(4, 4, 1, 1), lwd = 1, cex = 1, tck = -0.03, las = 1)

    for (i in 1:length(variableNames)) {

        # get data and unit
        data <- get(variableNames[i])
        unit <- get(paste(variableNames[i], "Unit", sep = ""))

        # plot characteristics
        my.xlab <- paste(variableNames[i], "Phase bias (Model minus FLUXNET)", sep = " ")
        my.ylab <- "Number of Sites"
        subplot <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")

        graphics::hist(data$phase, breaks = seq(0, 6, 1), col = "grey", border = "grey", xlab = my.xlab, ylab = my.ylab, main = NA)
        graphics::legend("topleft", subplot[i])
        graphics::legend("topright", unit, bty = "n")
        graphics::box()
    }
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}

