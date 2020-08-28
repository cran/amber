################################################################################
#' Barplots of reference and model data with same unit
#'
#' @description This function produces barplots of annual mean values for a selection of variables and reference data.
#' The user may choose to plot both, absolute and fractional values. A typical application
#' is a barplot for variables related to the surface energy balance.
#' @param mod.path.list A list with paths for each model run, e.g. mod.path.list <- list(mod01.path, mod02.path, mod03.path).
#' The respective folders must each contain the file 'metricsTable' produced by the function \link{scores.tables}.
#' @param modelIDs An R object with the different model run IDs, e.g. c('CLASSIC.CRUJRAv2', 'CLASSIC.GSWP3W5E5', 'CLASSIC.CRUNCEP')
#' @param variableNames Names of variables that should be plotted, e.g.  c('HFG', 'HFLS', 'HFSS', 'RNS')
#' @param referenceNames Names of reference data that should be plotted, e.g. c('FLUXNET', 'CLASS')
#' @param showFractions Logical. If set to TRUE, the Figure will show fractional values,
#' such as the fractions of net surface radiation that are converted into heat fluxes.
#' @param targetVariable A string that gives the variable for which fractions are
#' computed, e.g. 'RNS' for net surface radiation. Only relevant if showFractions = TRUE.
#' @param ofileName A string that gives the output file name, e.g. 'barplot.pdf'
#' @param plot.width Number that gives the plot width, e.g. 12
#' @param plot.height Number that gives the plot height, e.g. 8
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @param subcaption Logical. If TRUE then the subcaptions (a) and (b) are added.
#' @return Barplots of annual mean values annual mean values for a selection of variables and reference data.
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
#'
#' mod.path.list <- list(mod01.path, mod02.path)
#' modelIDs <- c('CLASSIC.CRUJRAv2', 'CLASSIC.GSWP3W5E5')
#' variableNames <- c('HFG', 'HFLS', 'HFSS', 'RNS')
#' referenceNames <- c('FLUXNET', 'CLASS')
#'
#' plotBars(mod.path.list, modelIDs, variableNames, referenceNames,
#' showFractions = TRUE, targetVariable = 'RNS', ofileName = 'barplot.pdf',
#' plot.width = 5, plot.height = 5, outputDir = FALSE)
#'
#' @export
plotBars <- function(mod.path.list, modelIDs, variableNames, referenceNames, showFractions = FALSE, targetVariable = "RNS", ofileName = "barplot.pdf",
    plot.width = 12, plot.height = 8, outputDir = FALSE, subcaption = FALSE) {

    # extract data for selected variables and reference data
    id <- expand.grid(variableNames, referenceNames)
    id <- data.frame(id)
    colnames(id) <- c("variable.name", "ref.id")

    nmod <- length(mod.path.list)  # number of model runs

    # Make a nested list myList[[i]][1] give reference data for model i myList[[i]][2] give model data for model i

    myList <- foreach::foreach(i = 1:nmod) %do% {
        data <- utils::read.table(file.path(mod.path.list[[i]], "metricsTable"))
        data <- subset(data, variable.name %in% id$variable.name)
        data <- subset(data, ref.id %in% id$ref.id)
        data <- merge(data, id, by = c("variable.name", "ref.id"), all = TRUE)

        variable <- data$variable.name
        unit <- stats::na.omit(data)$variable.unit[1]
        ref.id <- data$ref.id
        ref <- data$ref.mean.scalar
        mod <- data$mod.mean.scalar

        data <- data.frame(variable, unit, ref.id, ref, mod)
        unit <- data$unit[1]
        data <- subset(data, select = -c(unit), drop = TRUE)
        dataList <- split(data, ref.id, drop = TRUE)
        variableNames <- dataList[[1]]$variable

        fun.getValues <- function(x) {
            x <- subset(x, select = -c(variable, ref.id), drop = FALSE)
            return(x)
        }

        dataList.values <- lapply(dataList, FUN = fun.getValues)
        values <- do.call("cbind", dataList.values)

        # subset model
        values.mod <- subset(values, select = grepl(".mod", names(values)))
        values.mod <- as.matrix(values.mod)
        colnames(values.mod) <- gsub(".mod", "", colnames(values.mod))
        values.mod <- as.matrix(values.mod)

        # subset ref
        values.ref <- subset(values, select = grepl(".ref", names(values)))
        values.ref <- as.matrix(values.ref)
        colnames(values.ref) <- gsub(".ref", "", colnames(values.ref))
        values.ref <- as.matrix(values.ref)

        list(values.ref, values.mod, variableNames)
    }

    variableNames <- as.character(myList[[1]][[3]])

    # get min and max values to define ylim
    minAndMax <- foreach::foreach(i = 1:nmod) %do% {
        values <- unlist(myList[[i]][1:2])
        min <- min(values, na.rm = TRUE)
        max <- max(values, na.rm = TRUE)
        minAndMax <- c(min, max)
    }
    minAndMax <- do.call("cbind", minAndMax)
    min <- min(minAndMax)
    max <- max(minAndMax)

    min <- min - 0.1 * abs(min)
    if (min > 0) {
        min <- 0
    }
    max <- 1.1 * max
    my.ylim <- c(min, max)

    # Plot
    myColRef <- viridis::viridis(n = nrow(myList[[1]][[1]]), end = 0.75, option = "viridis")
    # myColRef <- 'grey'
    myColMod <- viridis::viridis(n = nmod, begin = 0.15, end = 0.85, option = "plasma")

    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", ofileName, sep = ""), width = plot.width, height = plot.height)
    }

    graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(5, 5, 2, 0.5), lwd = 1, cex = 1, tck = 0.03, las = 1, xpd = TRUE)

    values.ref <- myList[[1]][[1]]  #
    x <- graphics::barplot(values.ref, ylim = my.ylim, col = myColRef, ylab = latex2exp::TeX(unit), legend.text = FALSE, las = 1,
        border = NA, beside = TRUE)

    for (i in 1:nmod) {
        graphics::barplot(myList[[i]][[2]], col = NA, border = myColMod[i], beside = TRUE, add = TRUE, axes = FALSE, axisnames = FALSE)
    }
    graphics::barplot(values.ref, col = NA, border = "black", beside = TRUE, add = TRUE, axes = FALSE, axisnames = FALSE)
    graphics::legend("topleft", modelIDs, col = myColMod[1:nmod], pch = 16, bty = "n", cex = 0.5)
    graphics::legend("bottom", variableNames, col = myColRef, pch = 16, bty = "n", horiz = TRUE, inset = -0.3)
    if (subcaption == TRUE) {
        graphics::title("(a)", line = 0.5, adj = 0, font.main = 2, cex.main = 1.25)
    }
    graphics::box()

    if (outputDir != FALSE) {
        grDevices::dev.off()
    }

    #---------------------------------------------------------------------------

    # relative partitioning

    #---------------------------------------------------------------------------
    if (showFractions == TRUE) {

        myList <- foreach::foreach(i = 1:nmod) %do% {
            data <- utils::read.table(file.path(mod.path.list[[i]], "metricsTable"))
            data <- subset(data, variable.name %in% id$variable.name)
            data <- subset(data, ref.id %in% id$ref.id)
            data <- merge(data, id, by = c("variable.name", "ref.id"), all = TRUE)

            variable <- data$variable.name
            unit <- data$variable.unit
            ref.id <- data$ref.id
            ref <- data$ref.mean.scalar
            mod <- data$mod.mean.scalar

            data <- data.frame(variable, unit, ref.id, ref, mod)
            unit <- data$unit[1]
            data <- subset(data, select = -c(unit), drop = TRUE)
            dataList <- split(data, ref.id, drop = TRUE)
            variableNames <- dataList[[1]]$variable

            fun.getFractions <- function(x) {
                total <- subset(x, variable == targetVariable)
                total <- subset(total, select = -c(variable, ref.id), drop = FALSE)
                total <- as.matrix(total)
                x <- subset(x, select = -c(variable, ref.id), drop = FALSE)
                x.frac <- sweep(x, 2, total, FUN = "/")
                return(x.frac)
            }


            dataList.fractions <- lapply(dataList, FUN = fun.getFractions)
            fractions <- do.call("cbind", dataList.fractions)
            fractions <- as.matrix(round(fractions, 2))

            fractions.mod <- subset(fractions, select = grepl(".mod", colnames(fractions)))
            fractions.ref <- subset(fractions, select = grepl(".ref", colnames(fractions)))

            colnames(fractions.mod) <- gsub(".mod", "", colnames(fractions.mod))
            colnames(fractions.ref) <- gsub(".ref", "", colnames(fractions.ref))

            list(fractions.ref, fractions.mod, variableNames)
        }

        # Plot
        myColRef <- viridis::viridis(n = nrow(myList[[1]][[1]]), end = 0.75, option = "viridis")
        myColMod <- viridis::viridis(n = nmod, begin = 0.15, end = 0.85, option = "plasma")

        if (outputDir != FALSE) {
            grDevices::pdf(paste(outputDir, "/", "Fraction", ofileName, sep = ""), width = plot.width, height = plot.height)
        }

        graphics::par(mfrow = c(1, 1), font.main = 1, mar = c(5, 5, 2, 0.5), lwd = 1, cex = 1, tck = 0.03, las = 1, xpd = TRUE)

        values.ref <- myList[[1]][[1]]  #
        x <- graphics::barplot(values.ref, ylim = c(0, 1.1), col = myColRef, ylab = "Fraction (-)", legend.text = FALSE, las = 1,
            border = NA, beside = TRUE)

        for (i in 1:nmod) {
            graphics::barplot(myList[[i]][[2]], col = NA, border = myColMod[i], beside = TRUE, add = TRUE, axes = FALSE, axisnames = FALSE)
        }
        graphics::barplot(values.ref, col = NA, border = "black", beside = TRUE, add = TRUE, axes = FALSE, axisnames = FALSE)
        # graphics::legend('topleft', modelIDs, col = myColMod[1:nmod], pch = 16, bty = 'n')
        graphics::legend("bottom", as.character(variableNames), col = myColRef, pch = 16, bty = "n", horiz = TRUE, inset = -0.3)
        if (subcaption == TRUE) {
            graphics::title("(b)", line = 0.5, adj = 0, font.main = 2, cex.main = 1.25)
        }
        graphics::box()

        if (outputDir != FALSE) {
            grDevices::dev.off()
        }
    }

}
utils::globalVariables(c("%do%", "variable.name", "zonal.mean.ref"))
