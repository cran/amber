################################################################################ 
#' Summarize results in a table and a plot
#' @description This function merges all tables that have been created by the
#' functions
#' \link{scores.fluxnet.csv} or \link{scores.fluxnet.nc},
#' \link{scores.grid.notime},
#' \link{scores.grid.time}, and
#' \link{scores.runoff}.
#' @param plot.width Number that gives the plot width, e.g. 6
#' @param plot.height Number that gives the plot height, e.g. 5
#' @param myMargin An R object that gives the figure margins, e.g. c(4, 13, 3, 4)
#' @param inputDir A string that gives the input directory, e.g. '/home/project/study'.
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @return (1) Two tables in LaTeX format that shows the scores of all variables that were assessed (with and without mass weighting).
#' (2) Two figures in PDF format that show the same information as (1).
#' (3) Two NetCDF files that show the same information as (1).
#' (4) Five tables in LaTeX format that show the globally averaged statistical
#' metrics for calculating a score (without mass weighting).
#' (\eqn{S_{bias}}, \eqn{S_{rmse}}, \eqn{S_{phase}}, \eqn{S_{iav}}, \eqn{S_{dist}}).
#' @examples
#'
#' library(amber)
#' library(latex2exp)
#' library(viridis)
#' library(xtable)
#'
#' myInputDir <- paste(system.file('extdata', package = 'amber'), 'scores', sep = '/')
#' scores.tables(inputDir = myInputDir)
#'
#' @export
scores.tables <- function(plot.width = 6, plot.height = 5, myMargin = c(4, 13, 3, 4), inputDir = getwd(), outputDir = FALSE) {
    # summarize results in table
    my.list <- list.files(path = inputDir, pattern = "scorevalues_")
    not <- list.files(pattern = "all_scorevalues_")  # exclude these files
    my.list <- setdiff(my.list, not)
    my.files <- paste(inputDir, my.list, sep = "/")
    data <- lapply(my.files, utils::read.table)
    data <- do.call("rbind", data)
    weights <- rep(c("not.weighted", "weighted"), nrow(data)/2)
    data <- data.frame(weights, data)
    row.names(data) <- c()  # drop row names
    colnames(data) <- c("weights", "variable", "reference", "$S_{bias}$", "$S_{rmse}$", "$S_{phase}$", "$S_{iav}$", "$S_{dist}$", 
        "$S_{overall}$")
    data.nw <- subset(data, weights == "not.weighted")
    data.nw <- subset(data.nw, select = -c(weights))
    data.w <- subset(data, weights == "weighted")
    data.w <- subset(data.w, select = -c(weights))
    rownames(data.nw) <- c()  # omit rownames
    rownames(data.w) <- c()  # omit rownames
    
    data.nw.latex <- data.nw  # used for making a table in LaTeX
    data.w.latex <- data.w  # used for making a table in LaTeX
    
    # convert to LaTeX
    data.nw <- xtable::xtable(data.nw.latex)
    data.w <- xtable::xtable(data.w.latex)
    xtable::caption(data.nw) <- "Scores (without mass weighting)"
    xtable::caption(data.w) <- "Scores (with mass weighting)"
    if (outputDir != FALSE) {
        xtable::print.xtable(data.nw, include.rownames = FALSE, label = "tab:scores_nw", type = "latex", file = "score_without_massweighting.tex", 
            caption.placement = "top", sanitize.text.function = function(x) {
                x
            })
    }
    if (outputDir != FALSE) {
        xtable::print.xtable(data.w, include.rownames = FALSE, label = "tab:scores_w", type = "latex", file = "score_with_massweighting.tex", 
            caption.placement = "top", sanitize.text.function = function(x) {
                x
            })
    }
    
    # summary table with globally averaged inputs for computing scores
    
    my.list <- list.files(path = inputDir, pattern = "scoreinputs_")
    my.files <- paste(inputDir, my.list, sep = "/")
    data <- lapply(my.files, utils::read.table)
    data <- do.call("rbind", data)
    colnames(data)
    colnames(data) <- c("name", "variable", "reference", "unit", "$v_{mod}$", "$v_{ref}$", "$bias$", "$\\sigma_{ref}$", 
        "$\\epsilon_{bias}$ (-)", "$S_{bias}$ (-)", "$rmse$", "$crmse$", "$\\sigma_{ref}$", "$\\epsilon_{rmse}$ (-)", "$S_{rmse}$ (-)", 
        "$max_{cmod}$", "$max_{cref}$", "$\\theta$ (months)", "$S_{phase}$ (-)", "$iav_{mod}$", "$iav_{ref}$", "$\\epsilon_{iav}$ (-)", 
        "$S_{iav}$", "$\\sigma_{\\overline{v_{mod}}}$", "$\\sigma_{\\overline{v_{ref}}}$", "$\\sigma$ (-)", "$R$ (-)", "$S_{dist}$ (-)")
    rownames(data) <- c()  # omit rownames
    
    info <- data[1:4]
    scoreinputs.Sbias <- cbind(info, data[5:10])
    scoreinputs.Srmse <- cbind(info, data[11:15])
    scoreinputs.Sphase <- cbind(info, data[16:19])
    scoreinputs.Siav <- cbind(info, data[20:23])
    scoreinputs.Sdist <- cbind(info, data[24:28])
    
    rownames(scoreinputs.Sbias) <- c()  # omit rownames
    rownames(scoreinputs.Srmse) <- c()  # omit rownames
    rownames(scoreinputs.Sphase) <- c()  # omit rownames
    rownames(scoreinputs.Siav) <- c()  # omit rownames
    rownames(scoreinputs.Sdist) <- c()  # omit rownames
    
    # convert to LaTeX
    scoreinputs.Sbias <- xtable::xtable(scoreinputs.Sbias)
    scoreinputs.Srmse <- xtable::xtable(scoreinputs.Srmse)
    scoreinputs.Sphase <- xtable::xtable(scoreinputs.Sphase)
    scoreinputs.Siav <- xtable::xtable(scoreinputs.Siav)
    scoreinputs.Sdist <- xtable::xtable(scoreinputs.Sdist)
    
    xtable::caption(scoreinputs.Sbias) <- "Globally averaged statistical metrics for calculating $S_{bias}$ (without mass weighting)"
    xtable::caption(scoreinputs.Srmse) <- "Globally averaged statistical metrics for calculating $S_{rmse}$ (without mass weighting)"
    xtable::caption(scoreinputs.Sphase) <- "Globally averaged statistical metrics for calculating $S_{phase}$ (without mass weighting)"
    xtable::caption(scoreinputs.Siav) <- "Globally averaged statistical metrics for calculating $S_{iav}$ (without mass weighting)"
    xtable::caption(scoreinputs.Sdist) <- "Globally averaged statistical metrics for calculating $S_{dist}$ (without mass weighting)"
    
    if (outputDir != FALSE) {
        xtable::print.xtable(scoreinputs.Sbias, include.rownames = FALSE, label = "tab:global_stats", type = "latex", file = "scoreinputs.Sbias.tex", 
            caption.placement = "top", sanitize.text.function = function(x) {
                x
            })
    }
    if (outputDir != FALSE) {
        xtable::print.xtable(scoreinputs.Srmse, include.rownames = FALSE, label = "tab:global_stats", type = "latex", file = "scoreinputs.Srmse.tex", 
            caption.placement = "top", sanitize.text.function = function(x) {
                x
            })
    }
    if (outputDir != FALSE) {
        xtable::print.xtable(scoreinputs.Sphase, include.rownames = FALSE, label = "tab:global_stats", type = "latex", file = "scoreinputs.Sphase.tex", 
            caption.placement = "top", sanitize.text.function = function(x) {
                x
            })
    }
    if (outputDir != FALSE) {
        xtable::print.xtable(scoreinputs.Siav, include.rownames = FALSE, label = "tab:global_stats", type = "latex", file = "scoreinputs.Siav.tex", 
            caption.placement = "top", sanitize.text.function = function(x) {
                x
            })
    }
    if (outputDir != FALSE) {
        xtable::print.xtable(scoreinputs.Sdist, include.rownames = FALSE, label = "tab:global_stats", type = "latex", file = "scoreinputs.Sdist.tex", 
            caption.placement = "top", sanitize.text.function = function(x) {
                x
            })
    }
    
    # Score maps make a score map without mass weighting get row and column labels
    data <- data.nw
    my.col.lab <- colnames(data[3:ncol(data)])
    my.col.lab <- latex2exp::TeX(my.col.lab)
    my.row.lab <- paste(data$variable, data$reference, sep = "-")
    # convert to raster
    data <- data[3:ncol(data)]
    colnames(data) <- c()
    data <- as.matrix(data)
    data <- raster::raster(data)
    my.width <- 2
    raster::extent(data) <- c(0, my.width * ncol(data), 0.5, nrow(data) + 0.5)
    # colors
    min <- 0
    max <- 1
    interval <- 0.1
    my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
    my.labels <- round(seq(min, max, interval), 3)  # location of labels
    my.col <- viridis::plasma(n = length(my.breaks) - 1, direction = -1)
    # misc
    legend.bar.text <- "score (-)"
    my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
    my.legend.args <- list(text = legend.bar.text, side = 2, font = 1, line = 1, cex = 1)
    my.filename <- "score_without_massweighting"
    my.title <- "Scores"  # (no mass weighting)
    
    # save as netcdf file
    
    if (outputDir != FALSE) {
        raster::writeRaster(data, filename = paste(outputDir, "/", my.filename, ".nc", sep = ""), format = "CDF", varname = "scores", 
            overwrite = TRUE)
    }
    
    # plot
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar))
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, ".pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(1, 1), font.main = 1, mar = myMargin, lwd = 1, cex = 1)
    raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = my.title, axes = FALSE)
    raster::text(data, digits = 2, cex = 0.7)
    graphics::axis(side = 2, at = rev(seq(1, nrow(data), 1)), labels = my.row.lab, las = 2)
    graphics::axis(side = 1, at = seq(1, my.width * ncol(data), my.width), labels = my.col.lab, las = 2)
    
    raster::plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args, 
        legend.width = 1.5, legend.shrink = 1, font = 1)
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
    
    # make a score map with mass weighting
    
    # get row and column labels
    data <- data.w
    my.col.lab <- colnames(data[3:ncol(data)])
    my.col.lab <- latex2exp::TeX(my.col.lab)
    my.row.lab <- paste(data$variable, data$reference, sep = "-")
    # convert to raster
    data <- data[3:ncol(data)]
    colnames(data) <- c()
    data <- as.matrix(data)
    data <- raster::raster(data)
    my.width <- 2
    raster::extent(data) <- c(0, my.width * ncol(data), 0.5, nrow(data) + 0.5)
    # colors
    min <- 0
    max <- 1
    interval <- 0.1
    my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
    my.labels <- round(seq(min, max, interval), 3)  # location of labels
    my.col <- viridis::plasma(n = length(my.breaks) - 1, direction = -1)
    # misc
    legend.bar.text <- "score (-)"
    my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
    my.legend.args <- list(text = legend.bar.text, side = 2, font = 1, line = 1, cex = 1)
    my.filename <- "score_with_massweighting"
    my.title <- "Scores (with mass weighting)"
    
    # save as netcdf file
    if (outputDir != FALSE) {
        raster::writeRaster(data, filename = paste(outputDir, "/", my.filename, ".nc", sep = ""), format = "CDF", varname = "scores", 
            overwrite = TRUE)
    }
    
    # plot
    if (outputDir != FALSE) {
        grDevices::pdf(paste(outputDir, "/", my.filename, ".pdf", sep = ""), width = plot.width, height = plot.height)
    }
    graphics::par(mfrow = c(1, 1), font.main = 1, mar = myMargin, lwd = 1, cex = 1)
    raster::plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = my.title, axes = FALSE)
    raster::text(data, digits = 2, cex = 0.7)
    graphics::axis(side = 2, at = rev(seq(1, nrow(data), 1)), labels = my.row.lab, las = 2)
    graphics::axis(side = 1, at = seq(1, my.width * ncol(data), my.width), labels = my.col.lab, las = 2)
    
    raster::plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args, 
        legend.width = 1.5, legend.shrink = 1, font = 1)
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}

