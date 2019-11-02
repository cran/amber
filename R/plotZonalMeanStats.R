################################################################################
#' Plot zonal mean plots of AMBER results (bias, bias scores, etc)
#' @description This function plots results from \link{zonalMeanStats}, i.e. zonal mean values of model and reference data and the zonal mean bias, centralized root-mean-squre error, phase, inter-annual variability, and corresponding scores.
#' @param zonalMeanStats A string that gives the name and location of the zonalMeanStats file produced by \link{zonalMeanStats}, e.g. '/home/project/study/zonalMeanStats'.
#' @param zonalMeanStatsUnits A string that gives the name and location of the zonalMeanStatsUnits file produced by \link{zonalMeanStats}, e.g. '/home/project/study/zonalMeanStatsUnits'.
#' @param lat.range Latitudinal range of ticks and labels on the horizontal axis, e.g. c(-50, 80).
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @return Figures that show zonal mean values of model and reference data, centralized root-mean-squre error, phase, inter-annual variability, and corresponding
#' scores for each variable and globally gridded reference data set.
#'
#' @examples
#'
#' library(amber)
#' library(foreach)
#' library(latex2exp)
#' library(ncdf4)
#' library(raster)
#'
#' zonalMeanStats <- system.file('extdata/zonalMeanStats', 'zonalMeanStats', package = 'amber')
#' zonalMeanStatsUnits <- system.file('extdata/zonalMeanStats',
#'  'zonalMeanStatsUnits', package = 'amber')
#' plotZonalMeanStats(zonalMeanStats, zonalMeanStatsUnits, outputDir = FALSE)
#'
#' @export
plotZonalMeanStats <- function(zonalMeanStats, zonalMeanStatsUnits, lat.range = c(-50, 80), outputDir = FALSE) {

    units <- utils::read.table(zonalMeanStatsUnits)
    data <- utils::read.table(zonalMeanStats)

    zone <- data$zone
    names <- colnames(data)
    names <- names[-1]  # drop column zone

    # select all biases
    bias.score <- grep("bias.score", names)
    bias.score <- names[bias.score]
    variables <- gsub(pattern = ".bias.score", replacement = "", bias.score)
    variables <- data.frame(variables)
    # merge variable names and units
    variables <- merge(variables, units, by = "variables")
    names <- toString(variables[, 1])
    names <- unlist(strsplit(names, ", "))
    units <- variables[, 2]

    # plot
    for (i in 1:length(names)) {
        oldpar <- graphics::par(mfrow = c(1, 2))
        on.exit(graphics::par(oldpar))
        # replace dot with dash in filename to avoid problems when using file in LaTeX
        variable.name <- gsub(pattern = ".", replacement = "-", names[i], fixed = TRUE)
        my.filename <- paste("zonalMeanStats-", variable.name, ".pdf", sep = "")
        plot.width <- 6
        plot.height <- 9
        if (outputDir != FALSE) {
            grDevices::pdf(paste(outputDir, "/", my.filename, sep = ""), width = plot.width, height = plot.height)
        }
        graphics::par(mfrow = c(6, 1), font.main = 1, mar = c(0, 0, 0, 0), oma = c(5, 5, 2, 5), lwd = 1, cex = 1, xpd = NA,
            tck = 0.03, las = 1)

        myTitle <- gsub(pattern = ".", replacement = " ", names[i], fixed = TRUE)
        # mean
        my.ylab <- "mean"
        my.unit <- latex2exp::TeX(units[i])
        colname <- paste(names[i], ".mean", sep = "")
        colname.mod <- paste(names[i], ".mod.mean", sep = "")
        colname.ref <- paste(names[i], ".ref.mean", sep = "")
        mean.mod <- data[colname.mod]
        mean.ref <- data[colname.ref]
        zone <- data[1]
        mean <- data.frame(data[1], mean.mod[1], mean.ref[1])
        colnames(mean) <- c("zone", "mean.mod", "mean.ref")
        y.min <- min(mean$mean.mod, mean$mean.ref, na.rm = TRUE)
        y.max <- max(mean$mean.mod, mean$mean.ref, na.rm = TRUE)
        my.ylim <- c(y.min, y.max)
        plot(mean$zone, mean$mean.mod, type = "l", main = NA, ylim = my.ylim, xlab = NA, ylab = "mean", xaxt = "n", yaxt = "n",
            panel.first = graphics::rect(-23, y.min, 23, y.max, density = NULL, angle = 0, col = "grey90", border = NA))
        graphics::lines(mean$zone, mean$mean.ref, col = "red")
        graphics::legend("topleft", "(a)")
        graphics::legend("topright", c("model", "reference"), col = c("black", "red"), pch = 16, horiz = TRUE, bty = "n")
        graphics::axis(1, at = seq(lat.range[1], lat.range[2], 10), labels = FALSE, tcl = 0.3)
        graphics::axis(2, labels = FALSE, tcl = 0.3)
        graphics::axis(3, at = seq(lat.range[1], lat.range[2], 10), labels = FALSE, tcl = 0.3)
        graphics::axis(4, labels = TRUE, tcl = 0.3)
        graphics::mtext(my.unit, 4, line = 3, las = 3)
        graphics::mtext(myTitle, 3, line = 1, las = 1)

        # bias
        my.ylab <- "bias"
        my.unit <- latex2exp::TeX(units[i])
        colname <- paste(names[i], ".bias", sep = "")
        bias <- data[colname]
        zone <- data[1]
        bias <- data.frame(data[1], bias[1])
        colnames(bias) <- c("zone", "bias")
        plot(bias$zone, bias$bias, type = "l", main = NA, ylab = my.ylab, xlab = NA, xaxt = "n", panel.first = graphics::rect(-23,
            min(bias$bias, na.rm = TRUE), 23, max(bias$bias, na.rm = TRUE), density = NULL, angle = 0, col = "grey90", border = NA))
        graphics::axis(1, at = seq(lat.range[1], lat.range[2], 10), labels = FALSE, tcl = 0.3)
        graphics::axis(4, labels = FALSE, tcl = 0.3)
        graphics::legend("topleft", "(b)")
        graphics::lines(lat.range, c(0, 0), lty = 2)
        graphics::mtext(my.unit, 4, line = 1, las = 3)

        # crmse
        my.ylab <- "crmse"
        my.unit <- latex2exp::TeX(units[i])
        colname <- paste(names[i], ".crmse", sep = "")
        crmse <- data[colname]
        zone <- data[1]
        crmse <- data.frame(data[1], crmse[1])
        colnames(crmse) <- c("zone", "crmse")
        plot(crmse$zone, crmse$crmse, type = "l", main = NA, xlab = NA, ylab = "crmse", xaxt = "n", yaxt = "n", panel.first = graphics::rect(-23,
            min(crmse$crmse, na.rm = TRUE), 23, max(crmse$crmse, na.rm = TRUE), density = NULL, angle = 0, col = "grey90",
            border = NA))
        graphics::legend("topleft", "(c)")
        graphics::axis(1, at = seq(lat.range[1], lat.range[2], 10), labels = FALSE, tcl = 0.3)
        graphics::axis(2, labels = FALSE, tcl = 0.3)
        graphics::axis(4, labels = TRUE, tcl = 0.3)
        graphics::mtext(my.unit, 4, line = 3, las = 3)

        # phase
        my.ylab <- "phase bias"
        my.unit <- "months"
        colname <- paste(names[i], ".phase", sep = "")
        phase <- data[colname]
        zone <- data[1]
        phase <- data.frame(data[1], phase[1])
        colnames(phase) <- c("zone", "phase")
        plot(phase$zone, phase$phase, type = "l", main = NA, xlab = NA, ylab = my.ylab, xaxt = "n", panel.first = graphics::rect(-23,
            min(phase$phase, na.rm = TRUE), 23, max(phase$phase, na.rm = TRUE), density = NULL, angle = 0, col = "grey90",
            border = NA))
        graphics::legend("topleft", "(d)")
        graphics::axis(1, at = seq(lat.range[1], lat.range[2], 10), labels = FALSE, tcl = 0.3)
        graphics::axis(4, labels = FALSE, tcl = 0.3)
        graphics::mtext(my.unit, 4, line = 1, las = 3)

        # iav
        my.ylab <- "iav"
        my.unit <- latex2exp::TeX(units[i])
        colname <- paste(names[i], ".iav", sep = "")
        colname.mod <- paste(names[i], ".mod.iav", sep = "")
        colname.ref <- paste(names[i], ".ref.iav", sep = "")
        iav.mod <- data[colname.mod]
        iav.ref <- data[colname.ref]
        zone <- data[1]
        iav <- data.frame(data[1], iav.mod[1], iav.ref[1])
        colnames(iav) <- c("zone", "iav.mod", "iav.ref")
        y.min <- min(iav$iav.mod, iav$iav.ref, na.rm = TRUE)
        y.max <- max(iav$iav.mod, iav$iav.ref, na.rm = TRUE)
        my.ylim <- c(y.min, y.max)
        plot(iav$zone, iav$iav.mod, type = "l", main = NA, ylim = my.ylim, xlab = NA, ylab = "iav", xaxt = "n", yaxt = "n",
            panel.first = graphics::rect(-23, y.min, 23, y.max, density = NULL, angle = 0, col = "grey90", border = NA))
        graphics::lines(iav$zone, iav$iav.ref, col = "red")
        graphics::legend("topleft", "(e)")
        graphics::legend("topright", c("model", "reference"), col = c("black", "red"), pch = 16, horiz = TRUE, bty = "n")
        graphics::axis(1, at = seq(lat.range[1], lat.range[2], 10), labels = FALSE, tcl = 0.3)
        graphics::axis(2, labels = FALSE, tcl = 0.3)
        graphics::axis(4, labels = TRUE, tcl = 0.3)
        graphics::mtext(my.unit, 4, line = 3, las = 3)

        # scores
        my.ylab <- "score (-)"
        my.unit <- "dimensionless"
        colname.bias.score <- paste(names[i], ".bias.score", sep = "")
        colname.rmse.score <- paste(names[i], ".rmse.score", sep = "")
        colname.phase.score <- paste(names[i], ".phase.score", sep = "")
        colname.iav.score <- paste(names[i], ".iav.score", sep = "")

        bias.score <- data[colname.bias.score]
        rmse.score <- data[colname.rmse.score]
        phase.score <- data[colname.phase.score]
        iav.score <- data[colname.iav.score]

        zone <- data[1]
        scores <- data.frame(data[1], bias.score[1], rmse.score[1], phase.score[1], iav.score[1])
        colnames(scores) <- c("zone", "bias.score", "rmse.score", "phase.score", "iav.score")

        plot(scores$zone, scores$bias.score, type = "l", ylim = c(0, 1), xlab = "Degrees Latitude", xaxt = "n", ylab = my.ylab,
            panel.first = graphics::rect(-23, 0, 23, 1, density = NULL, angle = 0, col = "grey90", border = NA))

        graphics::lines(scores$zone, scores$rmse.score, col = "red")
        graphics::lines(scores$zone, scores$phase.score, col = "blue")
        graphics::lines(scores$zone, scores$iav.score, col = "darkgreen")
        graphics::legend("topleft", "(f)")
        graphics::legend("bottomleft", expression("s"[bias], "s"[rmse], "s"[phase], "s"[iav]), col = c("black", "red", "blue",
            "darkgreen"), pch = 16, horiz = TRUE, bty = "n")
        graphics::axis(1, at = seq(lat.range[1], lat.range[2], 10), labels = TRUE, tcl = 0.3)
        graphics::axis(4, labels = FALSE, tcl = 0.3)
        graphics::mtext(my.unit, 4, line = 1, las = 3)
        grDevices::dev.off()
    }

}

utils::globalVariables("%do%")

