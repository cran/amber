################################################################################
#' Compares scores from two model runs
#' @description This function compares scores from two different model runs.
#' @param mod01.path A string that gives the path where the output from \link{scores.tables} is stored (model 1)
#' @param mod02.path A string that gives the path where the output from \link{scores.tables} is stored (model 2)
#' @param mod01.id A string that gives the name of a model, e.g. 'CanESM2'
#' @param mod02.id A string that gives the name of a model, e.g. 'CanESM5'
#' @param score.weights R object that gives the weights of each score (\eqn{S_{bias}}, \eqn{S_{rmse}}, \eqn{S_{phase}}, \eqn{S_{iav}}, \eqn{S_{dist}})
#' that are used for computing the overall score, e.g. c(1,2,1,1,1)
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @return A figure in PDF format that shows the (a) skill score and (b) skill
#' score difference when compared against a different model run. White numbers
#' indicate that score differences are not statistically significant at the
#' 5 percent level using the two-sided Wilcox significance test.
#' @examples
#'
#' library(amber)
#' library(raster)
#'
#' mod01.path <- system.file('extdata/model01', package = 'amber')
#' mod02.path <- system.file('extdata/model02', package = 'amber')
#' mod01.id <- 'Model01'
#' mod02.id <- 'Model02'
#' score.weights <- c(1,2,1,1,1)
#' scores.compare(mod01.path, mod02.path, mod01.id, mod02.id,  score.weights)
#'
#' @export
scores.compare <- function(mod01.path, mod02.path, mod01.id, mod02.id, score.weights = c(1, 2, 1, 1, 1), outputDir = FALSE) {

    # get variable and reference data name
    my.list <- list.files(path = mod01.path, pattern = "allscorevalues-")
    split <- strsplit(my.list, "-")
    split <- unlist(split)
    split <- matrix(split, nrow = 3)
    split <- t(split)
    split <- as.data.frame(split)
    split <- split[c(2, 3)]  # select columns
    colnames(split) <- c("variable", "ref.id")
    id <- paste(split$variable, split$ref.id, sep = "-")

    # Make data frames that give mean score values
    data.list <- lapply(paste(mod01.path, my.list, sep = "/"), utils::read.table)
    mod01.list <- data.list  # used later in significance test

    # calculate column mean for each element of the list
    data <- lapply(data.list, colMeans, na.rm = TRUE)

    # convert list to matrix
    data <- do.call("rbind", data)
    df <- data.frame(split, data)

    # function for computing overall score
    fun.overall.score <- function(scores) {
        mask <- scores - scores + 1
        w <- mask * score.weights
        overall.score <- sum(w * scores, na.rm = TRUE)/sum(w, na.rm = TRUE)
        return(overall.score)
    }

    # compute overall score and add to data frame
    scores <- df[, 3:7]
    overall.scores <- apply(scores, 1, fun.overall.score)
    df <- data.frame(df, overall.scores)
    mod01.scores <- df

    # get data from second run
    data.list <- lapply(paste(mod02.path, my.list, sep = "/"), utils::read.table)
    mod02.list <- data.list  # used later in significance test

    # calculate column mean for each element of the list
    data <- lapply(data.list, colMeans, na.rm = TRUE)

    # convert list to matrix
    data <- do.call("rbind", data)
    df <- data.frame(split, data)

    # compute overall score and add to data frame
    scores <- df[, 3:7]
    overall.scores <- apply(scores, 1, fun.overall.score)
    df <- data.frame(df, overall.scores)
    mod02.scores <- df

    # Compare both runs using Wilcox significance test
    datalist <- list()
    for (i in 1:length(mod01.list)) {
        # mod01
        mod01 <- mod01.list[i]
        mod01 <- data.frame(mod01)
        mod01.bias <- stats::na.exclude(mod01$bias.score)
        mod01.rmse <- stats::na.exclude(mod01$rmse.score)
        mod01.phase <- stats::na.exclude(mod01$phase.score)
        mod01.iav <- stats::na.exclude(mod01$iav.score)
        mod01.dist <- stats::na.exclude(mod01$dist.score)
        mod01.overall <- c(mod01.bias, mod01.rmse, mod01.phase, mod01.iav, mod01.dist[1])

        # mod02
        mod02 <- mod02.list[i]
        mod02 <- data.frame(mod02)
        mod02.bias <- stats::na.exclude(mod02$bias.score)
        mod02.rmse <- stats::na.exclude(mod02$rmse.score)
        mod02.phase <- stats::na.exclude(mod02$phase.score)
        mod02.iav <- stats::na.exclude(mod02$iav.score)
        mod02.dist <- stats::na.exclude(mod02$dist.score)
        mod02.overall <- c(mod02.bias, mod02.rmse, mod02.phase, mod02.iav, mod02.dist[1])

        # p-value
        bias.p <- ifelse(length(mod01.bias) > 0, stats::wilcox.test(mod01.bias, mod02.bias, alternative = c("two.sided"))$p.value,
            NA)
        rmse.p <- ifelse(length(mod01.rmse) > 0, stats::wilcox.test(mod01.rmse, mod02.rmse, alternative = c("two.sided"))$p.value,
            NA)
        phase.p <- ifelse(length(mod01.phase) > 0, stats::wilcox.test(mod01.phase, mod02.phase, alternative = c("two.sided"))$p.value,
            NA)
        iav.p <- ifelse(length(mod01.iav) > 0, stats::wilcox.test(mod01.iav, mod02.iav, alternative = c("two.sided"))$p.value, NA)
        dist.p <- 1  # since S_dist is a single number for entire globe
        overall.p <- ifelse(length(mod01.overall) > 0, stats::wilcox.test(mod01.overall, mod02.overall, alternative = c("two.sided"))$p.value,
            NA)
        pvalues <- data.frame(bias.p, rmse.p, phase.p, iav.p, dist.p, overall.p)
        datalist[[i]] <- pvalues
    }
    pvalues <- do.call(rbind, datalist)
    pvalues <- round(pvalues, 2)
    rownames(pvalues) <- id

    # convert to raster
    data <- mod01.scores
    data <- data[-c(1, 2)]  # drop first two columns, i.e. variable and ref.id
    data <- as.matrix(data)
    data <- raster::raster(data)
    mod01 <- data

    # mod02
    data <- mod02.scores
    data <- data[-c(1, 2)]  # drop first two columns, i.e. variable and ref.id
    data <- as.matrix(data)
    data <- raster::raster(data)
    mod02 <- data

    # delta
    delta <- mod02 - mod01

    # pvalues
    data <- pvalues
    data <- as.matrix(data)
    data <- raster::raster(data)
    pvalues <- data
    pmask <- pvalues
    pmask[pvalues > 0.05] <- NA
    pmask <- pmask - pmask + 1
    delta.sig <- delta * pmask

    # get statistically significant differences

    # prepare inputs for plot
    my.col.lab <- latex2exp::TeX(c("$S_{bias}$", "$S_{rmse}$", "$S_{phase}$", "$S_{iav}$", "$S_{dist}$", "$S_{overall}$"))
    my.row.lab <- id

    # raster extent
    my.width <- 2
    data <- delta
    my.extent <- c(0, my.width * ncol(data), 0.5, nrow(data) + 0.5)
    raster::extent(delta) <- my.extent
    raster::extent(delta.sig) <- my.extent
    raster::extent(mod02) <- my.extent

    # plot
    oldpar <- graphics::par(mfrow = c(1,2))
    on.exit(graphics::par(oldpar))
    if (outputDir != FALSE) {
        grDevices::pdf("scores_compare.pdf", width = 8.1, height = 7)
    }
    graphics::par(mfrow = c(1, 2), font.main = 1, oma = c(0, 12, 0, 0), mar = c(10, 0, 3, 1), lwd = 1, cex = 1, family = "Helvetica")

    # (a) plot mod02

    data <- mod02
    # legend
    mmi.delta <- intFun.min.max.int(data)
    min <- mmi.delta[1]
    max <- mmi.delta[2]
    interval <- mmi.delta[3]
    my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
    my.labels <- round(seq(min, max, interval), 3)  # location of labels
    my.col <- viridis::plasma(n = length(my.breaks) - 1, direction = -1)
    legend.bar.text <- "score (-)"
    my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
    my.legend.args <- list(text = legend.bar.text, side = 3, font = 1, line = 1, cex = 1)

    plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = paste("(a)", mod02.id, sep = " "), axes = FALSE)
    raster::text(data, digits = 2, cex = 0.7)
    axis(side = 2, at = rev(seq(1, nrow(data), 1)), labels = my.row.lab, las = 2)
    axis(side = 1, at = seq(1, my.width * ncol(data), my.width), labels = my.col.lab, las = 2)
    plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args, legend.width = 1.5,
        legend.shrink = 1, font = 1, horizontal = TRUE)

    # (b) plot delta, i.e. mod02-mod01

    data <- delta
    # legend
    mmi.delta <- intFun.min.max.int.bias(data)
    min <- mmi.delta[1]
    max <- mmi.delta[2]
    interval <- mmi.delta[3]
    my.breaks <- round(seq(min, max, interval), 3)  # location of color breaks
    my.labels <- round(seq(min, max, interval), 3)  # location of labels
    my.col <- scico::scico(n = length(my.breaks) - 1, palette = "vik")
    # legend.bar.text <- expression(paste(Delta, 'score (-)', sep = ''))
    legend.bar.text <- "score difference (-)"
    my.axis.args <- list(at = my.labels, labels = my.labels, cex.axis = 1)
    my.legend.args <- list(text = legend.bar.text, side = 3, font = 1, line = 1, cex = 1)

    # plot
    plot(data, col = my.col, breaks = my.breaks, legend = FALSE, main = paste("(b)", mod02.id, "-", mod01.id, sep = " "), axes = FALSE)
    raster::text(delta, digits = 2, cex = 0.7, col = "white")
    plot(delta.sig, col = my.col, breaks = my.breaks, legend = FALSE, axes = FALSE, add = TRUE)
    ngc <- raster::ncell(delta.sig)  # number of grid cells
    nas <- raster::cellStats(is.na(delta.sig), "sum")  # number of NA's
    if (ngc > nas)
        raster::text(delta.sig, digits = 2, cex = 0.7)  # plot values if data has values
    # axis(side=2, at=rev(seq(1,nrow(data),1)), labels = my.row.lab, las=2)
    axis(side = 1, at = seq(1, my.width * ncol(data), my.width), labels = my.col.lab, las = 2)
    plot(data, legend.only = TRUE, col = my.col, breaks = my.breaks, axis.args = my.axis.args, legend.args = my.legend.args, legend.width = 1.5,
        legend.shrink = 1, font = 1, horizontal = TRUE)
    if (outputDir != FALSE) {
        grDevices::dev.off()
    }
}
