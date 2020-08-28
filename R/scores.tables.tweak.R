################################################################################
#' Tweak summary table
#' @description This function allows the user to tweak the summary table computed
#' by \link{scores.tables}. Contrary to \link{scores.tables}, this function can be used
#' to create a single summary table that includes the most important metrics only.
#' The user can specify what variables to include and in what order they should appear.
#' @param myVariables An R object with variable names of variables that should be included in table, e.g. c('GPP', 'RECO', 'NEE')
#' @param myCaption A string that is used as table caption, e.g. 'Globally averaged statistical metrics'.
#' @param inputDir A string that gives the input directory, e.g. '/home/project/study'.
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.
#' @return One table in LaTeX format that shows a subset of statistical metrics
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
#' myInputDir <- paste(system.file('extdata', package = 'amber'), 'scores', sep = '/')
#' myVariables <- c('GPP', 'LAI', 'ALBS')
#' scores.tables.tweak(myVariables = myVariables, inputDir = myInputDir)
#' @export
scores.tables.tweak <- function(myVariables, myCaption = "Globally averaged statistical metrics", inputDir = getwd(), outputDir = FALSE) {

    # summary table with globally averaged inputs for computing scores
    my.list <- list.files(path = inputDir, pattern = "scoreinputs_")
    my.files <- paste(inputDir, my.list, sep = "/")
    data <- lapply(my.files, utils::read.table)
    data <- do.call("rbind", data)
    colnames(data)

    myOrder <- seq(1, length(myVariables), 1)
    myOrder <- data.frame(myVariables, myOrder)
    colnames(myOrder) <- c("variable.name", "order")
    data <- merge(data, myOrder, by = "variable.name")
    data <- data[order(data$order, data$ref.id), ]
    data <- subset(data, select = -c(order))

    colnames(data) <- c("Name", "Variable", "Reference", "Unit", "$v_{mod}$", "$v_{ref}$", "Bias", "Bias (\\%)", "$\\sigma_{ref}$",
        "$\\epsilon_{bias}$ (-)", "$S_{bias}$ (-)", "$rmse$", "$crmse$", "$\\sigma_{ref}$", "$\\epsilon_{rmse}$ (-)", "$S_{rmse}$ (-)",
        "$max_{cmod}$", "$max_{cref}$", "$\\theta$ (months)", "$S_{phase}$ (-)", "$iav_{mod}$", "$iav_{ref}$", "$\\epsilon_{iav}$ (-)",
        "$S_{iav}$", "$\\sigma_{\\overline{v_{mod}}}$", "$\\sigma_{\\overline{v_{ref}}}$", "$\\sigma$ (-)", "$R$ (-)", "$S_{dist}$ (-)")
    rownames(data) <- c()  # omit rownames


    # Make a table that only includes selected reference data and metrics
    metricsTable <- data[c(1, 3, 4, 5, 6, 7, 8, 13, 19, 23, 27, 28)]  # selection of variables



    rownames(metricsTable) <- c()  # omit rownames

    # convert to LaTeX
    metricsTable <- xtable::xtable(metricsTable)

    xtable::caption(metricsTable) <- myCaption

    if (outputDir != FALSE) {
        xtable::print.xtable(metricsTable, include.rownames = FALSE, label = "tab:global_stats", type = "latex", file = "metricsTable.tex",
            caption.placement = "top", sanitize.text.function = function(x) {
                x
            })
    }
}

