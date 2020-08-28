################################################################################
#' Table of globally summed values and corresponding biases.
#' @description This function produces a table of globally summed values and
#' corresponding biases when comparing model and reference data. The table provides
#' the variable name, the reference data name, globally summed values from the model
#' and reference data, the absolute bias, the relative bias, units, and time period.
#' The inputs consist of netCDF files produced by \link{scores.grid.time} and
#' \link{scores.grid.notime}. The global sums are based on grid cells that model
#' and reference data have in common.
#' @param mod.path.list A list with paths for each model run, e.g. mod.path.list <- list(mod01.path, mod02.path, mod03.path).
#' @param modelIDs An R object with the different model run IDs, e.g. c('CLASSIC.CRUJRAv2', 'CLASSIC.GSWP3W5E5', 'CLASSIC.CRUNCEP')
#' @param variableNames Variable names, e.g. c('GPP', 'CSOIL', 'BURNT')
#' @param unitLaTeX Desired output units in LaTeX format, e.g. c('PgC yr$^{-1}$', 'PgC', '10$^6$ km$^2$ yr^{-1}')
#' @param conversionFactor Factors that convert the unit of the input files into the desired output unit,
#' e.g. c((10^(15) / 365.25)^(-1), 10^(-12), 12*10^(-14)).
#' Examples of typical conversions include (i) carbon fluxes from gC/day to PgC/yr,
#' (ii) carbon stocks from kgC to PgC, and fractional area burnt from percentage per month to 10^6 km2 per year.
#' Conversion from gC/day to PgC/yr: 1 PgC/yr = 10^15 gC/yr = 10^15 / 365.25 gC/day. Taking the reciprocal leads to the conversion factor, i.e. (10^(15) / 365.25)^(-1).
#' Conversion from kgC to PgC: 1 PgC = 10^15 gC = 10^12 kgC. The resulting conversion factor is 10^(-12).
#' Conversion from percentage of grid cell per month to 10^6 km2 per year: 10^6 km2 per year = 10^12 m2 per year = 1/12 * 10^12 m2 per month.
#' To account for the conversion from percentage to fraction, we need to multiply by 100. The resulting conversion factor is 12*10^(-14).
#' @param outputDir A string that gives the output directory, e.g. '/home/project/study'. The output table in LaTeX format will only be written if the user specifies an output directory.
#' @return A table with globally summed values and corresponding biases.
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
#' mod01.path <- paste(system.file('extdata', package = 'amber'), 'model01', sep = '/')
#' mod02.path <- paste(system.file('extdata', package = 'amber'), 'model02', sep = '/')
#'
#' mod.path.list <- list(mod01.path, mod02.path)
#' modelIDs <- c('CLASSIC.CRUJRAv2', 'CLASSIC.GSWP3W5E5')
#'
#' variableNames <- c('GPP', 'BURNT')
#' unitLaTeX <- c('PgC yr$^{-1}$', '10$^6$ km$^2$ yr^{-1}')
#' conversionFactor <- c((10^(15) / 365.25)^(-1), 12*10^(-14))
#'
#' globalSumsTable(mod.path.list, modelIDs, variableNames, unitLaTeX, conversionFactor,
#' outputDir = FALSE)
#'
#' @export
globalSumsTable <- function(mod.path.list, modelIDs, variableNames, unitLaTeX, conversionFactor, outputDir = FALSE) {

    nmod <- length(mod.path.list)  # number of model runs
    variableName <- variableNames
    varUniCon <- data.frame(variableName, unitLaTeX, conversionFactor)

    # make a list of netCDF files from the first model folder
    eachVariable <- foreach::foreach(v = variableNames) %do% {
        myPattern <- paste(v, "*-mean.nc", sep = "")
        nc.files <- list.files(path = mod.path.list[[1]], pattern = utils::glob2rx(myPattern))
    }
    nc.files <- unlist(eachVariable)

    #

    #---------------------------------------------------------------------------

    # process those files for each model

    #---------------------------------------------------------------------------

    eachModel <- foreach::foreach(m = 1:nmod) %do% {

        eachFile <- foreach::foreach(i = 1:length(nc.files)) %do% {
            # get name
            nc.file <- unlist(nc.files[i])
            nc.file <- paste(mod.path.list[m], nc.file, sep = "/")

            # get data
            data <- raster::raster(nc.file)

            # get ref ID and detect wheter model or reference data
            fileName <- nc.files[i]
            fileName <- base::unlist(base::strsplit(fileName, "-"))
            refID <- fileName[2]
            modOrRef <- fileName[3]

            # get meta data
            meta <- names(data)
            meta <- base::unlist(base::strsplit(meta, "_"))
            variableName <- meta[2]
            start <- round(as.numeric(meta[5]), 0)
            end <- round(as.numeric(meta[7]), 0)
            period <- paste(start, end, sep = "-")

            # get values
            data <- raster::raster(nc.file)

            # compute global sum
            area <- raster::area(data) * 1000 * 1000  # area of grid cell in m^2
            dataXarea <- data * area
            globalSum <- sum(raster::values(dataXarea), na.rm = TRUE)
            data.frame(variableName, refID, globalSum, period, modOrRef)
        }
        data <- do.call(rbind, eachFile)
        myOrder <- seq(1, nrow(data), 1)
        data <- data.frame(myOrder, data)
        data <- merge(data, varUniCon, by = "variableName")

        # convert unit
        globalSum <- data$globalSum * data$conversionFactor
        globalSum <- round(globalSum, 2)

        data <- subset(data, select = -c(globalSum, conversionFactor))
        data <- data.frame(data, globalSum)

        mod <- subset(data, modOrRef == "mod")
        ref <- subset(data, modOrRef == "ref")

        data <- merge(mod, ref, by = c("variableName", "refID", "period", "unitLaTeX"))
        data <- data[order(data$myOrder.x), ]

        gs.mod <- data$globalSum.x
        gs.ref <- data$globalSum.y

        bias.abs <- gs.mod - gs.ref
        bias.rel <- (gs.mod - gs.ref)/abs(gs.ref) * 100
        bias.abs <- round(bias.abs, 2)
        bias.rel <- round(bias.rel, 2)

        variable <- data$variableName
        model <- gs.mod
        reference <- gs.ref
        unit <- data$unitLaTeX
        period <- data$period
        ref.id <- data$refID
        model.id <- modelIDs[m]
        mod.id <- rep(model.id, nrow(data))
        myOrder <- seq(1, nrow(data))

        data <- data.frame(myOrder, variable, ref.id, mod.id, reference, model, bias.abs, bias.rel, unit, period)
        data.globalSumsTable <- data

        colnames(data) <- c("myOrder", "Variable", "Ref. ID", "Model ID", "Ref.", "Model", "Bias", "Bias ($\\%$)", "Unit", "Period")
        globalSumsTable <- data
    }
    data <- do.call("rbind", eachModel)
    data <- data[order(data$myOrder), ]
    globalSumsTable <- subset(data, select = -c(myOrder))

    globalSumsTable <- xtable::xtable(globalSumsTable, digits = 2)
    xtable::caption(globalSumsTable) <- "Globally summed mean values and corresponding biases"

    if (outputDir != FALSE) {
        utils::write.table(data.globalSumsTable, "globalSumsTable")  # not in LaTeX format
        xtable::print.xtable(globalSumsTable, include.rownames = FALSE, label = "tab:globalSums", type = "latex", file = "globalSumsTable.tex",
            caption.placement = "top", sanitize.text.function = function(x) {
                x
            })

    }
    return(data)
}
utils::globalVariables("%do%")
utils::globalVariables("v")
utils::globalVariables("m")
