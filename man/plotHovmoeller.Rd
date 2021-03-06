% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotHovmoeller.R
\name{plotHovmoeller}
\alias{plotHovmoeller}
\title{Plot Hovmoeller diagrams that show monthly climatological mean values and biases}
\usage{
plotHovmoeller(plot.me, long.name, mod.id, ref.id, myBin = 20,
  gridCellWidth = 2, plot.width = 4, plot.height = 5.2,
  my.ylim = c(-100, 100), outputDir = FALSE)
}
\arguments{
\item{plot.me}{A list that is produced by \link{scores.grid.time}}

\item{long.name}{A string that gives the full name of the variable, e.g. 'Gross primary productivity'}

\item{mod.id}{A string that identifies the source of the reference data set, e.g. 'CLASSIC'}

\item{ref.id}{A string that identifies the source of the reference data set, e.g. 'MODIS'}

\item{myBin}{An integer number that defines the latitudinal range used for computing the zonal mean.
For instance, a value of 10 implies that a zonal mean is computed for every 10 degrees latitude.}

\item{gridCellWidth}{A number that is used as a factor to adjust the width of grid cells, e.g. 1.}

\item{plot.width}{Number that gives the plot width, e.g. 8}

\item{plot.height}{Number that gives the plot height, e.g. 8}

\item{my.ylim}{An R object with the latitudinal range that should be plotted, e.g. c(-40, 65).}

\item{outputDir}{A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.}
}
\value{
Figures in PDF format.
}
\description{
This function plots Hovmoeller diagrams of monthly climatological
mean values and biases computed by \link{scores.grid.time}.
}
\examples{

\donttest{
# Global plots on a regular grid

library(amber)
library(classInt)
library(doParallel)
library(foreach)
library(Hmisc)
library(latex2exp)
library(ncdf4)
library(parallel)
library(raster)
library(rgdal)
library(rgeos)
library(scico)
library(sp)
library(stats)
library(utils)
library(viridis)
library(xtable)

long.name <- 'Gross primary productivity'
nc.mod <- system.file('extdata/modelRegular', 'gpp_monthly.nc', package = 'amber')
nc.ref <- system.file('extdata/referenceRegular', 'gpp_GBAF_128x64.nc', package = 'amber')
mod.id <- 'CLASSIC' # define a model experiment ID
ref.id <- 'GBAF' # give reference dataset a name
unit.conv.mod <- 86400*1000 # optional unit conversion for model data
unit.conv.ref <- 86400*1000 # optional unit conversion for reference data
variable.unit <- 'gC m$^{-2}$ day$^{-1}$' # unit after conversion (LaTeX notation)

# Short version using default settings:
plot.me <- scores.grid.time(long.name, nc.mod, nc.ref, mod.id, ref.id, unit.conv.mod,
unit.conv.ref, variable.unit)
plotHovmoeller(plot.me, long.name, mod.id, ref.id)

} #donttest
}
