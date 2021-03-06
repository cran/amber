% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlationMatrixFluxnet.R
\name{correlationMatrixFluxnet}
\alias{correlationMatrixFluxnet}
\title{Correlation matrix for statistical metrics computed by AMBER for FLUXNET data}
\usage{
correlationMatrixFluxnet(metric, inputDir, outputDir = FALSE,
  fileNames = c("GPP_FLUXNET", "HFLS_FLUXNET", "HFSS_FLUXNET",
  "NEE_FLUXNET", "RECO_FLUXNET", "RNS_FLUXNET"),
  significanceLevel = 0.01, plot.width = 8, plot.height = 6.8,
  plot.margin = c(10, 10, 1, 4))
}
\arguments{
\item{metric}{A string that indicates for what statistical metric the correlation matrix should be computed.
Options are 'mod.mean', ref.mean', bias', 'crmse', 'phase', 'bias.score', 'crmse.score', 'phase.score', or 'iav.score'.}

\item{inputDir}{A string that gives the location of text files produced by \link{scores.fluxnet.csv}, e.g. '/home/project/study'.}

\item{outputDir}{A string that gives the output directory, e.g. '/home/project/study'. The output will only be written if the user specifies an output directory.}

\item{fileNames}{An object of strings that give the filenames that should be included. The default is c('GPP_FLUXNET', 'HFLS_FLUXNET', 'HFSS_FLUXNET', 'NEE_FLUXNET', 'RECO_FLUXNET', 'RNS_FLUXNET')}

\item{significanceLevel}{A number that gives the desired significance level of a correlation, e.g. 0.01}

\item{plot.width}{A number that gives the plot width, e.g. 8}

\item{plot.height}{A number that gives the plot height, e.g. 6.8}

\item{plot.margin}{An R object that gives the plot margin, e.g. c(10, 10, 1, 4)}
}
\value{
A list with the Spearman correlation coefficient and corresponding p-values, and a Figure of the correlation matrix
}
\description{
This function produces a correlation matrices for mean values, \eqn{bias},
\eqn{crmse}, \eqn{phase}, and corresponding scores. The input data consist of text files
produced by \link{scores.fluxnet.csv}.
}
\examples{
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

inputDir <- paste(system.file('extdata', package = 'amber'), 'scores', sep = '/')
correlationMatrixFluxnet(metric = 'bias', inputDir)

}
