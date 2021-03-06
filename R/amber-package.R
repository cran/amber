#' Overview of AMBER functions
#'
#' The Canadian Land Surface Scheme Including Biogeochemical Cycles (CLASSIC) is
#' the land surface component of the Canadian Earth System Model (CanESM) (Melton et al, 2020).
#' The Automated Model Benchmarking R package (AMBER) evaluates the ability of CLASSIC
#' to reproduce land surface processes by comparing model outputs against quasi-observational
#' data sets derived from remote sensing products, eddy covariance flux tower measurements,
#' and stream flow measurements. To summarize model performance across different
#' statistical metrics AMBER employs a skill score system that was originally developed
#' by the International Land Model Benchmarking (ILAMB) framework (Collier et al., 2018).
#' AMBER was created to tailor the ILAMB skill score approach for CLASSIC
#' model outputs. While AMBER was tested for CLASSIC only it may also work for other models.
#'
#' @section amber functions:
#' The functions provided by AMBER can be grouped into three categories:
#' \itemize{
#'
#' \item Functions that compute skill scores and other metrics for a single variable.
#' This includes:
#' \code{\link{scores.fluxnet.csv}},
#' \code{\link{scores.fluxnet.nc}},
#' \code{\link{scores.fluxnet.site}},
#' \code{\link{scores.functional.response}},
#' \code{\link{scores.grid.notime}},
#' \code{\link{scores.grid.time}},
#' \code{\link{scores.runoff}}, and
#' \code{\link{scores.site.notime}}.
#'
#' \item Functions that visualize model output for a single variable.
#' This includes:
#' \code{\link{plotEnsembleHovmoeller}},
#' \code{\link{plotEnsembleMean}},
#' \code{\link{plotGrid}},
#' \code{\link{plotHovmoeller}},
#' \code{\link{plotNcIrreg}},
#' \code{\link{plotNc}},
#' \code{\link{plotZonalMeans}},
#' \code{\link{plotZonalMeanStats}},
#' \code{\link{seasonalCycleIrreg}},
#' \code{\link{seasonalCycle}},
#' \code{\link{zonalMeanIrreg}},
#' \code{\link{zonalMean}}, and
#' \code{\link{zonalMeanStats}}.
#'
#' \item Functions that visualize summary statistics across multiple variables.
#' This includes:
#' \code{\link{scores.compare.benchmarks}},
#' \code{\link{scores.compare.ensemble}},
#' \code{\link{scores.compare}},
#' \code{\link{scores.tables}},
#' \code{\link{scores.tables.tweak}},
#' \code{\link{correlationMatrixDiff}},
#' \code{\link{correlationMatrixFluxnet}},
#' \code{\link{correlationMatrix}},
#' \code{\link{globalSumsTable}},
#' \code{\link{metrics.compare}},
#' \code{\link{plotBars}}, and
#' \code{\link{plotFluxnetStats}}.
#'
#'  The suffix \emph{Irreg} indicates that this function may be applied to data
#'  that is on an irregular grid. In the case of CLASSIC, this applies to high-resolution
#'  simulations that are conducted for the Canadian domain. The string \emph{ensemble}
#'  implies that the function is designed to handle multiple model runs.
#'
#' }
#'
#' @section Skill scores:
#' The performance of a model is expressed through scores that range from zero to
#' one, where increasing values imply better performance. These scores are usually
#' computed in five steps:
#' \itemize{
#' \item \eqn{(i)} computation of a statistical metric,
#' \item \eqn{(ii)} nondimensionalization,
#' \item \eqn{(iii)} conversion to unit interval,
#' \item \eqn{(iv)} spatial integration, and
#' \item \eqn{(v)} averaging scores computed from different statistical metrics.
#' }
#' The latter includes the bias, root-mean-square error, phase shift,
#' inter-annual variability, and spatial distribution. The equations for computing
#' the bias score (\eqn{S_{bias}}) are:
#'
#' \eqn{(i) \ bias(\lambda, \phi)=\overline{v_{mod}}(\lambda, \phi)-\overline{v_{ref}}(\lambda, \phi)}
#'
#' \eqn{(ii) \ \varepsilon_{bias}=|bias(\lambda, \phi)|/\sigma_{ref}(\lambda, \phi)}
#'
#' \eqn{(iii) \ s_{bias}(\lambda, \phi)=e^{-\varepsilon_{bias}(\lambda, \phi)}}
#'
#' \eqn{(iv) \ S_{bias}=\overline{\overline{s_{bias}}}}
#'
#' The equations for computing the root mean square error score (\eqn{S_{rmse}}) are:
#'
#' \eqn{(i) \ crmse(\lambda, \phi) =\sqrt{\frac{1}{t_{f}-t_{0}}\int_{t_{0}}^{t_{f}}[(v_{mod}(t,\lambda, \phi)-\overline{v_{mod}}(\lambda, \phi))-(v_{ref}(t,\lambda, \phi)-\overline{v_{ref}}(\lambda, \phi))]^{2}dt}}
#'
#' \eqn{(ii) \ \varepsilon_{rmse}(\lambda, \phi)=crmse(\lambda, \phi)/\sigma_{ref}(\lambda, \phi)}
#'
#' \eqn{(iii) \ s_{rmse}(\lambda, \phi)=e^{-\varepsilon_{rmse}(\lambda, \phi)}}
#'
#' \eqn{(iv) \ S_{rmse}=\overline{\overline{s_{rmse}}}}
#'
#' The equations for computing the phase score (\eqn{S_{phase}}) are:
#'
#' \eqn{(i) \ \theta(\lambda, \phi)=\max(c_{mod}(t,\lambda, \phi))-\max(c_{ref}(t,\lambda, \phi))}
#'
#' \eqn{(ii) \ \textrm{not applicable, as units are consistent across all variables}}
#'
#' \eqn{(iii) \ s_{phase}(\lambda, \phi)=\frac{1}{2}[1+\cos(\frac{2\pi\theta(\lambda, \phi)}{365})]}
#'
#' \eqn{(iv) \ S_{phase}=\overline{\overline{s_{phase}}}}
#'
#' The equations for computing the inter-annual variability score (\eqn{S_{iav}}) are:
#'
#' \eqn{(i) \ iav_{ref}(\lambda, \phi)=\sqrt{\frac{1}{t_{f}-t_{0}}\int_{t_{0}}^{t_{f}}(v_{ref}(t,\lambda, \phi)-c_{ref}(t,\lambda, \phi))^{2}dt}}
#'
#' \eqn{(i) \ iav_{mod}(\lambda, \phi)=\sqrt{\frac{1}{t_{f}-t_{0}}\int_{t_{0}}^{t_{f}}(v_{mod}(t,\lambda, \phi)-c_{mod}(t,\lambda, \phi))^{2}dt}}
#'
#' \eqn{(ii) \ \varepsilon_{iav}=|(iav_{mod}(\lambda, \phi)-iav_{ref}(\lambda, \phi))|/iav_{ref}(\lambda, \phi)}
#'
#' \eqn{(iii) \ s_{iav}(\lambda, \phi)=e^{-\varepsilon_{iav}(\lambda, \phi)}}
#'
#' \eqn{(iv) \ S_{iav}=\overline{\overline{s_{iav}}}}
#'
#' The equations for computing the spatial distribution score (\eqn{S_{dist}}) are:
#'
#' \eqn{(i) \ \sigma=\sigma_{\overline{v_{mod}}}/\sigma_{\overline{v_{ref}}}}
#'
#' \eqn{(ii)} and \eqn{ (iii) \ \textrm{not applicable}}
#'
#' \eqn{(iv) \ S_{dist}=2(1+R)/(\sigma+\frac{1}{\sigma})^{2}}
#'
#' where \eqn{\overline{v_{mod}}(\lambda, \phi)} and
#' \eqn{\overline{v_{ref}}(\lambda, \phi)} are the mean values in time \eqn{t}
#' of a variable {v} as a function of longitude \eqn{\lambda} and latitude
#' \eqn{\phi} for model and reference data, respectively,
#' \eqn{t_0} and \eqn{t_f} are the initial and final time step,
#' \eqn{\sigma_{\overline{v_{mod}}}} and \eqn{\sigma_{\overline{v_{ref}}}} are the
#' standard deviation of the time mean values from the model and reference data,
#' and \eqn{R} is the spatial correlation coefficient of
#' \eqn{\overline{v_{ref}}(\lambda, \phi)} and
#' \eqn{\overline{v_{mod}}(\lambda, \phi)}.
#'
#' Score values are then combined to derive a single overall score for each output variable:
#'
#' \eqn{(v) \ S_{overall}=\frac{S_{bias}+2S_{rmse}+S_{phase}+S_{iav}+S_{dist}}{1+2+1+1+1}.}
#'
#' Note that \eqn{S_{rmse}} is weighted by a factor of two, which emphasizes its importance.
#'
#' @section References:
#' Collier, Nathan, Forrest M. Hoffman, David M. Lawrence, Gretchen Keppel-Aleks,
#' Charles D. Koven, William J. Riley, Mingquan Mu, and James T. Randerson. 2018.
#' “The International Land Model Benchmarking (ILAMB) System: Design, Theory, and
#' Implementation.” Journal of Advances in Modeling Earth Systems 10 (11): 2731–54.
#'
#' Melton, Joe R., Vivek K. Arora, Eduard Wisernig-Cojoc, Christian Seiler,
#' Matthew Fortier, Ed Chan, and Lina Teckentrup. 2020. “CLASSIC v1.0:
#' The Open-Source Community Successor to the Canadian Land Surface Scheme (CLASS)
#' and the Canadian Terrestrial Ecosystem Model (CTEM) - Part 1: Model Framework
#' and Site-Level Performance.” https://doi.org/10.5194/gmd-13-2825-2020.
#'
#' @docType package
#' @name amber-package
NULL
# > NULL
