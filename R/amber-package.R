#' Overview of the functions in the amber package
#'
#' The Canadian Land Surface Scheme Including Biogeochemical Cycles (CLASSIC) is
#' the land surface component of the Canadian Earth System Model (CanESM) (Melton and Arora, 2016).
#' The Automated Model Benchmarking (AMBER) package evaluates the ability of CLASSIC
#' to reproduce land surface processes by comparing model outputs against quasi-observational
#' data sets derived from remote sensing products, eddy covariance flux tower measurements,
#' and stream flow measurements. To summarize model performance across different
#' statistical metrics AMBER employs a skill score system that was originally developed
#' by the International Land Model Benchmarking (ILAMB) framework (Collier et al., 2018).
#' The amber package was created to tailor the ILAMB skill score approach for CLASSIC
#' model outputs.
#'
#' @section amber functions:
#' The functions provided by the amber-package can be grouped into two categories:
#' \itemize{
#' \item Functions that compute skill scores. This includes all functions
#'  with the prefix \emph{scores}. Different functions apply to different reference data,
#'   such as spatially gridded data (\code{\link{scores.grid.time}}), eddy covariance flux tower measurements
#'    (\code{\link{scores.fluxnet.csv}}), or stream flow measurements (\code{\link{scores.runoff}}).
#'
#'  \item Functions for data visualization. This includes
#'  \code{\link{plotGrid}},
#'  \code{\link{plotNc}},
#'  \code{\link{plotNcIrreg}},
#'  \code{\link{seasonalCycle}},
#'  \code{\link{seasonalCycleIrreg}},
#'  \code{\link{zonalMean}}, and
#'  \code{\link{zonalMeanIrreg}}.
#'  The suffix \emph{Irreg} indicates that this function may be applied to data
#'  that is on an irregular grid. In the case of CLASSIC, this applies to high-resolution
#'  simulations that are conducted for the Canadian domain.
#' }
#'
#' @section Skill scores:
#' The performance of a model is expressed through scores that range from zero to
#' one, where increasing values imply better performance. These scores are computed
#' in five steps:
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
#' Melton, J. R., and V. K. Arora. 2016. “Competition between Plant Functional
#' Types in the Canadian Terrestrial Ecosystem Model (CTEM) v. 2.0.” Geoscientific
#' Model Development 9 (1): 323–61.
#'
#' @docType package
#' @name amber-package
NULL
# > NULL
