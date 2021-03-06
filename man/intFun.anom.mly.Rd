% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intFun.R
\name{intFun.anom.mly}
\alias{intFun.anom.mly}
\title{Monthly anomalies}
\usage{
intFun.anom.mly(mly, clim)
}
\arguments{
\item{mly}{An R object with a monthly time series}

\item{clim}{An R object with climatological monthly mean data}
}
\value{
A raster object with monthly anomalies
}
\description{
This function computes monthly anomalies.
}
\examples{

# make some data
month <- seq(1,12,1)
month <- rep(month,10)
mod <- runif(length(month), 0,10)
mod <- data.frame(month, mod)
# make an index
index <- list(mod$month)
clim.mly <- apply(mod, 2, function(x) {tapply(x, index, mean, na.rm=TRUE)})
mod.anom <- intFun.anom.mly(mod, clim.mly)

}
\keyword{internal}
