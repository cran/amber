% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intFun.R
\name{intFun.grid.phase}
\alias{intFun.grid.phase}
\title{Time distance in months}
\usage{
intFun.grid.phase(x)
}
\arguments{
\item{x}{A raster object of type (a)}
}
\value{
A raster object of type (b)
}
\description{
Consider two rasters that show the month of the maximum value of
a variable during the climatological mean annual cycle. Calculating the
absolute difference between both rasters may yield values that range from
0 to 11 (type a). A more meaningful range that considers the circular nature
of the annual cycle would range from 0 to 6 (type b). This function converts
a raster of type (a) to type (b).
}
\examples{

library(raster)
# make some data
data <- runif(100,0,11)
data <- matrix(data, ncol=10)
data <- raster::raster(data)
# convert raster of type (a) to type (b)
data <- calc(data,intFun.grid.phase)
plot(data)
plot(data); text(data)

}
\keyword{internal}
