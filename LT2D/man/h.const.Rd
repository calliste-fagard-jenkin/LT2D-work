% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2DLTfunctions mixture.r
\name{h.const}
\alias{h.const}
\title{Detection hazard function \code{h.const} prob(detect | available at x,y)}
\usage{
h.const(y, x, b = 1)
}
\arguments{
\item{y}{Forward distance}

\item{x}{perpendicular distance}

\item{b}{single value parameter vector (giving probability of detection).}
}
\value{
probability of detection given that an animal is available at location x,y
}
\description{
This function is for constant detectability throughout
}
\examples{
h.const(0.5,0.5,b=1)
}
\seealso{
\code{\link{h1}} \code{\link{h2}} \code{\link{h.exp2}} \code{\link{h.okamura}}
}
