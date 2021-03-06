% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2DLTfunctions mixture.r
\name{h1.to.HB}
\alias{h1.to.HB}
\title{Parameter conversion for Hayes+Buckland hazard rate model.}
\usage{
h1.to.HB(b)
}
\arguments{
\item{b}{parameter vector of the \link{\code{h1}}.}
}
\value{
parameter vector c(a1,b1) for Hayes and Buckland (1983) hazard rate model
g(x)=1-exp(-a1*x^(-(b1-1))).
}
\description{
Converts parameters from the 2-dimensional a*r^{-b} form of the Hayes and Buckland (1983)
hazard rate model used in \link{\code{h1}}, to parameters of the perpendicular distance form
of the model used in conventional distance sampling.
}
\examples{
h1.to.HB(c(-7, 1))
}
\seealso{
\code{\link{h1}}
}
