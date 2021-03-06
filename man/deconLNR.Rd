% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconLNR.R
\name{deconLNR}
\alias{deconLNR}
\title{Deconvolution using logNormReg package to run linear mean model and log error
 model}
\usage{
deconLNR(Y, X, bg = 0, weights = NULL, epsilon = NULL, maxit = 1000)
}
\arguments{
\item{Y}{p-length expression vector or p * N expression matrix - the actual
(linear-scale) data}

\item{X}{p * K Training matrix}

\item{bg}{scalar or matrix of expected background counts per data point.}

\item{weights}{The same as the weights argument used by lm}

\item{epsilon}{optional,  a very small non-zero number to use as a lower
threshold to make fits well-behaved}

\item{maxit}{Maximum number of iterations. Default 1000.}
}
\value{
a list: beta (estimate), sigmas (covariance matrix of estimate,
derived by inverting the hessian from lognlm)
}
\description{
Calls lognlm() to optimize the model.
}
\keyword{internal}
