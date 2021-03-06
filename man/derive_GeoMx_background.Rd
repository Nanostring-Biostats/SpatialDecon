% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_GeoMx_background.R
\name{derive_GeoMx_background}
\alias{derive_GeoMx_background}
\title{Derive background at the scale of the normalized data for GeoMx data}
\usage{
derive_GeoMx_background(norm, probepool, negnames)
}
\arguments{
\item{norm}{Matrix of normalized data, genes in rows and segments in columns.
Must include negprobes, and must have rownames.}

\item{probepool}{Vector of probe pool names for each gene, aligned to the
rows of "norm".}

\item{negnames}{Names of all negProbes in the dataset. Must be at least one
neg.name within each probe pool.}
}
\value{
A matrix of expected background values, in the same scale and
 dimensions as the "norm" argument.
}
\description{
Estimates per-datapoint background levels from a GeoMx experiment.
In studies with two or more probe pools, different probes will have different
background levels. This function provides a convenient way to account for
this phenomenon.
}
\examples{
data(mini_geomx_dataset)
# estimate background:
mini_geomx_dataset$bg <- derive_GeoMx_background(
  norm = mini_geomx_dataset$normalized,
  probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
  negnames = "NegProbe"
)
}
