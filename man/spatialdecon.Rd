% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialdecon.R
\name{spatialdecon}
\alias{spatialdecon}
\title{Mixed cell deconvolution of spatiall-resolved gene expression data}
\usage{
spatialdecon(
  norm,
  bg,
  X = NULL,
  raw = NULL,
  wts = NULL,
  resid_thresh = 3,
  lower_thresh = 0.5,
  align_genes = TRUE,
  is_pure_tumor = NULL,
  n_tumor_clusters = 10,
  cell_counts = NULL,
  cellmerges = NULL,
  maxit = 1000
)
}
\arguments{
\item{norm}{p-length expression vector or p * N expression matrix - the
actual (linear-scale) data}

\item{bg}{Same dimension as norm: the background expected at each data point.}

\item{X}{Cell profile matrix. If NULL, the safeTME matrix is used.}

\item{raw}{Optional for using an error model to weight the data points.
p-length expression vector or p * N expression matrix - the raw
(linear-scale) data}

\item{wts}{Optional, a matrix of weights.}

\item{resid_thresh}{A scalar, sets a threshold on how extreme individual data
points' values
can be (in log2 units) before getting flagged as outliers and set to NA.}

\item{lower_thresh}{A scalar. Before log2-scale residuals are calculated,
both observed and fitted
values get thresholded up to this value. Prevents log2-scale residuals from
becoming extreme in
points near zero.}

\item{align_genes}{Logical. If TRUE, then Y, X, bg, and wts are row-aligned
by shared genes.}

\item{is_pure_tumor}{A logical vector denoting whether each AOI consists of
pure tumor. If specified,
then the algorithm will derive a tumor expression profile and merge it with
the immune profiles matrix.}

\item{n_tumor_clusters}{Number of tumor-specific columns to merge into the
cell profile matrix.
 Has an impact only when is_pure_tumor argument is used to indicate pure
  tumor AOIs.
 Takes this many clusters from the pure-tumor AOI data and gets the average
 expression profile in each cluster.  Default 10.}

\item{cell_counts}{Number of cells estimated to be within each sample. If
provided alongside norm_factors,
 then the algorithm will additionally output cell abundance esimtates on the
 scale of cell counts.}

\item{cellmerges}{A list object holding the mapping from beta's cell names to
combined cell names. If left
NULL, then defaults to a mapping of granular immune cell definitions to
 broader categories.}

\item{maxit}{Maximum number of iterations. Default 1000.}
}
\value{
a list:
\itemize{
\item beta: matrix of cell abundance estimates, cells in rows and
observations in columns
\item sigmas: covariance matrices of each observation's beta estimates
\item p: matrix of p-values for H0: beta == 0
\item t: matrix of t-statistics for H0: beta == 0
\item se: matrix of standard errors of beta values
\item prop_of_all: rescaling of beta to sum to 1 in each observation
\item prop_of_nontumor: rescaling of beta to sum to 1 in each observation,
excluding tumor abundance estimates
\item cell.counts: beta rescaled to estimate cell numbers, based on
prop_of_all and nuclei count
\item beta.granular: cell abundances prior to combining closely-related
cell types
\item sigma.granular: sigmas prior to combining closely-related cell types
\item cell.counts.granular: cell.counts prior to combining closely-related
cell types
\item resids: a matrix of residuals from the model fit.
 (log2(pmax(y, lower_thresh)) - log2(pmax(xb, lower_thresh))).
\item X: the cell profile matrix used in the decon fit.
}
}
\description{
Runs the spatialdecon algorithm with added optional functionalities.
Workflow is:
\enumerate{
\item compute weights from raw data
\item Estimate a tumor profile and merge it into the cell profiles matrix
\item run deconvolution once
\item remove poorly-fit genes from first round of decon
\item re-run decon with cleaned-up gene set
\item combine closely-related cell types
\item compute p-values
\item rescale abundance estimates, to proportions of total, proportions of
 immune, cell counts
}
}
\examples{
data(mini_geomx_dataset)
data(safeTME)
data(safeTME.matches)
# estimate background:
mini_geomx_dataset$bg <- derive_GeoMx_background(
  norm = mini_geomx_dataset$normalized,
  probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
  negnames = "NegProbe"
)
# run basic decon:
res0 <- spatialdecon(
  norm = mini_geomx_dataset$normalized,
  bg = mini_geomx_dataset$bg,
  X = safeTME
)
# run decon with bells and whistles:
res <- spatialdecon(
  norm = mini_geomx_dataset$normalized,
  bg = mini_geomx_dataset$bg,
  X = safeTME,
  cellmerges = safeTME.matches,
  cell_counts = mini_geomx_dataset$annot$nuclei,
  is_pure_tumor = mini_geomx_dataset$annot$AOI.name == "Tumor"
)
}
