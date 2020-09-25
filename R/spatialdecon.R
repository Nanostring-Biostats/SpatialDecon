# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
# Copyright (C) 2020, NanoString Technologies, Inc.
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.
# 530 Fairview Avenue N
# Seattle, WA 98109
# Tel: (888) 358-6266
# pdanaher@nanostring.com


#' Mixed cell deconvolution of spatiall-resolved gene expression data
#'
#' Runs the spatialdecon algorithm with added optional functionalities.
#' Workflow is:
#' \enumerate{
#' \item compute weights from raw data
#' \item Estimate a tumor profile and merge it into the cell profiles matrix
#' \item run deconvolution once
#' \item remove poorly-fit genes from first round of decon
#' \item re-run decon with cleaned-up gene set
#' \item combine closely-related cell types
#' \item compute p-values
#' \item rescale abundance estimates, to proportions of total, proportions of
#'  immune, cell counts
#' }
#'
#' @param norm p-length expression vector or p * N expression matrix - the
#' actual (linear-scale) data
#' @param bg Same dimension as norm: the background expected at each data point.
#' @param X Cell profile matrix. If NULL, the safeTME matrix is used.
#' @param raw Optional for using an error model to weight the data points.
#'  p-length expression vector or p * N expression matrix - the raw
#'  (linear-scale) data
#' @param wts Optional, a matrix of weights.
#' @param resid_thresh A scalar, sets a threshold on how extreme individual data
#'  points' values
#'  can be (in log2 units) before getting flagged as outliers and set to NA.
#' @param lower_thresh A scalar. Before log2-scale residuals are calculated,
#'  both observed and fitted
#'  values get thresholded up to this value. Prevents log2-scale residuals from
#'  becoming extreme in
#'  points near zero.
#' @param align_genes Logical. If TRUE, then Y, X, bg, and wts are row-aligned
#'  by shared genes.
#' @param is_pure_tumor A logical vector denoting whether each AOI consists of
#'  pure tumor. If specified,
#'  then the algorithm will derive a tumor expression profile and merge it with
#'  the immune profiles matrix.
#' @param cell_counts Number of cells estimated to be within each sample. If
#' provided alongside norm_factors,
#'  then the algorithm will additionally output cell abundance esimtates on the
#'  scale of cell counts.
#' @param cellmerges A list object holding the mapping from beta's cell names to
#'  combined cell names. If left
#'  NULL, then defaults to a mapping of granular immune cell definitions to
#'   broader categories.
#' @param n_tumor_clusters Number of tumor-specific columns to merge into the
#' cell profile matrix.
#'  Has an impact only when is_pure_tumor argument is used to indicate pure
#'   tumor AOIs.
#'  Takes this many clusters from the pure-tumor AOI data and gets the average
#'  expression profile in each cluster.  Default 10.
#' @param maxit Maximum number of iterations. Default 1000.
#' @return a list:
#' \itemize{
#' \item beta: matrix of cell abundance estimates, cells in rows and
#' observations in columns
#' \item sigmas: covariance matrices of each observation's beta estimates
#' \item p: matrix of p-values for H0: beta == 0
#' \item t: matrix of t-statistics for H0: beta == 0
#' \item se: matrix of standard errors of beta values
#' \item prop_of_all: rescaling of beta to sum to 1 in each observation
#' \item prop_of_nontumor: rescaling of beta to sum to 1 in each observation,
#' excluding tumor abundance estimates
#' \item cell.counts: beta rescaled to estimate cell numbers, based on
#' prop_of_all and nuclei count
#' \item beta.granular: cell abundances prior to combining closely-related
#' cell types
#' \item sigma.granular: sigmas prior to combining closely-related cell types
#' \item cell.counts.granular: cell.counts prior to combining closely-related
#' cell types
#' \item resids: a matrix of residuals from the model fit.
#'  (log2(pmax(y, lower_thresh)) - log2(pmax(xb, lower_thresh))).
#' \item X: the cell profile matrix used in the decon fit.
#' }
#' @examples
#' data(mini_geomx_dataset)
#' # estimate background:
#' mini_geomx_dataset$bg <- derive_GeoMx_background(
#'   norm = mini_geomx_dataset$normalized,
#'   probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'   negnames = "NegProbe"
#' )
#' # run basic decon:
#' res0 <- spatialdecon(
#'   norm = mini_geomx_dataset$normalized,
#'   bg = mini_geomx_dataset$bg,
#'   X = safeTME
#' )
#' # run decon with bells and whistles:
#' res <- spatialdecon(
#'   norm = mini_geomx_dataset$normalized,
#'   bg = mini_geomx_dataset$bg,
#'   X = safeTME,
#'   cellmerges = safeTME.matches,
#'   cell_counts = mini_geomx_dataset$annot$nuclei,
#'   is_pure_tumor = mini_geomx_dataset$annot$AOI.name == "Tumor"
#' )
#' @export
spatialdecon <- function(norm, bg, X = NULL,
                         raw = NULL, wts = NULL,
                         resid_thresh = 3, lower_thresh = 0.5,
                         align_genes = TRUE,
                         is_pure_tumor = NULL, n_tumor_clusters = 10,
                         cell_counts = NULL,
                         cellmerges = NULL,
                         maxit = 1000) {

    #### preliminaries ---------------------------------

    # check formatting:
    if (!is.matrix(norm)) {
        stop("norm should be a matrix")
    }
    if ((length(X) > 0) & (!is.matrix(X))) {
        stop("X should be a matrix")
    }
    if ((length(raw) > 0) & (!is.matrix(raw))) {
        stop("raw must be numeric")
    }
    if ((length(wts) > 0) & (!is.matrix(wts))) {
        stop("wts must be numeric")
    }
    if ((length(cell_counts) > 0) & (!is.numeric(cell_counts))) {
        stop("cell_counts must be numeric")
    }


    if (length(bg) == 1) {
        bg <- matrix(bg, nrow(norm), ncol(norm),
            dimnames = list(rownames(norm), colnames(norm))
        )
    }

    # prep training matrix:
    if (length(X) == 0) {
        X <- SpatialDecon::safeTME
    }
    sharedgenes <- intersect(rownames(norm), rownames(X))
    if (length(sharedgenes) == 0) {
        stop("no shared gene names between norm and X")
    }
    if (length(sharedgenes) < 100) {
        stop(paste0(
            "Only ", length(sharedgenes),
            " genes are shared between norm and X - this may not be enough
                to support accurate deconvolution."
        ))
    }

    # calculate weights based on expected SD of counts
    # wts = replace(norm, TRUE, 1)
    if (length(raw) > 0) {
        weight.by.TIL.resid.sd <-
            length(intersect(colnames(X), colnames(SpatialDecon::safeTME))) > 10
        wts <- deriveWeights(norm,
            raw = raw, error.model = "dsp",
            weight.by.TIL.resid.sd = weight.by.TIL.resid.sd
        )
    }

    #### if pure tumor AOIs are specificed, get tumor expression profile --------
    if (sum(is_pure_tumor) > 0) {

        # derive tumor profiles and merge into X:
        # (derive a separate profile for each tissue)
        X <- mergeTumorIntoX(
            norm = norm,
            bg = bg,
            pure_tumor_ids = is_pure_tumor,
            X = X[sharedgenes, ],
            K = n_tumor_clusters
        )

        sharedgenes <- intersect(rownames(norm), rownames(X))
    }


    #### Run decon  -----------------------------------
    res <- algorithm2(
        Y = norm[sharedgenes, ],
        bg = bg[sharedgenes, ],
        X = X[sharedgenes, ],
        weights = wts[sharedgenes, ],
        maxit = maxit
    )


    #### combine closely-related cell types ------------------------------------

    if (length(cellmerges) > 0) {
        tempconv <- convertCellTypes(
            beta = res$beta,
            matching = cellmerges,
            stat = sum,
            na.rm = FALSE,
            sigma = res$sigmas
        )
        # overwrite original beta with merged beta:
        res$beta.granular <- res$beta
        res$sigma.granular <- res$sigmas
        res$sigmas <- NULL
        res$beta <- tempconv$beta
        res$sigma <- tempconv$sigma
    }


    #### compute p-values -------------------------------------------
    tempbeta <- res$beta
    tempse <- tempp <- tempt <- tempbeta * NA
    for (i in seq_len(ncol(tempse))) {
        tempse[, i] <- suppressWarnings(sqrt(diag(res$sigma[, , i])))
    }
    tempt <- (tempbeta / tempse)
    tempp <- 2 * (1 - stats::pt(tempt, df = length(sharedgenes) - ncol(X) - 1))
    res$p <- tempp
    res$t <- tempt
    res$se <- tempse


    #### rescale abundance estimates --------------------------------
    # (to proportions of total, proportions of immune, cell counts)

    # proportions:
    res$prop_of_all <- sweep(res$beta, 2, colSums(res$beta), "/")
    nontumorcellnames <- rownames(res$beta)[!grepl("tumor", rownames(res$beta))]
    res$prop_of_nontumor <- sweep(
        res$beta[nontumorcellnames, ], 2,
        colSums(res$beta[nontumorcellnames, ]), "/"
    )

    # on scale of cell counts:
    if (length(cell_counts) > 0) {
        res$cell.counts <- convertCellScoresToCounts(
            beta = res$beta,
            nuclei.counts = cell_counts,
            omit.tumor = TRUE
        )
        if (exists("res$beta.granular") > 0) {
            res$cell.counts.granular <- convertCellScoresToCounts(
                beta = res$beta.granular,
                nuclei.counts = cell_counts,
                omit.tumor = TRUE
            )
        }
    }

    # add other pertinent info to res:
    res$X <- X[rownames(res$resids), ]
    return(res)
}
