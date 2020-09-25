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


#' Estimate a tumor-specific profile and merge it with the pre-specified cell
#'  profile matrix (X)
#'
#' Given the input of "tumor-only" AOI's, estimates an collection of
#'  tumor-specific
#' expression profiles and merges them with the immune cell expression
#' training matrix.
#' The process:
#' \enumerate{
#' \item log2/normalized data from tumor-only AOIs is clustered with hclust,
#' and cutree() is used to define clusters.
#' \item 2. Each cluster's geomean profile is merged into the immune cell
#' profile matrix.
#' }
#'
#' @param norm matrix of normalized data
#' @param bg matrix of expected background, on the scale of norm.
#' @param pure_tumor_ids Vector identifying columns of norm that are pure tumor.
#'  Can be indices, logicals or column names.
#' @param X The training matrix
#' @param K the number of clusters to fit
#' @return an updated X matrix with new columns, "tumor.1", "tumor.2", ...
#' @examples
#' data(mini_geomx_dataset)
#' mini_geomx_dataset$bg <- derive_GeoMx_background(
#'   norm = mini_geomx_dataset$normalized,
#'   probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'   negnames = "NegProbe"
#' )
#' safeTME.with.tumor <- mergeTumorIntoX(
#'   norm = mini_geomx_dataset$norm,
#'   bg = mini_geomx_dataset$bg,
#'   pure_tumor_ids = mini_geomx_dataset$annot$AOI.name == "Tumor",
#'   X = safeTME,
#'   K = 3
#' )
#' @export
mergeTumorIntoX <- function(norm, bg, pure_tumor_ids, X, K = 10) {

    # round up 0 values in norm:
    min.nonzero <- min(norm[norm > 0], na.rm = TRUE)
    norm <- pmax(norm, min.nonzero)

    # subset data to only the pure tumor IDs:
    norm <- norm[, pure_tumor_ids, drop = FALSE]
    bg <- bg[, pure_tumor_ids, drop = FALSE]

    # bg-subtract:
    norm <- pmax(norm - bg, min(norm) / 20)

    # fix K if too big:
    if (ncol(norm) < K) {
        K <- ncol(norm)
    }

    # case 1: want to use every column in norm as a separate profile:
    #  (includes case of just one column in norm)
    if (K == ncol(norm)) {
        tumorX <- norm
    }

    # case 2: if many tumor AOIs, get profiles for K clusters of data:
    if (K < ncol(norm)) {
        # cluster and cut:
        h <- stats::hclust(stats::dist(t(log2(norm))))
        cut <- stats::cutree(h, k = K)
        # get clusters' geomean profiles:
        tumorX <- c()
        for (cid in unique(cut)) {
            tumorX <- cbind(
                tumorX,
                exp(rowMeans(log(norm[, cut == cid, drop = FALSE])))
            )
        }
        colnames(tumorX) <- paste0("tumor.", seq_len(ncol(tumorX)))
    }

    # align tumorX with X:
    sharedgenes <- intersect(rownames(tumorX), rownames(X))
    tumorX <- tumorX[sharedgenes, ]
    X <- X[sharedgenes, ]

    # rescale tumor X:
    meanq90 <- max(mean(apply(X, 2, stats::quantile, 0.9)), 1e-3)
    tumorq90s <- apply(tumorX, 2, stats::quantile, 0.9)
    tumorX <- sweep(tumorX, 2, tumorq90s, "/") * meanq90

    # merge:
    out <- cbind(X, tumorX)
    return(out)
}
