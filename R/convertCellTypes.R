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


#' Convert results from an arbitrary training matrix to our standardized cell
#' types.
#'
#' Takes betas from a decon run, using any cell type names whatsoever, and maps
#'  them back to the
#' "official" cell types we use. Allows for multiple rows of beta to map to the
#'  same official cell type,
#' in which case those rows will be added up.
#'
#' @param beta K * n matrix of estimated beta values (cell type abundances)
#' @param matching A list object holding the mapping from beta's cell names to
#'  official cell names.
#'  See str(safeTME.matches) for an example.
#' @param stat The function used to combine related cell types. Defaults to sum.
#' @param na.rm Whether to ignore NAs while computing stat
#' @param sigma A list of covariance matrices of beta estimates, in the format
#'  output by spatialdecon.
#' @return A list with two elements:
#' \itemize{
#' \item beta: a matrix of cell abundances, with specified cell types added
#' together
#' \item sigma: an array of covariance matrices for each observation's beta
#' vector
#' }
convertCellTypes <- function(beta, matching, stat = sum,
                             na.rm = FALSE, sigma = NULL) {
    # format matching list as a matrix to take a linear combination of beta:
    A <- matrix(0, length(matching), nrow(beta),
        dimnames = list(names(matching), rownames(beta))
    )
    for (name in names(matching)) {
        cellnames <- matching[[name]]
        A[name, cellnames] <- 1
    }

    # apply A transformation to beta:
    beta2 <- A %*% beta

    # if Sigma provided, get vcov of beta2:
    if (length(sigma) > 0) {
        if (length(dim(sigma)) == 2) {
            sigma2 <- A %*% sigma %*% t(A)
        }
        if (length(dim(sigma)) == 3) {
            sigma2 <- array(NA,
                dim = c(nrow(A), nrow(A), dim(sigma)[3]),
                dimnames = list(rownames(A), rownames(A), dimnames(sigma)[[3]])
            )
            for (i in seq_len(dim(sigma)[3])) {
                sigma2[, , i] <- A %*% sigma[, , i] %*% t(A)
            }
        }
    }

    # if no Sigma, just return transformed beta:
    if (length(sigma) == 0) {
        return(beta2)
    }
    # if there is a sigma, return beta and the sigma:
    if (length(sigma) > 0) {
        out <- list(beta = beta2, sigma = sigma2)
        return(out)
    }
}
