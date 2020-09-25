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


#' Convert abundance measurements to cell counts
#'
#' Converts cell abundance scores to cell counts, under the assumption that the
#'  observation with
#'  the greatest sum of cell abundance scores is 100% modelled cells.
#' @param beta Matrix of cell abundance scores, with cells in rows and
#'  observations in columns. The
#'  assumption is that this matrix is from well-normalized data.
#' @param nuclei.counts Optional. A vector of total nuclei counts. If provided,
#' the function will
#'  output not only cells.per.100 but also total cells.
#' @param omit.tumor Logical. If FALSE, any rows of beta with "tumor" in their
#' name will be omitted.
#' @return A list with two elements, each a rescaled version of beta.
#' cells.per.100 gives estimated
#'  percents of total, and cell.counts is cells.per.100 * nuclei.counts.
convertCellScoresToCounts <- function(beta, nuclei.counts = NULL,
                                      omit.tumor = FALSE) {
    # strip tumor rows if called for:
    if (omit.tumor) {
        beta <- beta[!grepl("tumor", rownames(beta)), , drop = FALSE]
    }

    # calc max abundance scores:
    max.total.abundance <- max(colSums(beta))

    # calculate rescaled scores:
    out <- list()
    out$cells.per.100 <- beta / max.total.abundance * 100
    if (length(nuclei.counts) == ncol(beta)) {
        out$cell.counts <- sweep(out$cells.per.100, 2, nuclei.counts, "*") / 100
    }
    return(out)
}
