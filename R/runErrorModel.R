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



#' Apply error model to estimate technical SD from raw counts
#'
#' Based on raw counts, uses past data to estimate each raw count's log-scale
#' SD from technical noise.
#' Specifies different error models for different platforms.
#'
#' @param counts vector or matrix of raw counts
#' @param platform String specifying which platform was used to create
#' "rawCounts". Default to "dsp".
#'  Other options include "ncounter", "rsem" and "quantile".
#' @return a matrix of log2-scale SDs
runErrorModel <- function(counts, platform = "general") {
    if (platform == "ncounter") {
        sds <- counts * 0 + 0.1
        sds <- replace(sds, counts < 200, 0.2)
        sds <- replace(sds, counts < 100, 0.3)
        sds <- replace(sds, counts < 75, 0.4)
        sds <- replace(sds, counts < 50, 0.5)
        sds <- replace(sds, counts < 40, 0.7)
        sds <- replace(sds, counts < 30, 1)
        sds <- replace(sds, counts < 20, 3)
    }


    if (platform == "rsem") {
        sds <- counts * 0 + 0.5930982
        sds <- replace(sds, log2(counts) < 9.5, 0.6458475)
        sds <- replace(sds, log2(counts) < 8.5, 0.7847597)
        sds <- replace(sds, log2(counts) < 7.5, 1.0576471)
        sds <- replace(sds, log2(counts) < 6.5, 1.2990917)
        sds <- replace(sds, log2(counts) < 5.5, 1.5061735)
        sds <- replace(sds, log2(counts) < 4.5, 1.6930872)
        sds <- replace(sds, log2(counts) < 3.5, 1.7894239)
    }

    if (platform == "dsp") {
        predictsd.dsp <- function(rawcounts) {
            m <- log2(pmax(rawcounts, 1e-3))
            meanvec <- c(-1e-6, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, Inf)
            sdvec <- c(
                1.5, 1.383, 1.191, 0.800, 0.48, 0.301, 0.301,
                0.301, 0.263, 0.235, 0.235
            )

            s <- replace(m, TRUE, sdvec[1])
            for (i in seq_len(length(meanvec) - 1)) {
                s <- replace(s, m >= meanvec[i], sdvec[i + 1])
            }
            return(s)
        }

        if (is.vector(counts)) {
            sds <- vapply(
                X = counts,
                FUN = predictsd.dsp,
                FUN.VALUE = numeric(length(counts))
            )
        }
        if (is.matrix(counts)) {
            sds <- predictsd.dsp(counts)
        }
    }


    if (platform == "quantile") {
        if (is.vector(counts)) {
            quantile <- rank(counts) / length(counts)
        }
        if (is.matrix(counts)) {
            quantile <- matrix(rank(counts), nrow(counts)) / length(counts)
        }

        sds <- quantile * 0 + 0.1
        sds <- replace(sds, quantile < 0.2, 0.2)
        sds <- replace(sds, quantile < 0.15, 0.3)
        sds <- replace(sds, quantile < 0.1, 0.4)
        sds <- replace(sds, quantile < 0.05, 0.5)
        sds <- replace(sds, quantile < 0.01, 1)
    }
    return(sds)
}
