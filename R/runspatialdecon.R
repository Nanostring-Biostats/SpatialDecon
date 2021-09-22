#' Run spatialdecon
#' 
#' Runs spatialdecon from an S4 object
#' @param object An S4 object such as a Seurat object or a GeoMxSet object
#' @param ... Arguments passed to spatialdecon
setGeneric("runspatialdecon", signature = "object",
           function(object, ...) standardGeneric("runspatialdecon"))


#' Run spatialdecon on a Seurat object
#' 
#' A wrapper for applying spatialdecon to the Spatial data element in a Seurat object. 
#' Unlike spatialdecon, which expects a normalized data matrix, this function operates 
#' on raw counts. Scaling for total cells 
#' @param object A seurat object. Must include a "Spatial" element in the "assays" slot.
#' @param X Cell profile matrix. If NULL, the safeTME matrix is used.
#' @param bg Expected background counts. Either a scalar applied equally to 
#'  all points in the count matrix, or a matrix with the same dimensions 
#'  as the count matrix in GetAssayData(object, assay = "Spatial").
#'  Recommended to use a small non-zero value, default of 0.1.
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
#' @return if not given cellmerges and cell_counts, a list
#' including the following items:
#'  \itemize{
#'    \item beta: matrix of cell abundance estimates, cells in rows and
#'                observations in columns
#'    \item p: matrix of p-values for H0: beta == 0
#'    \item t: matrix of t-statistics for H0: beta == 0
#'    \item se: matrix of standard errors of beta values
#'    \item prop_of_all: rescaling of beta to sum to 1 in each observation
#'    \item prop_of_nontumor: rescaling of beta to sum to 1 in each observation,
#'                            excluding tumor abundance estimates
#'    \item yhat: a matrix of fitted values
#'    \item resids: a matrix of residuals from the model fit.
#'                  (log2(pmax(y, lower_thresh)) - log2(pmax(xb, lower_thresh))).
#'    \item X: the cell profile matrix used in the decon fit.
#'    \item sigmas: covariance matrices of each observation's beta estimates
#'}
#'
#'
#' if given cellmerges, the list will additionally include
#' the following items
#'  \itemize{
#'    \item beta.granular: cell abundances prior to combining closely-related
#'                         cell types
#'    \item sigma.granular: sigmas prior to combining closely-related cell types
#'}
#'
#' if given cell_counts, the list will additionally include
#'  the following items
#'  \itemize{
#'    \item cell.counts: beta rescaled to estimate cell numbers, based on
#'                       prop_of_all and nuclei count
#'}
#'
#' if given both cellmerges and cell_counts, the list will
#'  additionally include the following items
#'  \itemize{
#'    \item cell.counts.granular: cell.counts prior to combining closely-related
#'                                cell types
#'}
#' @importFrom SeuratObject GetAssayData
#' @importClassesFrom GeomxTools NanoStringGeoMxSet
#' @export
#' @examples
#' # get mouse brain data:
#' library(SeuratData)
#' brain <- SeuratData::LoadData("stxBrain", type = "anterior1")
#' # get cell profile matrix:
#' ref <- download_profile_matrix("Mouse_Brain")
#' res <- runspatialdecon(brain, X = ref)
#' str(res)
setMethod("runspatialdecon", "Seurat", function(
  object,
  X = NULL,
  bg = 0.1,  
  wts = NULL,
  resid_thresh = 3, lower_thresh = 0.5,
  align_genes = TRUE,
  is_pure_tumor = NULL, 
  n_tumor_clusters = 10,
  cell_counts = NULL,
  cellmerges = NULL,
  maxit = 1000) {
  
  # check that it's a spatial assay:
  if (!is.element("Spatial", names(object@assays))) {
    stop("Expecting \'Spatial\' element in assays slot")
  }
  
  # prep components:
  raw <- as.matrix(SeuratObject::GetAssayData(object, assay = "Spatial"))
  # make bg a matrix if only a scalar was specified:
  if (length(bg) == 1) {
    bg <- 0 * raw + bg
  }
  
  # run spatialdecon:
  res <- spatialdecon(norm = raw, 
                      bg = bg, 
                      X = X,
                      raw = raw, 
                      wts = wts,
                      resid_thresh = resid_thresh, 
                      lower_thresh = lower_thresh,
                      align_genes = align_genes,
                      is_pure_tumor = is_pure_tumor, 
                      n_tumor_clusters = n_tumor_clusters,
                      cell_counts = cell_counts,
                      cellmerges = cellmerges,
                      maxit = maxit) 
  
  return(res)
})

#' Run spatialdecon on a NanostringGeomxSet object
#'
#' A wrapper for applying spatialdecon to the Spatial data element in a NanostringGeomxSet object.
#'
#' @param object A NanostringGeomxSet object.
#' @param norm_elt normalized data element in assayData
#' @param raw_elt raw data element in assayData
#' @param X Cell profile matrix. If NULL, the safeTME matrix is used.
#' @param bg Expected background counts. Either a scalar applied equally to
#'  all points in the count matrix, or a matrix with the same dimensions
#'  as the count matrix in assayDataElement(object , elt = norm_elt).
#'  Recommended to use a small non-zero value, default of 0.1.
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

#' @importFrom Biobase assayData
#' @importClassesFrom GeomxTools NanoStringGeoMxSet
#'
#' @return if not given cellmerges and cell_counts, a valid GeoMx S4 object
#' including the following items
#'  \itemize{
#'    \item In pData
#'    \itemize{
#'       \item beta: matrix of cell abundance estimates, cells in rows and
#'                observations in columns
#'       \item p: matrix of p-values for H0: beta == 0
#'       \item t: matrix of t-statistics for H0: beta == 0
#'       \item se: matrix of standard errors of beta values
#'       \item prop_of_all: rescaling of beta to sum to 1 in each observation
#'       \item prop_of_nontumor: rescaling of beta to sum to 1 in each observation,
#'                            excluding tumor abundance estimates
#'       \item sigmas: covariance matrices of each observation's beta estimates
#'      }
#'    \item In assayData
#'    \itemize{
#'       \item yhat: a matrix of fitted values
#'       \item resids: a matrix of residuals from the model fit.
#'                  (log2(pmax(y, lower_thresh)) - log2(pmax(xb, lower_thresh))).
#'    }
#'    \item In experimentData
#'    \itemize{
#'       \item SpatialDeconMatrix: the cell profile matrix used in the decon fit.
#'    }
#'}
#'
#'
#' if given cellmerges, the valid GeoMx S4 object will additionally include
#' the following items
#'  \itemize{
#'    \item In pData
#'    \itemize{
#'       \item beta.granular: cell abundances prior to combining closely-related
#'                         cell types
#'       \item sigma.granular: sigmas prior to combining closely-related cell types
#'    }
#'}
#'
#' if given cell_counts, the valid GeoMx S4 object will additionally include
#'  the following items
#'  \itemize{
#'  \item In pData
#'    \itemize{
#'       \item cell.counts: beta rescaled to estimate cell numbers, based on
#'                       prop_of_all and nuclei count
#'    }
#'}
#'
#' if given both cellmerges and cell_counts, the valid GeoMx S4 object will
#'  additionally include the following items
#'  \itemize{
#'  \item In pData
#'    \itemize{
#'       \item cell.counts.granular: cell.counts prior to combining closely-related
#'                                cell types
#'    }
#'}
#'
#' @examples
#'
#' library(GeomxTools)
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data", package = "GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' 
#' demoData <- shiftCountsOne(demoData)
#' demoData <- aggregateCounts(demoData)
#' demoData <- runspatialdecon(object = demoData, 
#'                             norm_elt = "exprs_norm",
#'                             raw_elt = "exprs")
#'
#' @export
#'

setMethod("runspatialdecon", "NanoStringGeoMxSet",
          function(object, X = NULL, norm_elt = NULL, raw_elt = NULL,
                   wts = NULL,
                   resid_thresh = 3, lower_thresh = 0.5,
                   align_genes = TRUE,
                   is_pure_tumor = NULL, n_tumor_clusters = 10,
                   cell_counts = NULL,
                   cellmerges = NULL,
                   maxit = 1000){
            
            # check that norm and raw data exists:
            if(is.null(norm_elt)){
              stop("norm_elt must be set")
            }
            if(is.null(raw_elt)){
              stop("raw_elt must be set")
            }
            if (!is.element(norm_elt, names(object@assayData))) {
              stop(paste(norm_elt, "is not an element in assaysData slot"))
            }
            if (!is.element(raw_elt, names(object@assayData))) {
              stop(paste(raw_elt, "is not an element in assaysData slot"))
            }
          
            # prep components:
            norm <- as.matrix(Biobase::assayDataElement(object , elt = norm_elt))
            raw <- as.matrix(Biobase::assayDataElement(object , elt = raw_elt))
            
            # estimate background
            bg <- derive_GeoMx_background(norm = norm,
                                          # access the probe pool information from the feature metadata
                                          probepool = Biobase::fData(object)$Module,
                                          # access the names of the negative control probes
                                          negnames = Biobase::fData(object)$TargetName[Biobase::fData(object)$Negative])
            
            # run spatialdecon:
            result <- spatialdecon(norm = norm,
                                   bg = bg,
                                   X = X,
                                   raw = raw,
                                   wts = wts,
                                   resid_thresh = resid_thresh, 
                                   lower_thresh = lower_thresh,
                                   align_genes = align_genes,
                                   is_pure_tumor = is_pure_tumor, 
                                   n_tumor_clusters = n_tumor_clusters,
                                   cell_counts = cell_counts,
                                   cellmerges = cellmerges,
                                   maxit = maxit)
            
            
            # append results to the object
            # beta
            Biobase::pData(object)[["beta"]] <- matrix(NA, nrow = ncol(object), ncol = nrow(result$beta),
                                                       dimnames = list(Biobase::sampleNames(object), rownames(result$beta)))
            Biobase::pData(object)[["beta"]][colnames(result$beta), ] <- t(result$beta)
            
            
            # p
            Biobase::pData(object)[["p"]] <- matrix(NA, nrow = ncol(object), ncol = nrow(result$p),
                                                    dimnames = list(Biobase::sampleNames(object), rownames(result$p)))
            Biobase::pData(object)[["p"]][colnames(result$p), ] <- t(result$p)
            
            
            # t
            Biobase::pData(object)[["t"]] <- matrix(NA, nrow = ncol(object), ncol = nrow(result$t),
                                                    dimnames = list(Biobase::sampleNames(object), rownames(result$t)))
            Biobase::pData(object)[["t"]][colnames(result$t), ] <- t(result$t)
            
            
            # se
            Biobase::pData(object)[["se"]] <- matrix(NA, nrow = ncol(object), ncol = nrow(result$se),
                                                     dimnames = list(Biobase::sampleNames(object), rownames(result$se)))
            Biobase::pData(object)[["se"]][colnames(result$se), ] <- t(result$se)
            
            
            # prop_of_all
            Biobase::pData(object)[["prop_of_all"]] <- matrix(NA, nrow = ncol(object), ncol = nrow(result$prop_of_all),
                                                              dimnames = list(Biobase::sampleNames(object), rownames(result$prop_of_all)))
            Biobase::pData(object)[["prop_of_all"]][colnames(result$prop_of_all), ] <- t(result$prop_of_all)
            
            
            # prop_of_nontumor
            Biobase::pData(object)[["prop_of_nontumor"]] <- matrix(NA, nrow = ncol(object), ncol = nrow(result$prop_of_nontumor),
                                                                   dimnames = list(Biobase::sampleNames(object), rownames(result$prop_of_nontumor)))
            Biobase::pData(object)[["prop_of_nontumor"]][colnames(result$prop_of_nontumor), ] <- t(result$prop_of_nontumor)
            
            
            # add yhat matrix
            yhat_mat <- matrix(NA, nrow = nrow(object), ncol = ncol(object),
                               dimnames = dimnames(object))
            yhat_mat[rownames(result$yhat), colnames(result$yhat)] <- result$yhat
            Biobase::assayDataElement(object, "yhat") <- yhat_mat
            
            
            # add resids matrix
            resids_mat <- matrix(NA, nrow = nrow(object), ncol = ncol(object),
                                 dimnames = dimnames(object))
            resids_mat[rownames(result$resids), colnames(result$resids)] <- result$resids
            Biobase::assayDataElement(object, "resids") <- resids_mat
            
            
            # add X matrix
            Biobase::fData(object)[["X"]] <- matrix(NA, nrow = nrow(object), ncol = ncol(result$X),
                                                    dimnames = list(Biobase::featureNames(object), colnames(result$X)))
            Biobase::fData(object)[["X"]][rownames(result$X), ] <- result$X
            
            
            # transpose 3-d array
            trans <- function(array){
              NewArray <- array(NA, dim = c(dim(array)[3], dim(array)[2], dim(array)[1]))
              
              for(or in(1:dim(array)[1])){
                for(oc in(1:dim(array)[2])){
                  NewArray[,oc,or] <- array[or,oc,]
                }
              }
              array <- NewArray
            }
            
            if (is.null(cellmerges)) {
              # add sigmas matrix
              sigmas_mat <- trans(result$sigmas)
              sigmas_dimname <- list(Biobase::sampleNames(object),
                                     paste0("var", seq_len((dim(result$sigmas)[1]))),
                                     paste0("var", seq_len((dim(result$sigmas)[2]))))
              dimnames(sigmas_mat) <- sigmas_dimname
              
              Biobase::pData(object)[["sigmas"]] <- array(NA, dim = c(ncol(object),
                                                                      dim(result$sigmas)[1],
                                                                      dim(result$sigmas)[2]),
                                                          dimnames = sigmas_dimname)
              Biobase::pData(object)[["sigmas"]][unlist(dimnames(sigmas_mat)[1]),, ] <- sigmas_mat
            } else {
              # add sigma matrix
              sigma_mat <- trans(result$sigma)
              sigma_dimname <- list(Biobase::sampleNames(object),
                                    paste0("var", seq_len((dim(result$sigma)[1]))),
                                    paste0("var", seq_len((dim(result$sigma)[2]))))
              dimnames(sigma_mat) <- sigma_dimname
              
              Biobase::pData(object)[["sigma"]] <- array(NA, dim = c(ncol(object),
                                                                     dim(result$sigma)[1],
                                                                     dim(result$sigma)[2]),
                                                         dimnames = sigma_dimname)
              Biobase::pData(object)[["sigma"]][unlist(dimnames(sigma_mat)[1]),, ] <- sigma_mat
              
              
              # add beta.granular matrix
              Biobase::pData(object)[["beta.granular"]] <- matrix(NA, nrow = ncol(object),
                                                                  ncol = nrow(result$beta.granular),
                                                                  dimnames = list(Biobase::sampleNames(object),
                                                                                  rownames(result$beta.granular)))
              Biobase::pData(object)[["beta.granular"]][colnames(result$beta.granular), ] <- t(result$beta.granular)
              
              
              # add sigma.granular matrix
              sigma.granular_mat <- trans(result$sigma.granular)
              sigma.granular_dimname <- list(Biobase::sampleNames(object),
                                             paste0("var", seq_len((dim(result$sigma.granular)[1]))),
                                             paste0("var", seq_len((dim(result$sigma.granular)[2]))))
              dimnames(sigma.granular_mat) <- sigma.granular_dimname
              
              Biobase::pData(object)[["sigma.granular"]] <- array(NA, dim = c(ncol(object),
                                                                              dim(result$sigma.granular)[1],
                                                                              dim(result$sigma.granular)[2]),
                                                                  dimnames = sigma.granular_dimname)
              Biobase::pData(object)[["sigma.granular"]][unlist(dimnames(sigma.granular_mat)[1]),, ] <- sigma.granular_mat
            }
            
            
            if (length(cell_counts) > 0) {
              if (length(result$cell.counts) == 2){
                # create cell.counts list
                Biobase::pData(object)[["cell.counts"]] <- data.frame(matrix(NA,
                                                                             nrow = ncol(object),
                                                                             ncol = 2,
                                                                             dimnames = list(Biobase::sampleNames(object),
                                                                                             c("cells.per.100",                                                                                           "cells.counts"))))
                # add cells.per.100 matrix
                Biobase::pData(object)[["cell.counts"]][["cells.per.100"]] <- matrix(NA, nrow = ncol(object),
                                                                                     ncol = nrow(result$cell.counts$cells.per.100),
                                                                                     dimnames = list(Biobase::sampleNames(object),
                                                                                                     rownames(result$cell.counts$cells.per.100)))
                Biobase::pData(object)[["cell.counts"]][["cells.per.100"]][colnames(result$cell.counts$cells.per.100), ] <- t(result$cell.counts$cells.per.100)
                
                
                # add cells.counts matrix
                Biobase::pData(object)[["cell.counts"]][["cells.counts"]] <- matrix(NA, nrow = ncol(object),
                                                                                    ncol = nrow(result$cell.counts$cells.counts),
                                                                                    dimnames = list(Biobase::sampleNames(object),
                                                                                                    rownames(result$cell.counts$cells.counts)))
                Biobase::pData(object)[["cell.counts"]][["cells.counts"]][colnames(result$cell.counts$cells.counts), ] <- t(result$cell.counts$cells.counts)
              } else {
                # create cell.counts list
                Biobase::pData(object)[["cell.counts"]] <- data.frame(matrix(NA,
                                                                             nrow = ncol(object),
                                                                             ncol = 1,
                                                                             dimnames = list(Biobase::sampleNames(object),
                                                                                             "cells.per.100")))
                
                # add cells.per.100 matrix
                Biobase::pData(object)[["cell.counts"]][["cells.per.100"]] <- matrix(NA, nrow = ncol(object),
                                                                                     ncol = nrow(result$cell.counts$cells.per.100),
                                                                                     dimnames = list(Biobase::sampleNames(object),
                                                                                                     rownames(result$cell.counts$cells.per.100)))
                Biobase::pData(object)[["cell.counts"]][["cells.per.100"]][colnames(result$cell.counts$cells.per.100), ] <- t(result$cell.counts$cells.per.100)
              }
            }
            
            
            if (length(result$cell.counts.granular) > 0) {
              
              if(length(result$cell.counts.granular) == 2){
                # create cell.counts.granular list
                Biobase::pData(object)[["cell.counts.granular"]] <- data.frame(matrix(NA,
                                                                                      nrow = ncol(object),
                                                                                      ncol = 2,
                                                                                      dimnames = list(Biobase::sampleNames(object),
                                                                                                      c("cells.per.100",
                                                                                                        "cells.counts"))))
                # add cells.per.100 matrix
                Biobase::pData(object)[["cell.counts.granular"]][["cells.per.100"]] <- matrix(NA, nrow = ncol(object),
                                                                                              ncol = nrow(result$cell.counts.granular$cells.per.100),
                                                                                              dimnames = list(Biobase::sampleNames(object),
                                                                                                              rownames(result$cell.counts.granular$cells.per.100)))
                Biobase::pData(object)[["cell.counts.granular"]][["cells.per.100"]][colnames(result$cell.counts.granular$cells.per.100), ] <- t(result$cell.counts.granular$cells.per.100)
                
                
                # add cells.count matrix
                Biobase::pData(object)[["cell.counts.granular"]][["cells.counts"]] <- matrix(NA, nrow = ncol(object),
                                                                                             ncol = nrow(result$cell.counts.granular$cells.counts),
                                                                                             dimnames = list(Biobase::sampleNames(object),
                                                                                                             rownames(result$cell.counts.granular$cells.counts)))
                Biobase::pData(object)[["cell.counts.granular"]][["cells.counts"]][colnames(result$cell.counts.granular$cells.per.100), ] <- t(result$cell.counts.granular$cells.counts)
              } else {
                # create cell.counts.granular list
                Biobase::pData(object)[["cell.counts.granular"]] <- data.frame(matrix(NA,
                                                                                      nrow = ncol(object),
                                                                                      ncol = 1,
                                                                                      dimnames = list(Biobase::sampleNames(object),
                                                                                                      "cells.per.100")))
                # add cells.per.100 matrix
                Biobase::pData(object)[["cell.counts.granular"]][["cells.per.100"]] <- matrix(NA, nrow = ncol(object),
                                                                                              ncol = nrow(result$cell.counts.granular$cells.per.100),
                                                                                              dimnames = list(Biobase::sampleNames(object),
                                                                                                              rownames(result$cell.counts.granular$cells.per.100)))
                Biobase::pData(object)[["cell.counts.granular"]][["cells.per.100"]][colnames(result$cell.counts.granular$cells.per.100), ] <- t(result$cell.counts.granular$cells.per.100)
              }
            }
            
            object@experimentData@other$SpatialDeconMatrix <- result$X
            
            return(object)
            
          })


