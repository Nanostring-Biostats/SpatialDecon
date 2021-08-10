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
#' @importFrom SeuratObject GetAssayData
#' @export
#' @example 
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
                      resid_thresh = resid_thresh, lower_thresh = lower_thresh,
                      align_genes = align_genes,
                      is_pure_tumor = is_pure_tumor, n_tumor_clusters = n_tumor_clusters,
                      cell_counts = cell_counts,
                      cellmerges = cellmerges,
                      maxit = maxit) 
  
  return(res)
})


