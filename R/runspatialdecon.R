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
#' @return Appends spatialdecon results to the object as \code{object@misc$spatialdecon}
#' @importFrom SeuratObject GetAssayData
#' @export
#' @example 
#' # get mouse brain data:
#' library(SeuratData)
#' brain <- SeuratData::LoadData("stxBrain", type = "anterior1")
#' # get cell profile matrix:
#' ref <- download_profile_matrix("Mouse_Brain")
#' brain <- runspatialdecon(brain, X = ref)
#' str(brain@misc$spatialdecon)
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
                      wts = NULL,
                      resid_thresh = resid_thresh, lower_thresh = lower_thresh,
                      align_genes = align_genes,
                      is_pure_tumor = is_pure_tumor, n_tumor_clusters = n_tumor_clusters,
                      cell_counts = cell_counts,
                      cellmerges = cellmerges,
                      maxit = maxit) 
  
  # append results to seurat object:
  object@misc$spatialdecon <- res
  return(object)
})




#' Run spatialdecon on a GeoMxSet object
#' 
#' A wrapper for applying spatialdecon to a GeoMxSet object. 
#' @param object A GeoMxSet object. 
#' @param X Cell profile matrix. If NULL, the safeTME matrix is used.
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
#' @return Appends spatialdecon results to the object as \code{object@___________}
#' @importFrom ________ __________  
#' @export
#' @example 
#' ________________
#'
#'
#'
setMethod("runspatialdecon", "GeoMxSet", function() {
  
  # check that the object is in proper geomxset format:
  
  
  # estimate background
  bg <- derive_GeoMx_background(norm = _____,
                                probepool = _____,
                                negnames = _____)
  
  # run spatialdecon:
  res <- spatialdecon(norm = ______, 
                      bg = bg, 
                      X = X,
                      raw = ______, 
                      wts = NULL,
                      resid_thresh = resid_thresh, lower_thresh = lower_thresh,
                      align_genes = align_genes,
                      is_pure_tumor = is_pure_tumor, n_tumor_clusters = n_tumor_clusters,
                      cell_counts = cell_counts,
                      cellmerges = cellmerges,
                      maxit = maxit) 
  
  # append to object:
  __________
})