#' Run collapseCellTypes
#' 
#' Runs collapseCellTypes from an S4 object
#' @param object An S4 object such as a GeoMxSet object
#' @param ... Arguments passed to collapseCellTypes
#' @return A reshaped deconvolution result object
#' @examples
#' library(GeomxTools)
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data", package = "GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' 
#' demoData <- shiftCountsOne(demoData)
#' target_demoData <- aggregateCounts(demoData)
#' 
#' target_demoData <- normalize(target_demoData, "quant")
#'                 
#' # run basic decon:
#' res0 <- runspatialdecon(object = target_demoData,
#'                         norm_elt = "exprs_norm",
#'                         raw_elt = "exprs")
#'
#' # run reverse decon:
#' target_demoData <- runReverseDecon(object = target_demoData,
#'                                    norm_elt = "exprs_norm",
#'                                    beta = pData(res0)$beta)
#'

setGeneric("runCollapseCellTypes", signature = "object",
           function(object, ...) standardGeneric("runCollapseCellTypes"))

#' Collapse related cell types within a deconvolution result
#'
#' Given the input of an SpatialDecon result output and a list of which cell
#' types to combine,
#'  returns a reshaped deconvolution result object with the specified cell
#'  types merged.
#' @param object A NanostringGeomxSet object returned by the SpatialDecon algorithm
#' @param matching A list object holding the mapping from beta's cell names to
#' official cell names.
#'  See str(safeTME.matches) for an example.
#' @return A reshaped deconvolution result object
#' @importFrom stats pnorm
#' @importFrom Biobase assayData
#' @examples
#' library(GeomxTools)
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data", package = "GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' 
#' demoData <- shiftCountsOne(demoData)
#' target_demoData <- aggregateCounts(demoData)
#' 
#' target_demoData <- normalize(target_demoData, "quant")
#'                 
#' # run basic decon:
#' res0 <- runspatialdecon(object = target_demoData,
#'                         norm_elt = "exprs_norm",
#'                         raw_elt = "exprs")
#'
#' # run reverse decon:
#' target_demoData <- runReverseDecon(object = target_demoData,
#'                                    norm_elt = "exprs_norm",
#'                                    beta = pData(res0)$beta)
#'
#' @export
#'
setMethod("runCollapseCellTypes",  "NanoStringGeoMxSet",
          function(object, matching = NULL){
              
              if(is.null(matching)){
                  stop("matching must be set")
              }
              
              fit <- list(beta=t(Biobase::pData(object)$beta),
                          p=t(Biobase::pData(object)$p),
                          t=t(Biobase::pData(object)$t),
                          se=t(Biobase::pData(object)$se),
                          prop_of_all=t(Biobase::pData(object)$prop_of_all),
                          prop_of_nontumor=t(Biobase::pData(object)$prop_of_nontumor),
                          sigma=Biobase::pData(object)$sigmas)
              
              temp <- collapseCellTypes(fit = fit, 
                                        matching = matching)
              
              Biobase::pData(object)$beta <- t(temp$beta)
              Biobase::pData(object)$p <- t(temp$p)
              Biobase::pData(object)$t <- t(temp$t)
              Biobase::pData(object)$se <- t(temp$se)
              Biobase::pData(object)$prop_of_all <- t(temp$prop_of_all)
              Biobase::pData(object)$prop_of_nontumor <- t(temp$prop_of_nontumor)
              Biobase::pData(object)$sigmas <- temp$sigma
              
              return(object)
          })
