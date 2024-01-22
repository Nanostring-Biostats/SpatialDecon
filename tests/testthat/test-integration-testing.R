context("spatialdecon")

library(testthat)


# strategy here:
# spatialdecon calls almost every function in the package
# we'll run it, then dissect its results to check that every piece of the
# puzzle worked we'll also test individual functions as a sanity check

#### load test data ----------------------
rm(list = ls())
load("testdata.RData")
load("expectedtestresults.RData")
data('safeTME')
data("safeTME.matches")
sharedgenes <- intersect(rownames(safeTME), rownames(snr))



#### test subsidiary functions individually -----------------------------------

test_that("deriveWeights is as expected", {
  wts <- deriveWeights(
    norm = snr,
    raw = raw,
    error.model = "dsp"
  )
  expect_true(all(abs(wts.test - wts) < 1e-3))
})

# merging tumor profiles:
test_that("mergeTumorIntoX is as expected", {
  set.seed(0)
  mergedX <- mergeTumorIntoX(
    norm = snr[sharedgenes, ],
    bg = replace(snr, TRUE, 1)[sharedgenes, ],
    pure_tumor_ids = annot$AOI.name == "Tumor",
    X = safeTME[sharedgenes, ],
    K = 5
  )
  expect_true(all(abs(mergedX.test - mergedX) < 1e-3))
})



# spatialdecon:
set.seed(0)
gres <- spatialdecon(
  norm = snr[sharedgenes, ],
  X = safeTME[sharedgenes, ],
  bg = replace(snr[sharedgenes, ], TRUE, 1),
  resid_thresh = 3, lower_thresh = 0.5
)

test_that("spatialdecon is as expected: beta", {
  expect_true(cor(as.vector(ires.test$beta), as.vector(gres$beta)) > 0.999)
})

test_that("spatialdecon is as expected: sigma", {
  expect_true(cor(as.vector(ires.test$sigmas), as.vector(gres$sigmas)) > 0.999)
})

test_that("spatialdecon is as expected: yhat", {
  expect_true(cor(as.vector(ires.test$yhat), as.vector(gres$yhat)) > 0.999)
})

test_that("spatialdecon is as expected: resids", {
  expect_true(all(abs(replace(ires.test$resids, is.na(ires.test$resids), 0) -
    replace(gres$resids, is.na(gres$resids), 0)) < 1e-2))
})

test_that("spatialdecon is as expected: p", {
  expect_true(all(abs(ires.test$p - gres$p) < 1e-2))
})



#### test SpatialDecon with TILs options ---------------------------------------
res <- spatialdecon(
  norm = snr,
  raw = raw,
  bg = replace(snr, TRUE, 1),
  cellmerges = safeTME.matches,
  is_pure_tumor = annot$AOI.name == "Tumor",
  cell_counts = annot$nuclei,
  n_tumor_clusters = 5
)

test_that("cell matching is as expected", {
  expect_equal(rownames(res.test$beta), rownames(res$beta))
})

test_that("spatialdeconTILs returned results as expected", {
  cols <- c("beta", "yhat", "sigma", "p", "t", "se", "prop_of_all", "prop_of_nontumor", 
            "cell.counts", "beta.granular", "sigma.granular", "cell.counts.granular", 
            "resids", "X")
  
  expect_true(all(cols %in% names(res)))
  expect_true(all(names(res) %in% cols))
})


test_that("spatialdeconTILs is as expected: beta", {
  expect_true(all(abs(res.test$beta - res$beta) < 1e-2))
})

test_that("spatialdeconTILs is as expected: sigma", {
  expect_true(all(abs(res.test$sigma - res$sigma) < 1e-2))
})

test_that("spatialdeconTILs is as expected: yhat", {
  expect_true(mean(abs(res.test$yhat - res$yhat)) < 1e-3)
})

test_that("spatialdeconTILs is as expected: resids", {
  expect_true(all(abs(replace(res.test$resids, is.na(res.test$resids), 0) -
    replace(res$resids, is.na(res$resids), 0)) < 1e-2))
})

test_that("spatialdeconTILs is as expected: p", {
  expect_true(all(abs(res.test$p - res$p) < 1e-2))
})

test_that("spatialdeconTILs is as expected: props", {
  expect_true(all(abs(res.test$prob_of_all - res$prob_of_all) < 1e-2))
  expect_true(all(abs(res.test$prob_of_nontumor - res$prob_of_nontumor) < 1e-2))
})

test_that("spatialdeconTILs is as expected: props of all", {
  expect_true(all(abs(res.test$prob_of_all - res$prob_of_all) < 1e-2))
  expect_true(all(abs(res.test$prob_of_nontumor - res$prob_of_nontumor) < 1e-2))
})

test_that("spatialdeconTILs is as expected: t", {
  expect_true(all(abs(res.test$t - res$t) < 1e-2))
})

test_that("spatialdeconTILs is as expected: se", {
  expect_true(all(abs(res.test$se - res$se) < 1e-2))
})

test_that("spatialdeconTILs is as expected: beta.granular", {
  expect_true(all(abs(res.test$beta.granular[,-7] - res$beta.granular[,-7]) < 2e-2))
})

### test reverse decon:
rdres <- suppressWarnings(reverseDecon(
  norm = snr,
  beta = res.test$beta,
  epsilon = 1
))

test_that("reverseDecon is as expected", {
    expect_true(all(abs(rdres.test$resids - rdres$resids) < 1e-1))
    expect_true(all(abs(rdres.test$yhat - rdres$yhat) < 1))
    expect_true(all(abs(rdres.test$coefs - rdres$coefs) < 100))
    expect_true(all(abs(rdres.test$cors - rdres$cors) < 1e-2, na.rm = TRUE))
    expect_true(all(abs(rdres.test$resid.sd - rdres$resid.sd) < 1e-2))
})




### test plotting functions
test_that("florets does not error", {
    expect_error(
        florets(
            x = annot$x, y = annot$y,
            b = res$beta.granular[!grepl("tumor", rownames(res$beta.granular)), ],
            legendwindow = TRUE
        ),
        NA
    )
})

test_that("TIL_barplot does not error", {
    expect_error(
        TIL_barplot(mat = res$beta, draw_legend = TRUE),
        NA
    )
})


### test matrix download:
test_that("matrix download works", {
  downloaded.X <- download_profile_matrix(species = "Human", 
                                          age_group = "Adult",
                                          matrixname = "Liver_HCA")
  expect_true(is.matrix(downloaded.X))
  expect_true(is.list(cellGroups))
  expect_true(all(colnames(downloaded.X) %in% unlist(cellGroups)))
  expect_true(all(downloaded.X == profile_matrix))
  
})
### test matrix creation
test_that("matrix creation works", {
  cellNames <- paste0("Cell", 1:1500)
  geneNames <- paste0("Gene", 1:1500)
  mtx <- matrix(data=sample(size = length(cellNames)*length(geneNames),
                            replace = TRUE,
                            x = 0:100,
                            prob = c(0.6784, rep(0.0075, 15), rep(0.005, 25),
                                   rep(0.002, 25),  rep(0.001, 35))),
                ncol = length(cellNames), nrow = length(geneNames),
                dimnames = list(geneNames,
                                cellNames))
  
  celltypes <- c("A", "B", "C", "D")
  cellAnnots <- as.data.frame(cbind(CellID=cellNames,
                                    cellType=sample(size = length(cellNames),
                                                    replace = TRUE,
                                                    x = celltypes,
                                                    prob = c(0.1, 0.4, 0.3, 0.2))))
  
  profile_matrix <- create_profile_matrix(mtx = mtx,
                                          cellAnnots = cellAnnots,
                                          cellTypeCol = "cellType",
                                          cellNameCol = "CellID",
                                          minGenes = 10,
                                          scalingFactor = 1,
                                          outDir = NULL
                                          )
  
  

  expect_true(is.matrix(profile_matrix))
  
  for(CT in celltypes){
    w2kp <- which(cellAnnots$cellType == CT)
    expect_true(all(rowMeans(mtx[,w2kp]) == profile_matrix[,CT]))
  }
  # expect_true(file.exists("Custom_profileMatrix.csv"))
  # 
  # unlink("Custom_profileMatrix.csv", force = TRUE)
})

### test collapseCellTypes:

# uncollapsed result:
res2 <- spatialdecon(
  norm = snr,
  raw = raw,
  bg = replace(snr, TRUE, 1),
  cellmerges = NULL,
  is_pure_tumor = annot$AOI.name == "Tumor",
  cell_counts = annot$nuclei,
  n_tumor_clusters = 5
)
# collapse them:
res2.collapsed <- collapseCellTypes(
  fit = res2,
  matching = safeTME.matches
)
# compare collapsed results from within spatialdecon vs. post-hoc:
test_that("collapseCellTypes works", {
  expect_true(all(abs(res2.collapsed$beta - res$beta) < 1e-2))
  expect_true(all(abs(res2.collapsed$sigmas - res$sigmas) < 1e-2))
  expect_true(all(abs(res2.collapsed$p - res$p) < 1e-2))
  expect_true(all(abs(res2.collapsed$t - res$t) < 1e-2))
  expect_true(all(abs(res2.collapsed$se - res$se) < 1e-2))
})


## test wrapper for seurat objects:
#make seurat object:
options(Seurat.object.assay.version = "v3")
seur <- SeuratObject::CreateSeuratObject(counts = raw, assay="Spatial")
test_that("runspatialdecon works on seurat objects - v3", {
  res <- runspatialdecon(seur)
  res2 <- spatialdecon(norm = raw, raw = raw, bg = 0.1)
  
  expect_true(is.matrix(res$beta)) # test beta is a matrix
  expect_true(is.matrix(res$yhat)) 
  expect_true(is.matrix(res$resids))
  expect_true(length(dim(res$sigmas)) == 3) # test sigmas is a 3d assar
  expect_true(is.matrix(res$p))
  expect_true(is.matrix(res$t))
  expect_true(is.matrix(res$se))
  expect_true(is.matrix(res$prop_of_all))
  expect_true(is.matrix(res$prop_of_nontumor))
  expect_true(is.matrix(res$X))
  
  expect_identical(res, res2)
})

options(Seurat.object.assay.version = "v5")
seur <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = raw, assay="Spatial"))
test_that("runspatialdecon works on seurat objects - v5", {
  res <- suppressWarnings(runspatialdecon(seur))
  
  res2 <- spatialdecon(norm = raw, raw = raw, bg = 0.1)
  
  expect_true(is.matrix(res$beta)) # test beta is a matrix
  expect_true(is.matrix(res$yhat)) 
  expect_true(is.matrix(res$resids))
  expect_true(length(dim(res$sigmas)) == 3) # test sigmas is a 3d assar
  expect_true(is.matrix(res$p))
  expect_true(is.matrix(res$t))
  expect_true(is.matrix(res$se))
  expect_true(is.matrix(res$prop_of_all))
  expect_true(is.matrix(res$prop_of_nontumor))
  expect_true(is.matrix(res$X))
  
  expect_identical(res, res2)
})

## test wrapper for GeoMxSet objects:
#make GeoMxSet object:

datadir <- system.file("extdata", "DSP_NGS_Example_Data", package = "GeomxTools")
demoData <- readRDS(file.path(datadir, "/demoData.rds"))

demoData <- GeomxTools::shiftCountsOne(demoData)
demoData <- GeomxTools::aggregateCounts(demoData)

demoData <- GeomxTools::normalize(demoData, "quant")

norm <- as.matrix(demoData@assayData$exprs_norm)

bg <- derive_GeoMx_background(norm = norm,
                              # access the probe pool information from the feature metadata
                              probepool = demoData@featureData@data$Module,
                              # access the names of the negative control probes
                              negnames = demoData@featureData@data$TargetName[demoData@featureData@data$Negative])

test_that("derive_GeoMx_background is in correct format", {
  expect_true(all(dim(bg) == dim(norm)))
  expect_true(all(colnames(bg) == colnames(norm)))
  expect_true(all(rownames(bg) == rownames(norm)))
})


res <- runspatialdecon(object = demoData, 
                       norm_elt = "exprs_norm",
                       raw_elt = "exprs")

test_that("runspatialdecon works on GeoMxSet objects", {
  res2 <- spatialdecon(norm = norm, bg = bg, 
                       raw = as.matrix(demoData@assayData$exprs))

  expect_true(is.matrix(res@phenoData@data$beta)) # test beta is a matrix
  expect_true(is.matrix(res@assayData$yhat)) 
  expect_true(is.matrix(res@assayData$resids))
  expect_true(length(dim(res@phenoData@data$sigmas)) == 3) # test sigmas is a 3d assar
  expect_true(is.matrix(res@phenoData@data$p))
  expect_true(is.matrix(res@phenoData@data$t))
  expect_true(is.matrix(res@phenoData@data$se))
  expect_true(is.matrix(res@phenoData@data$prop_of_all))
  expect_true(is.matrix(res@phenoData@data$prop_of_nontumor))
  expect_true(is.matrix(res@experimentData@other$SpatialDeconMatrix))
  
  expect_identical(res@phenoData@data$beta, t(res2$beta))
  expect_identical(res@assayData$yhat[match(rownames(res2$yhat), rownames(res@assayData$yhat)),], 
                   res2$yhat) 
  expect_identical(res@assayData$resids[match(rownames(res2$resids), rownames(res@assayData$resids)),], 
                   res2$resids)
  expect_identical(res@phenoData@data$p, t(res2$p))
  expect_identical(res@phenoData@data$t, t(res2$t))
  expect_identical(res@phenoData@data$se, t(res2$se))
  expect_identical(res@phenoData@data$prop_of_all, t(res2$prop_of_all))
  expect_identical(res@phenoData@data$prop_of_nontumor, t(res2$prop_of_nontumor))
  expect_identical(res@experimentData@other$SpatialDeconMatrix, res2$X)
})

sharedgenes <- intersect(rownames(safeTME), rownames(norm))

test_that("runmergeTumorIntoX works on GeoMxSet objects", {
    set.seed(0)
    mergedX <- mergeTumorIntoX(
      norm = norm,
      bg = bg,
      pure_tumor_ids = demoData@phenoData@data$`scan name` == "cw005 (PTL-10891) Slide1",
      X = safeTME[sharedgenes, ],
      K = 5
    )
    
    mergedX.test <- runMergeTumorIntoX(demoData, 
                                       norm_elt = "exprs_norm", 
                                       pure_tumor_ids = demoData@phenoData@data$`scan name` == "cw005 (PTL-10891) Slide1",
                                       X = safeTME[sharedgenes, ],
                                       K=5)
    
    expect_true(all(abs(mergedX.test - mergedX) < 1e-3))
})

test_that("runreverseDecon works on GeoMxSet objects", {
  rdres <- suppressWarnings(runReverseDecon(demoData, 
                                            norm_elt = "exprs_norm",
                                            beta = res@phenoData@data$beta,
                                            epsilon = 1))
  
  rdres.x <- suppressWarnings(reverseDecon(norm = norm,
                                           beta = t(res@phenoData@data$beta),
                                           epsilon = 1))
  
  
  expect_true(all(abs(rdres@featureData@data$coefs - rdres.x$coefs) < 1e-3))
  expect_true(all(abs(rdres@featureData@data$cors - rdres.x$cors) < 1e-3))
  expect_true(all(abs(rdres@featureData@data$resid.sd - rdres.x$resid.sd) < 1e-3))
  expect_true(all(abs(rdres@assayData$resids - rdres.x$resids) < 1e-3))
  expect_true(all(abs(rdres@assayData$yhat - rdres.x$yhat) < 1e-3))
})

test_that("runCollapseCellTypes works on GeoMxSet objects", {
  collapsed <- runCollapseCellTypes(object = res, 
                                    matching = safeTME.matches)
  
  expect_true(all(colnames(collapsed$beta) %in% names(safeTME.matches)))
  expect_lt(ncol(collapsed$beta), ncol(res$beta))
  expect_equal(nrow(collapsed$beta), nrow(res$beta))
})

