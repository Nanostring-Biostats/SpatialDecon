context("spatialdecon")


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
  expect_true(all(abs(ires.test$beta - gres$beta) < 1e-2))
})

test_that("spatialdecon is as expected: sigma", {
  expect_true(all(abs(ires.test$sigmas - gres$sigmas) < 1e-2))
})

test_that("spatialdecon is as expected: yhat", {
  expect_true(all(abs(ires.test$yhat - gres$yhat) < 1e-1))
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

### test reverse decon:
rdres <- suppressWarnings(reverseDecon(
  norm = snr,
  beta = res.test$beta,
  epsilon = 1
))

test_that("reverseDecon is as expected: ", {
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
  downloaded.X <- download_profile_matrix("Mouse_Brain")
  expect_true(is.matrix(downloaded.X))
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
                                        scalingFactor = 1)
  
  

  expect_true(is.matrix(profile_matrix))
  
  for(CT in celltypes){
    w2kp <- which(cellAnnots$cellType == CT)
    expect_true(all(rowMeans(mtx[,w2kp]) == profile_matrix[,CT]))
  }
  expect_true(file.exists("Custom_profileMatrix.csv"))
  
  unlink("Custom_profileMatrix.csv", force = TRUE)
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


