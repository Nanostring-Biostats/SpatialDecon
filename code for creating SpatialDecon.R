##### this script contains the code used to build the package.
rm(list=ls())

### strategy for shiny-in-a-package deployment: taken from https://www.r-bloggers.com/packaging-shiny-applications-a-deep-dive/
### code written to emulate: https://github.com/MangoTheCat/shinyAppDemo


#### package building code: from https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/


library("devtools")
library(roxygen2)
library(qpdf)
library(BiocCheck)
library("styler")
#library("formatR")
library("lintr")

setwd("SpatialDecon")
#usethis::use_tidy_style(indent_by = 4)
#styler::style_pkg(transformers = styler::tidyverse_style(indent_by = 4))

devtools::document()
devtools::test() # run unit tests
devtools::check()
setwd("..")
BiocCheck("SpatialDecon")
BiocCheckGitClone("SpatialDecon")
setwd("SpatialDecon")
devtools::build()
devtools::build(binary = TRUE)
setwd("..")
install("SpatialDecon")
#setwd("SpatialDecon")
library("SpatialDecon")


#########################################
#### Creating data objects for SpatialDecon   ####
#########################################

# create colors:
cellcols = c()
cellcols["CD4.T.cells"] = "red"
cellcols["CD8.T.cells"] = "firebrick"
cellcols["Treg"] = "#FF66FF"

cellcols["T.CD4.naive"] = "#CC0000"
cellcols["T.CD4.memory"] = "#FF0000"
cellcols["T.CD8.naive"] = "#FF6633"
cellcols["T.CD8.memory"] = "#FF9900"

cellcols["NK"] = "grey10"
cellcols["B"] = "darkblue"
cellcols["B.naive"] = "#000099"
cellcols["B.memory"] = "#0000FF"
cellcols["plasma"] = "#3399CC"
cellcols["pDC"] = "#00FFFF"
cellcols["pDCs"] = "#00FFFF"
cellcols["macrophages"] = "#006600"
cellcols["monocytes"] = "#33CC00"
cellcols["monocytes.C"] = "#66CC66"
cellcols["monocytes.NC.I"] = "#33CC00"
cellcols["mDCs"] = "#00FF00"
cellcols["neutrophils"] = "#9966CC"
cellcols["mast"] = "#FFFF00"
cellcols["fibroblasts"] = "#999999"
cellcols["endothelial.cells"] = "#996633"
cellcols["tumor"] = "#333333"

# see how they look:
tempb = rep(1, length(cellcols))
names(tempb) = names(cellcols)
barplot(tempb, las = 2, col = cellcols)



save(cellcols, file = "SpatialDecon/data/cellcols.RData")
# create test dataset:
#source("create test results.R")


