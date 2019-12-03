## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
knitr::opts_chunk$set(warning = FALSE)
library(moma)

## ----explore data--------------------------------------------------------
names(gbm.example)
dim(gbm.example$vipermat)
gbm.example$vipermat[1:3, 1:3]

## ----pathways------------------------------------------------------------
pathways <- list()
pathways[['cindy']] = gbm.example$cindy
pathways[['preppi']] = gbm.example$preppi

## ----moma object---------------------------------------------------------
momaObj <- moma.constructor(gbm.example$vipermat, gbm.example$rawsnp,
        gbm.example$rawcnv, gbm.example$fusions, pathways,
        gene.blacklist=gbm.example$mutSig,
        gene.loc.mapping=gbm.example$gene.loc.mapping)

## ----interactions, results = "hide", message = FALSE---------------------
momaObj$runDIGGIT(fCNV=gbm.example$fCNV)

momaObj$makeInteractions(cindy.only=FALSE)

momaObj$Rank(use.cindy=TRUE)

## ----clustering----------------------------------------------------------
clustering.solutions <- momaObj$Cluster()

# pick the 3 cluster solution and save it to the moma Object
momaObj$sample.clustering <- clustering.solutions[[2]]$clustering

## ----saturation, message=FALSE, warning=FALSE, results="hide"------------
momaObj$saturationCalculation()

## ----checkpoints---------------------------------------------------------
cluster1.checkpoint <- momaObj$checkpoints[[1]]
print (cluster1.checkpoint[1:10])

## ----make plots, warning=FALSE-------------------------------------------

momaObj$makeSaturationPlots()


## ----visualize plots, fig.height=5, fig.width=7--------------------------

momaObj$identity.plots$bar.plots[[1]]

momaObj$identity.plots$curve.plots[[1]]


## ----plot customization, fig.height=5, fig.width=7-----------------------
library(ggplot2)
momaObj$identity.plots$bar.plots[[1]] +
  ggtitle("GBM Subtype 1") +
  theme(axis.text.y = element_text(size=5)) + 
  theme_dark()


