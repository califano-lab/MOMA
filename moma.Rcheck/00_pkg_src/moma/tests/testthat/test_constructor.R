


context("MOMA constructor")
library(moma)


test_that("build gbm example", {

	pathways <- list()
	pathways[['cindy']] = gbm.example$cindy
	pathways[['preppi']] = gbm.example$preppi
	momaObj <- moma.constructor(gbm.example$vipermat, gbm.example$rawsnp,
        	gbm.example$rawcnv, gbm.example$fusions, pathways,
        	gene.blacklist=gbm.example$mutSig,
        	gene.loc.mapping=gbm.example$gene.loc.mapping)

	expect_equal(ncol(momaObj$viper), 40)
	expect_equal(class(momaObj$viper), "matrix")
	expect_equal(ncol(momaObj$mut), 40)
})
