# MOMA


Installation: 

	R CMD install moma_0.2.tar.gz

Test examples:

	The first example can be run from this directory, and includes GBM data:

		test/moma.gbm.R

Development Guide:

Roxygen2 is used for documetation. Use 'roxygen2::roxygenize()' to rebuild docs while in the 'moma' directory

Build check:
	
	R CMD build moma
	R CMD CHECK 	moma_1.0.0.tar.gz
	R CMD BiocCheck moma_1.0.0.tar.gz


Running unit tests:

	Install the package

	library(moma)
	library(devtools)
	devtools::test("moma")

