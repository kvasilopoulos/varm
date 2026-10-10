# Fast validation against reference packages; run from the package root:
#   Rscript validation/run-all.R
pkgload::load_all(".", quiet = TRUE)
testthat::test_dir("tests/testthat", package = "varm", load_package = "none", stop_on_failure = TRUE)
