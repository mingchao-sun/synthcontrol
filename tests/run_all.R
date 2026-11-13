# Run the whole test suite from the module folder:  Rscript tests/run_all.R
testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
