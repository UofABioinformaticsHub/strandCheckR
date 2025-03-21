test_that("GetWinFromBamFile Return Proper DataFrame", {
  file <- system.file("extdata", "ex1.bam", package="Rsamtools")
  win <- getStrandFromBamFile(file)
  require_cols <- c(
      "Seq","Start","End","NbPos","NbNeg","CovPos","CovNeg","MaxCoverage","File"
      )
  expect_true(all(require_cols %in% colnames(win)))
})
