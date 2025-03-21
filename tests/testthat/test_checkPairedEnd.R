test_that("checkPairedEnd works correctly", {
    fileP <- system.file("extdata", "ex1.bam", package="Rsamtools")
    p <- checkPairedEnd(fileP)
    fileS <- system.file("extdata","s1.sorted.bam", package = "strandCheckR")
    s <- checkPairedEnd(fileS)
    expect_true(p)
    expect_false(s)
})
