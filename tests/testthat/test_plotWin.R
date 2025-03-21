test_that("plotWin gives a ggplot object", {
    file <- system.file("extdata", "ex1.bam", package="Rsamtools")
    win <- getStrandFromBamFile(file)
    g <- plotWin(win, groupBy = "Type")
    expect_is(g,"ggplot")
})

