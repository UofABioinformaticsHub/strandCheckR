test_that("plotHist gives a ggplot object with correct data to plot", {
    file <- system.file("extdata", "ex1.bam", package="Rsamtools")
    win <- getStrandFromBamFile(file)
    g <- plotHist(win, groupBy = "Type")
    expect_true(is_ggplot(g))
    expect_equal(names(g@data),c("Type","PosProp","Coverage","ReadCountProp"))
})
