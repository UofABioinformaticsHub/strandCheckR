test_that(".summarizeHist works correctly", {
    files <- system.file(
        "extdata", c("s1.sorted.bam","s2.sorted.bam"), package="strandCheckR"
        )
    win <- getStrandFromBamFile(files)
    h <- .summarizeHist(win,groupBy = "File", normalizeBy = "File")
    expect_length(unique(h$File),2)
    expect_equal(sum(h$ReadCountProp),2)
    expect_true(
        all(unique(h$Coverage) %in% 
                c("(0,10]","(10,100]","(100,1000]","> 1000"))
        )
    expect_true(all(round(unique(h$PosProp)*100) %in% 0:100))
})

