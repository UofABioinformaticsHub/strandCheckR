test_that(".calculateStrandCoverage can work", {
    winPos <- list("Coverage"=Rle(as.integer(runif(2000,1,50))))
    winNeg <- list("Coverage"=Rle(as.integer(runif(2000,1,50))))
    w <- .calculateStrandCoverage(winPos, winNeg, winWidth = 1000, winStep = 100)
    expect_true(
        setequal(names(w),c("CovPos","CovNeg","MaxCoverage")) && 
            length(w$CovPos)==11 && length(w$CovNeg)==11 && 
            length(w$MaxCoverage)==11
        )
})
