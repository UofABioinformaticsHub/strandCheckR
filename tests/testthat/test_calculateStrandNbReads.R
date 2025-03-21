test_that(".calculateStrandNbReads can work", {
    winPos <- list("Win"=IRanges(start = 1, end = ceiling(runif(100,1,10))))
    winNeg <- list("Win"=IRanges(start = 1, end = ceiling(runif(100,1,11))))
    w <- .calculateStrandNbReads(winPos,winNeg)
    m <- max(end(winPos$Win),end(winNeg$Win))
    expect_true(
        setequal(names(w),c("NbPos","NbNeg")) && (length(w$NbPos)==m) && 
            (length(w$NbNeg)==m)
        )
})
