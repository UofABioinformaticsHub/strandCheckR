test_that("getWinOverlapEachIRange can work", {
    ranges <- IRanges(start = as.integer(runif(10,1,2000)),width = 100)
    w <- getWinOverlapEachIRange(ranges,readProp=0.5)
    s <- (start(w)-1)*100+1000
    s1 <- (start(w)-2)*100+1000
    e <- (end(w)-1)*100 + 1
    e1 <- (end(w))*100 + 1
    expect_true(
        all(s >= start(ranges)+50) && 
            all(start(w) == 1 | s1 < start(ranges)+50) && 
            all(e <= end(ranges)-50) && all(e1 > end(ranges)-50)
        )
})
