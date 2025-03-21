test_that(".sequenceInfoInPartition works correctly", {
    n <- round(runif(1,2,10))
    seqInfo <- data.frame(
        Sequence = paste0("seq",1:n), Length = round(runif(n,1e+6,1e+7)),
        NbReads = round(runif(n,1e+3,1e+4)), FirstBaseInPart = 0, 
        LastBaseInPart = 0, FirstReadInPart = 0, LastReadInPart = 0
        )
    winStep <- round(runif(1,50,150))
    winWidth <- winStep * 10
    seqInfo <- .sequenceInfoInPartition(seqInfo,winWidth,winStep)
    expect_equal(seqInfo$LastBaseInPart %% winStep,rep(0,n))
    expect_equal(seqInfo$LastBaseInPart[-n]+1,seqInfo$FirstBaseInPart[-1])    
    expect_equal(seqInfo$LastReadInPart[-n]+1,seqInfo$FirstReadInPart[-1])
    expect_equal(
        seqInfo$LastReadInPart-seqInfo$FirstReadInPart+1,seqInfo$NbReads
        )
})

