test_that(".getWinInSequence can work", {
    seqInfo <- data.frame(
        Sequence = c("seq1","seq2"), Length = c(1500,1500), NbReads = NA,  
        FirstBaseInPart = c(1,1501), LastBaseInPart = c(1500,3000), 
        FirstReadInPart = NA, LastReadInPart = NA
        )
    nbWin <- ceiling((3000 - 1000 + 1)/100)
    Win <- DataFrame(
        Type = Rle("",nbWin), Seq = Rle("",nbWin), Start = 1:nbWin, 
        End = Rle(0,nbWin), NbPos = ceiling(runif(nbWin,1,100)), 
        NbNeg = ceiling(runif(nbWin,1,100)), 
        CovPos = ceiling(runif(nbWin,100,1000)), 
        CovNeg = ceiling(runif(nbWin,100,1000)), 
        MaxCoverage = ceiling(runif(nbWin,10,20)), File = Rle("",nbWin)
        ) 
    w <- .getWinInSequence(Win,seqInfo)
    expect_equal(as.character(w$Seq),c(rep("seq1",6),rep("seq2",6)))
})
