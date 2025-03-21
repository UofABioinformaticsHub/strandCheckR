test_that(".keptReadFragment works properly", {
    file <- BamFile(
        system.file("extdata/s1.sorted.bam",package = "strandCheckR")
        )
    sbp <- ScanBamParam(
        what = c("pos","strand","cigar"),
        which = GRanges("10",ranges = IRanges(7800000,8000000))
        )
    readInfo <- scanBam(file, param = sbp)[[1]]
    winPosRecords <- getWinOverlapEachReadFragment(
        readInfo = readInfo, strand = "+", winWidth = 1000, 
        winStep = 100, readProp = 0.5, useCoverage = FALSE
        )
    winNegRecords <- getWinOverlapEachReadFragment(
        readInfo = readInfo, strand = "-", winWidth = 1000, 
        winStep = 100, readProp = 0.5, useCoverage = FALSE
        )
    probaWin <- .keptProbaWin(
        winPosAlignments = winPosRecords, winNegAlignments = winNegRecords, 
        winWidth = 1000, winStep = 100, threshold = 0.7, 
        pvalueThreshold = 0.05, errorRate = 0, mustKeepWin = NULL, 
        minCov = 0, maxCov = 0, getWin = FALSE, useCoverage = FALSE
        )
    keptPosRecord <- .keptReadFragment(
        fragments = winPosRecords$Win, keptProbaW = probaWin$Pos
        )
    keptNegRecord <- .keptReadFragment(
        fragments = winNegRecords$Win, keptProbaW = probaWin$Neg
        )
    pos <- unique(mcols(winPosRecords$Win)$alignment[keptPosRecord])
    neg <- unique(mcols(winNegRecords$Win)$alignment[keptNegRecord])
    expect_equal(as.character(unique(readInfo$strand[pos])),"+")
    expect_equal(as.character(unique(readInfo$strand[neg])),"-")
    
    pAlignments <- mcols(winPosRecords$Win)$alignment
    p <- which(pAlignments %in% pos)
    probaPos <- max(Views(probaWin$Pos, winPosRecords$Win[p]))
    minP <- min(vapply(
        pos, function(i){max(probaPos[which(pAlignments[p]==i)])}, numeric(1)
        )) 
    expect_gt(minP,0.7)
    
    nAlignments <- mcols(winNegRecords$Win)$alignment
    n <- which(nAlignments %in% neg)
    probaNeg <- max(Views(probaWin$Neg, winNegRecords$Win[n]))
    minN <- min(vapply(
        neg, function(i){max(probaNeg[which(nAlignments[n]==i)])}, numeric(1)
        )) 
    expect_gt(minN,0.7)
})
