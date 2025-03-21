test_that(".keptProbaWin can work", {
    file <- BamFile(
        system.file("extdata/s1.sorted.bam", package = "strandCheckR")
        )
    sbp <- ScanBamParam(
        what = c("pos","strand","cigar"), 
        which = GRanges("10",ranges = IRanges(7900000,8000000))
        )
    readInfo <- scanBam(file = file, param = sbp)[[1]]
    winPosRecords <- getWinOverlapEachReadFragment(
        readInfo = readInfo, strand = "+", winWidth = 1000, winStep = 100, 
        readProp = 0.5, useCoverage = FALSE
        )
    winNegRecords <- getWinOverlapEachReadFragment(
        readInfo = readInfo, strand = "-", winWidth = 1000, winStep = 100, 
        readProp = 0.5, useCoverage = FALSE
        )
    probaWin <- .keptProbaWin(
        winPosAlignments = winPosRecords, winNegAlignments = winNegRecords, 
        winWidth =  winWidth, winStep = 100, threshold = 0.7, 
        pvalueThreshold = 0.05, errorRate = 0, mustKeepWin = NULL, minCov = 0, 
        maxCov = 0, getWin = FALSE, useCoverage = FALSE
        )
    if (sum(probaWin$Pos>0) >0 ){
        expect_gt(min(probaWin$Pos[probaWin$Pos>0]),0.5)    
    }
    if (sum(probaWin$Neg>0) > 0){
        expect_gt(min(probaWin$Neg[probaWin$Neg>0]),0.5)
    }
})
