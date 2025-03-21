test_that("getWinOverlapEachReadFragment can work", {
    file <- system.file("extdata", "ex1.bam", package="Rsamtools")
    sbp <- ScanBamParam(
        what=c("pos","cigar","strand","flag"),
        which=GRanges("seq1",IRanges(start=1,end=1600))
        )
    readInfo <- scanBam(BamFile(file), param = sbp)[[1]]
    firstReadIndex <- which(floor(readInfo$flag/64)%%2 == 1)
    wP <- getWinOverlapEachReadFragment(
        readInfo = readInfo, strand = "+",winWidth = 1000, winStep =100, 
        readProp = 0.5,useCoverage = TRUE,subset = firstReadIndex
        )
    posAlFirstRead <- firstReadIndex[readInfo$strand[firstReadIndex]=="+"]
    posAlFirstRead <- posAlFirstRead[!is.na(posAlFirstRead)]
    expect_true(
        is(wP$Win,"IRanges") && 
            setequal(posAlFirstRead,unique(mcols(wP$Win)$alignment))
        )
})
