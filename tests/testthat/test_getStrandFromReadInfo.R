test_that("getStrandFromReadInfo return correct windows", {
    file <- system.file("extdata", "ex1.bam", package="Rsamtools")
    sbp <- ScanBamParam(
        what=c("pos","cigar","strand","flag"),
        which=GRanges("seq1",IRanges(start=1,end=1600))
        )
    readInfo <- scanBam(BamFile(file), param = sbp)[[1]]
    firstReadIndex <- which(floor(readInfo$flag/64)%%2 == 1)
    w1 <- getStrandFromReadInfo(readInfo,subset = firstReadIndex)
    expect_true(nrow(w1) == 7)
})
