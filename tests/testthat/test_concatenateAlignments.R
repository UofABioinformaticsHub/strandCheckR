test_that(".concatenateAlignments works correctly", {
    file <- BamFile(system.file("extdata","ex1.bam", package = "Rsamtools"))
    seqi <- seqinfo(file)
    sbp <- ScanBamParam(
        what = c("pos","cigar","strand"),
        which = GRanges(seqnames(seqi),IRanges(1,seqlengths(seqi))))  
    readInfo <- scanBam(file,param = sbp)
    seqInfo <- data.frame(
        Sequence = seqnames(seqi), Length = seqlengths(seqi),
        NbReads = vapply(readInfo, function(seq){length(seq$pos)}, integer(1)), 
        FirstBaseInPart = 0, LastBaseInPart = 0, FirstReadInPart = 0,
        LastReadInPart = 0
        )
    seqInfo <- .sequenceInfoInPartition(
        seqInfo = seqInfo, winWidth = 1000, winStep = 100)
    readInfoMerge <- .concatenateAlignments(
        readInfo = readInfo,seqInfo = seqInfo)
    expect_equal(
        c(readInfo[[1]]$cigar,readInfo[[2]]$cigar), readInfoMerge$cigar
        )
    expect_equivalent(
        readInfoMerge$pos[seqInfo$FirstReadInPart[2]:sum(seqInfo$NbReads)],
        readInfo[[2]]$pos +  seqInfo$LastBaseInPart[1]
        )
})
