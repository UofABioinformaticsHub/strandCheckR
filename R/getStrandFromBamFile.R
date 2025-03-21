#' @title Get the strand information of all windows from bam files
#'
#' @description Get the number of positive/negative reads/coverage of all
#' slding windows from the bam input files
#'
#' @param files the input bam files. Your bamfiles should be sorted and have
#' their index files located at the same path.
#' @param sequences character vector used to restrict analysed alignments to a
#' subset of chromosomes (i.e. sequences) within the provided bam file.
#' These correspond to chromosomes/scaffolds of the reference genome to which
#' the reads were mapped. If absent, the whole bam file will be read.
#' NB: This must match the chromosomes as defined in your reference genome.
#' If the reference chromosomes were specified using the 'chr' prefix, ensure
#' the supplied vector matches this specification.
#' @param mapqFilter every read that has mapping quality below
#' \code{mapqFilter} will be removed before any analysis.
#' @param yieldSize by default is 1e6, i.e. the bam file is read by block of
#' reads whose size is defined by this parameter. It is used to pass to same
#' parameter of the scanBam function.
#' @param winWidth the width of the sliding window, 1000 by default.
#' @param winStep the step length to sliding the window, 100 by default.
#' @param readProp A read is considered to be included in a window if at least
#' \code{readProp} of it is in the window. Specified as a proportion.
#' 0.5 by default.
#' @param paired if TRUE then the input bamfile will be considered as
#' paired-end reads. If missing, 100 thousands first reads will be inspected to
#' test if the input bam file in paired-end or single-end.
#'
#' @details
#' This function moves along the specified chromosomes (i.e. sequences) using a
#' sliding window approach, and counts the number of reads in each window which
#' align to the +/- strands of the reference genome.
#' As well as the number of reads, the total coverage for each strand is also
#' returned for each window, representing the total number of bases covered.
#'
#' Average coverage for the entire window can be simply calculated by dividing
#' the total coverage by the window size.
#'
#' @return a DataFrame object containing the number of positive/negative reads
#' and coverage of each window sliding across the bam file.
#' The returned DataFrame has 10 columns:
#'
#' Type: can be either SE if the input file contains single-end reads,
#' or R1/R2 if the input file contains paired-end reads.
#'
#' Seq: the reference sequence (chromosome/scaffold) that the reads were mapped
#' to.
#'
#' Start: the start position of the sliding window.
#'
#' End: the end position of the sliding window.
#'
#' NbPos/NbNeg: number of positive/negative reads that overlap the sliding
#' window.
#'
#' CovPos/CovNeg: number of bases coming from positive/negative reads that were
#' mapped in the sliding window.
#'
#' MaxCoverage: the maximum coverage within the sliding window.
#'
#' File: the name of the input file.
#'
#' @seealso \code{\link{filterDNA}}, \code{\link{plotHist}},
#' \code{\link{plotWin}}
#' @export
#' @importFrom IRanges Views
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom Rsamtools bamMapqFilter<- bamMapqFilter ScanBamParam BamFileList
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is
#' @examples
#' file <- system.file('extdata','s1.sorted.bam',package = 'strandCheckR')
#' win <- getStrandFromBamFile(file,sequences='10')
#' win

getStrandFromBamFile <- function(
    files, sequences, mapqFilter = 0, yieldSize = 1e+06, winWidth = 1000L,
    winStep = 100L, readProp = 0.5, paired
    )
{
    # Check the input is a BamFileList. Convert if necessary
    if (!is(files, "BamFileList")) {
        tryCatch(files <- BamFileList(files, yieldSize = yieldSize))
    }
    # Check valid mapqFilter value
    if (mapqFilter < 0 || !is.numeric(mapqFilter)) {
        stop("Invalid value for mapqFilter. Must be positive & numeric.")
    }

    # Check valide readProp
    if (!is.numeric(readProp) || readProp <= 0 || readProp > 1) {
        stop("Invalid value for readProp. Must be numeric and in (0,1].")
    }

    # Create a list to store the windows for each file
    allWin <- vector("list", length(files))

    # Read through each bam file
    for (b in seq_along(files)) {
        file <- files[[b]]
        # Check the paired status
        if (missing(paired)) {
            paired <- checkPairedEnd(file$path)
        }

        # Get the seqinfo object & all genomic information
        sq <- seqinfo(file)
        allRefSequences <- seqnames(sq)

        # Subset sequences if required
        if (!missing(sequences)) {
            sequenceId <- which(allRefSequences %in% sequences)
            stopifnot(length(sequenceId) > 0)
            sq <- sq[allRefSequences[sequenceId]]
            sequences <- seqnames(sq)
        } else {
            sequences <- allRefSequences
        }
        lengthSeq <- seqlengths(sq)

        # Initialise the data frames for window information to be returned
        allWin[[b]] <- list()

        # Initialise the data frames containing the sequence information in
        # each part
        seqInfo <- data.frame(
            Sequence = sequences, Length = lengthSeq, NbReads = 0,
            FirstBaseInPart = 0, LastBaseInPart = 0, FirstReadInPart = 0,
            LastReadInPart = 0
            )

        # what to scan from bam file
        scanWhat <- c("pos", "cigar", "strand")
        if (paired) scanWhat <- c(scanWhat, "flag")

        message("Reading file ", file$path)

        # Read the bam file in block defined by the yieldSize
        remainSequences <- sequences
        while (length(remainSequences) > 0) {
            # Get & Summarise the reads Return the reads from the bam file as
            # a list, with each element containing reads from a single seq
            scanGR <- GRanges(
                seqnames = remainSequences,
                ranges = IRanges(start = 1,end = lengthSeq[remainSequences])
                )
            sbp <- ScanBamParam(
                what = scanWhat,  which = scanGR, mapqFilter = mapqFilter
                )
            readInfo <- scanBam(file, param = sbp)
            # Get the number of records in each sequence and update seqInfo
            nbReads <- vapply(
                readInfo, function(seq){length(seq$pos)}, integer(1)
                )
            ind <- seqInfo$Sequence %in% remainSequences
            seqInfo$NbReads[ind] <- nbReads
            # Get the sequences that have been read
            idReadSeq <- which(nbReads > 0)
            if (length(idReadSeq)==0){
                remainSequences <- c()
            } else {
                readSeq <- remainSequences[idReadSeq]
                message("Read sequences ", paste(readSeq,collapse = " "))
                idP <- which(sequences %in% readSeq)
                partName <- readSeq[1]
                #Calculate the next remainSequences to read
                idNext <- idReadSeq[length(idReadSeq)] + 1
                if (idNext<=length(remainSequences)){
                    nextRange <- idNext:length(remainSequences)
                    remainSequences <- remainSequences[nextRange]
                } else {
                    remainSequences <- c()
                }
                # Calculate First/Last Base/Read in each part of the partition
                seqInfo[idP,] <- .sequenceInfoInPartition(
                    seqInfo[idP,], winWidth, winStep
                    )

                # Concatenate lists of mutiple sequences into one list
                readInfo <- .concatenateAlignments(
                    readInfo[idReadSeq], seqInfo[idP,]
                    )
                if (paired) {
                    firstReadIndex <- which(floor(readInfo$flag/64)%%2 == 1)
                    secondReadIndex <- which(floor(readInfo$flag/64)%%2 == 0)
                    if (length(firstReadIndex) == 0) {
                        subset <- list(NULL)
                        type <- "R2"
                    } else if (length(secondReadIndex) == 0) {
                        subset <- list(NULL)
                        type <- "R1"
                    } else {
                        subset <- list(R1 = firstReadIndex,
                            R2 = secondReadIndex)
                        type <- "R1"
                    }
                } else {
                    subset <- list(NULL)
                    type <- "SE"
                }
                for (s in seq_along(subset)) {
                    win <- getStrandFromReadInfo(
                        readInfo, winWidth, winStep, readProp, subset[[s]]
                        )
                    if (!is.null(win)){
                        win <- .getWinInSequence(
                            win, seqInfo[idP, ], winWidth, winStep
                            )
                        if (s == 1) {
                            win$Type <- Rle(type, nrow(win))
                            allWin[[b]][[partName]] <- win
                        } else {
                            win$Type <- Rle("R2", nrow(win))
                            allWin[[b]][[partName]] <-
                                rbind(allWin[[b]][[partName]], win)
                        }
                    } else{
                        if (s==1) allWin[[b]][[partName]] <- NULL
                    }
                }
            }
        }
        # combind the windows of all read parts
        if (length(allWin[[b]])==0){
            allWin[[b]] <- DataFrame(
                Type = Rle("",0), Seq = Rle("",0), Start = Rle(0,0),
                End = Rle(0,0), NbPos = Rle(0,0), NbNeg = Rle(0,0),
                CovPos = Rle(0,0), CovNeg = Rle(0,0), MaxCoverage = Rle(0,0)
                )
        } else{
            allWin[[b]] <- do.call(rbind, allWin[[b]])
        }
        allWin[[b]]$File <- Rle(file$path, nrow(allWin[[b]]))
    }
    allWin <- do.call(rbind, allWin)
    allWin$Start <- (allWin$Start - 1) * winStep + 1
    allWin$End <- allWin$Start + winWidth - 1
    return(allWin)
}
