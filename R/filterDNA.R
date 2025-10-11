#' @title Filter reads comming from double strand sequences from a bam File
#' @description Filter putative double strand DNA from a strand specific RNA-seq
#' using a window sliding across the genome.
#' @param file the input bam file to be filterd. Your bamfile should be sorted
#' and have an index file located at the same path.
#' @param destination the file path where the filtered output will be written
#' @param statFile the file to write the summary of the results
#' @param sequences the list of sequences to be filtered
#' @param mapqFilter every read that has mapping quality below \code{mapqFilter}
#' will be removed before any analysis.
#' If missing, the entire bam file will be read.
#' @param paired if TRUE then the input bamfile will be considered as paired-end
#' reads. If missing, 100 thousands first reads will be inspected to test if
#' the input bam file in paired-end or single-end.
#' @param yieldSize by default is 1e6, i.e. the bam file is read by block of
#' reads whose size is defined by this parameter. It is used to pass to same
#' parameter of the scanBam function.
#' @param winWidth the length of the sliding window, 1000 by default.
#' @param winStep the step length to sliding the window, 100 by default.
#' @param readProp a read is considered to be included in a window if at least
#' \code{readProp} of it is in the window. Specified as a proportion.
#' 0.5 by default.
#' @param threshold the strand proportion threshold to test whether to keep a
#' window or not. 0.7 by default
#' @param pvalueThreshold the threshold for the p-value in the test of keeping
#' windows. 0.05 by default
#' @param useCoverage if TRUE, then the strand information in each window
#' corresponds to the sum of coverage coming from positive/negative reads;
#' and not the number of positive/negative reads as default.
#' @param mustKeepRanges a GRanges object; all reads that map to those ranges
#' will be kept regardless the strand proportion of the windows containing them.
#' @param getWin if TRUE, the function will not only filter the bam file but
#' also return a data frame containing the information of all windows of the
#' original and filtered bam file.
#' @param minCov if \code{useCoverage=FALSE}, every window that has less than
#' \code{minCov} reads will be rejected regardless the strand proportion.
#' If \code{useCoverage=TRUE}, every window has max coverage least than
#' \code{minCov} will be rejected. 0 by default
#' @param maxCov if \code{useCoverage=FALSE}, every window that has more than
#' \code{maxCov} reads will be kept regardless the strand proportion.
#' If \code{useCoverage=TRUE}, every window with max coverage more than
#' \code{maxCov} will be kept.
#' If 0 then it doesn't have effect on selecting window. 0 by default.
#' @param errorRate the probability that an RNA read takes the false strand.
#' 0.01 by default.
#'
#' @details filterDNA reads a bam file containing strand specific RNA reads, and
#' filter reads coming from putative double strand DNA.
#' Using a window sliding across the genome, we calculate the positive/negative
#' proportion of reads in each window.
#' We then use logistic regression to estimate the strand proportion of reads in
#' each window, and calculate the p-value when comparing that to a given
#' threshold.
#' Let \eqn{\pi} be the strand proportion of reads in a window.
#'
#' Null hypothesis for positive window: \eqn{\pi \le threshold}.
#'
#' Null hypothesis for negative window: \eqn{\pi \ge 1-threshold}.
#'
#' Only windows with p-value <= \code{pvalueThreshold} are kept. For a kept
#' positive window, each positive read in this window is kept with the
#' probability (P-M)/P where P be the number of positive reads, and M be the
#' number of negative reads. That is because those M negative reads are
#' supposed to come from double-strand DNA, then there should be also M
#' postive reads among the P positive reads come from double-strand DNA. In
#' other words, there are only (P-M) positive reads come from RNA. Each
#' negative read is kept with the probability equalling the rate that an RNA
#' read of your sample has wrong strand, which  is \code{errorRate}.
#' Similar for kept negative windows.
#'
#' Since each alignment can be belonged to several windows, then the
#' probability of keeping an alignment is the maximum probability defined by
#' all windows that contain it.
#'
#' @return if \code{getWin} is TRUE: a DataFrame object which could also be
#' obtained by the function \code{getStrandFromBamFile}
#'
#' @seealso \code{\link{getStrandFromBamFile}}, \code{\link{plotHist}},
#' \code{\link{plotWin}}
#'
#' @examples
#' file <- system.file('extdata','s2.sorted.bam',package = 'strandCheckR')
#' out_bam <- tempfile(fileext = ".bam")
#' out_log <- tempfile(fileext = ".log")
#' filterDNA(file, sequences = '10', destination = out_bam, statFile = out_log)
#'
#' @importFrom GenomicRanges GRanges ranges
#' @importFrom Seqinfo seqinfo seqnames seqlengths
#' @importFrom IRanges IRanges
#' @importFrom methods is
#' @importFrom Rsamtools filterBam
#' @importFrom S4Vectors FilterRules
#'
#' @export
filterDNA <- function(
    file, destination, statFile = "out.stat", sequences, mapqFilter = 0, paired,
    yieldSize = 1e+06, winWidth = 1000L, winStep = 100L, readProp = 0.5,
    threshold = 0.7, pvalueThreshold = 0.05, useCoverage = FALSE,
    mustKeepRanges, getWin = FALSE, minCov = 0, maxCov = 0, errorRate = 0.01
    )
{
    startTime <- proc.time()

    # Check the input is a BamFile. Convert if necessary
    if (length(file) > 1) {
        message("Multiple files provided. Only the first will be filtered")
        file <- file[1]
    }
    if (!is(file, "BamFile")) {
        tryCatch(file <- BamFile(file, yieldSize = yieldSize))
    }

    # Check the destination path exists & has the suffix bam
    stopifnot(file.exists(dirname(destination)))
    destination <- ifelse(
        !grepl("\\.bam$", destination[1]), paste0(destination[1], ".bam"),
        destination[1]
        )

    # Check valid mapqFilter value
    if (mapqFilter < 0 || !is.numeric(mapqFilter)) {
        stop("Invalid value for mapqFilter. Must be positive & numeric.")
    }

    # More input checks
    stopifnot(is.logical(getWin))
    winWidth <- tryCatch(as.integer(winWidth))
    winStep <- tryCatch(as.integer(winStep))
    stopifnot(threshold > 0 & threshold < 1)
    stopifnot(pvalueThreshold > 0 & pvalueThreshold < 1)
    stopifnot(all(
        is.numeric(maxCov), is.numeric(minCov), is.numeric(errorRate),
        is.numeric(readProp)
        ))
    stopifnot(is.logical(useCoverage))

    # Check the paired status
    if (missing(paired))
        paired <- checkPairedEnd(file$path)

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


    # Define the file to write the summary of the results
    writeStats <- !is.null(statFile)
    if (writeStats) {
        stopifnot(file.exists(dirname(statFile)))
        message("Summary will be written to ", statFile)
        file.create(statFile)
    }

    # Get the sequence list of the ranges which must be kept
    if (!missing(mustKeepRanges)) {
        stopifnot(is(mustKeepRanges, "GRanges"))
        allSequencesMustKeep <- levels(seqnames(mustKeepRanges))
    }

    # Initialise the data frame containing the sequence information
    seqInfo <- data.frame(
        Sequence = sequences, Length = lengthSeq, NbReads = 0,
        FirstBaseInPart = 0, LastBaseInPart = 0, FirstReadInPart = 0,
        LastReadInPart = 0
        )

    # Define what to scan from bam file
    scanWhat <- c("pos", "cigar", "strand")
    if (paired) scanWhat <- c(scanWhat, "flag")

    # Create a list to store the results which record to keep for each sequence
    toKeepRecords <- list()

    # Initialise the data frames containing sliding window information to
    # be returned when getWin=TRUE
    if (getWin)
        allWin <- list()

    # Read the bam file in block defined by the yieldSize
    remainSequences <- sequences
    while (length(remainSequences) > 0) {
        # Get & Summarise the reads Return the reads from the bam file as a
        # list, with each element containing reads from a single seq
        scanGR <- GRanges(
            seqnames = remainSequences,
            ranges = IRanges(start = 1, end = lengthSeq[remainSequences])
            )
        sbp <- ScanBamParam(
            what = scanWhat,  which = scanGR, mapqFilter = mapqFilter
            )
        readInfo <- scanBam(file, param = sbp)
        # Get the number of records in each sequence
        nbReads <- vapply(readInfo, function(seq){length(seq$pos)}, integer(1))
        # Get the sequences that have been read
        idReadSeq <- which(nbReads > 0)
        if (length(idReadSeq) == 0) {
            remainSequences <- c()
        }
        else {
            ind <- seqInfo$Sequence %in% remainSequences
            seqInfo$NbReads[ind] <- nbReads
            readSeq <- remainSequences[idReadSeq]
            message("Read sequences ", paste(readSeq, collapse = " "))
            idP <- which(sequences %in% readSeq)
            partName <- readSeq[1]
            #Calculate the next remainSequences to read
            idNext <- idReadSeq[length(idReadSeq)] + 1
            if (idNext <= length(remainSequences)) {
                nextRange <- idNext:length(remainSequences)
                remainSequences <- remainSequences[nextRange]
            } else {
                remainSequences <- c()
            }

            # Calculate First/Last Base/Read in each part of the partition
            seqInfo[idP, ] <- .sequenceInfoInPartition(
                seqInfo[idP,], winWidth, winStep
                )

            # Concatenate lists of mutiple sequences into one list
            readInfo <- .concatenateAlignments(
                readInfo[idReadSeq], seqInfo[idP,]
                )
            # Calculate the windows that overlap mustKeepRanges
            mustKeepWin <- list()
            if (!missing(mustKeepRanges)) {
                if (length(intersect(allSequencesMustKeep, readSeq)) > 0) {
                    mustKeepWin <- getWinOverlapGRanges(
                        mustKeepRanges[seqnames(mustKeepRanges) %in% readSeq],
                        seqInfo[idP, ], winWidth, winStep
                        )
                }
            }
            # Get the index of R1 and R2 reads and process
            # each subset separately
            if (paired) {
                firstReadIndex <- which(floor(readInfo$flag/64) %% 2 == 1)
                secondReadIndex <- which(floor(readInfo$flag/64) %% 2 == 0)
                if (length(firstReadIndex) == 0) {
                    subset <- list(NULL)
                    type <- "R2"
                } else if (length(secondReadIndex) == 0) {
                    subset <- list(NULL)
                    type <- "R1"
                } else {
                    subset <- list(R1 = firstReadIndex, R2 = secondReadIndex)
                    type <- "R1"
                }
            } else {
                subset <- list(NULL)
                type <- "SE"
            }

            # Initialise the alignment indices to be kept
            keptRecords <- c()
            for (s in seq_along(subset)) {
                # Get the ids of sliding windows containing
                # each '+'/'-' read fragment
                winPosRecords <- getWinOverlapEachReadFragment(
                    readInfo, "+", winWidth, winStep, readProp = readProp,
                    useCoverage = (useCoverage || getWin), subset[[s]]
                    )
                winNegRecords <- getWinOverlapEachReadFragment(
                    readInfo, "-", winWidth, winStep, readProp = readProp,
                    useCoverage = (useCoverage || getWin), subset[[s]]
                    )

                # Calculate the keeping probability of each sliding window
                probaWin <- .keptProbaWin(
                    winPosRecords, winNegRecords, winWidth, winStep,
                    threshold, pvalueThreshold, errorRate, mustKeepWin,
                    minCov, maxCov, getWin = getWin, useCoverage = useCoverage
                    )

                # Calculate the '+'/'-' read fragments to be kept
                keptPosRecord <- .keptReadFragment(
                    winPosRecords$Win, probaWin$Pos
                    )
                keptNegRecord <- .keptReadFragment(
                    winNegRecords$Win, probaWin$Neg
                    )


                # Infer the index of kept alignments within the partition
                kept <- c(
                    unique(mcols(winPosRecords$Win)$alignment[keptPosRecord]),
                    unique(mcols(winNegRecords$Win)$alignment[keptNegRecord])
                    )

                # If getWin=TRUE, then return the strand information of sliding
                # windows from the orignial and filtered files
                if (getWin) {
                    # get the window information of filtered file
                    winA <- getStrandFromReadInfo(
                        readInfo, winWidth, winStep, readProp,subset = kept
                        )
                    # get the correct position of windows in each sequence
                    # of partition

                    win <- .getWinInSequence(
                        probaWin$Win, seqInfo[idP,], winWidth, winStep
                        )
                    winA <- .getWinInSequence(
                        winA, seqInfo[idP,], winWidth, winStep
                        )
                    # assign appropriate file name to each window data
                    win$File <- Rle(file$path,nrow(win))
                    winA$File <- Rle(destination,nrow(winA))
                    if (s == 1) {
                        win$Type <- Rle(type,nrow(win))
                        winA$Type <- Rle(type,nrow(winA))
                        allWin[[partName]] <- rbind(win, winA)
                    } else {
                        win$Type <- Rle("R2",nrow(win))
                        winA$Type <- Rle("R2",nrow(winA))
                        allWin[[partName]] <- rbind(
                            allWin[[partName]], rbind(win, winA)
                            )
                    }
                }
                # Add the current set to the vector of kept records
                keptRecords <- c(keptRecords, kept)
                # Tidy up the memory a little
                rm(winPosRecords, winNegRecords)
                rm(probaWin)
            }
            rm(readInfo)
            toKeepRecords[[partName]] <- rep(FALSE, sum(seqInfo$NbReads[idP]))
            toKeepRecords[[partName]][keptRecords] <- TRUE
            if (writeStats) {
                cat(
                    "Sequences: ", paste(readSeq, collapse = ","),
                    "\nNumber of reads: ",
                    sum(seqInfo$NbReads[idP]),"\nNumber of kept reads: ",
                    length(keptRecords),"\n", file = statFile, append = TRUE
                )
            }
            rm(keptRecords)
        }
    }
    # write the kept records into destination file
    scanGR <- GRanges(
        seqnames = sequences, IRanges(start = 1, end = seqInfo$Length)
        )
    filterBam(
        file = file, destination = destination,
        param = ScanBamParam(
            what = "rname", mapqFilter = mapqFilter, which = scanGR
            ),
        filter = FilterRules(
            list(
                Keep = function(x) {
                    seq <- as.character(x$rname[1])
                    return(toKeepRecords[[seq]])
                    }
                )
            )
        )
    nbKepReads <- vapply(toKeepRecords,sum,integer(1))
    if (writeStats) {
        cat("Summary:\n", file = statFile, append = TRUE)
        cat(
            "Number of original reads: ", sum(seqInfo$NbReads),
            "\nNumber of kept reads: ", sum(nbKepReads),
            "\nRemoval proportion: ",
            (sum(seqInfo$NbReads) - sum(nbKepReads))/sum(seqInfo$NbReads), "\n",
            file = statFile, append = TRUE
        )
        endTime <- proc.time()
        cat(
            "Total elapsed time: ", (endTime - startTime)[[3]]/60, " minutes\n",
            file = statFile, append = TRUE
        )
    }
    if (getWin) {
        allWin <- do.call(rbind, allWin)
        allWin$Start <- (allWin$Start - 1) * winStep + 1
        allWin$End <- allWin$Start + winWidth - 1
        return(allWin)
    }
}
