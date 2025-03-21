#' @title Get the sliding windows that overlap a GRanges object
#'
#' @description Get the sliding windows that overlap a GRanges object.
#'
#' @details
#' This finds the windows that overlaps the positive/negative strand of a
#' GRanges object. The GRanges object, which is \code{mustKeepRanges} in the
#' \code{filterDNA} method, defines the coordinates of the ranges in the
#' reference genome that all reads mapped to those ranges must be kept by the
#' filtering method \code{filterDNA}.
#' This method makes use of the method \code{getWinOverlapEachIRange} by
#' pretending each given range as the range of a read. Since the widths of
#' \code{x} are not necessarily the same (as normal read lengths), we
#' use \code{nbOverlapBases} to specify the minimum number of bases that a
#' window should overlap with a range of \code{x}, instead of using proprotion
#' as \code{readProp} in \code{getWinOverlapEachIRange}.
#'
#' @param x a GRanges object, which defines the coordinates of
#' the ranges in the reference genome that all reads mapped to those ranges
#' must be kept by the filtering method \code{filterDNA}.
#' @param seqInfo a data frame that contains some key information of the
#' sequences
#' @param winWidth the length of the sliding window, 1000 by default.
#' @param winStep the step length to sliding the window, 100 by default.
#' @param nbOverlapBases a window is considered to overlap with a range of
#' \code{x} if it overlaps with at least \code{nbOverlapBases} bases.
#' @return A list of two logical vectors (for positive and negative strand)
#' defining which windows that overlap with the given GRanges object.
#' @examples
#' library(GenomicRanges)
#' x <- GRanges(seqnames = "10",ranges = IRanges(start = c(10000,15000),
#' end=c(20000,30000)),strand = c("+","-"))
#' seqInfo <- data.frame("Sequence"=10,"FirstBaseInPart"=1)
#' getWinOverlapGRanges(x,seqInfo)
#' seqInfo <- data.frame("Sequence"=10,"FirstBaseInPart"=10000000)
#' getWinOverlapGRanges(x,seqInfo)
#'
#' @importFrom GenomicRanges start<- end<- ranges<- ranges
#' @importFrom IRanges start end coverage
#' @importFrom BiocGenerics strand
#' @importFrom GenomeInfoDb seqlevels
#' @export
getWinOverlapGRanges <- function(
    x, seqInfo, winWidth = 1000L, winStep = 100L, nbOverlapBases = 1
)
{
    # Check the correct columns are in the seqInfo df
    reqCols <- c("Sequence","FirstBaseInPart")
    if (!all(reqCols %in% names(seqInfo))) {
        stop("seqInfo must contain the column ", reqCols)
    }
    stopifnot(is.numeric(winWidth) || is.numeric(winStep))

    # Calculate start/end position of each sequence in the partition
    for (i in seq_along(seqInfo$Sequence)) {
        r <- which(as.vector(seqnames(x)) == seqInfo$Sequence[i])
        if (length(r) > 0) {
            end(ranges(x)[r]) <-
                end(ranges(x)[r]) + seqInfo$FirstBaseInPart[i] - 1
            start(ranges(x)[r]) <-
                start(ranges(x)[r]) + seqInfo$FirstBaseInPart[i] - 1
        }
    }

    # Calculate the windows that overlap the '+' ranges of x
    mustKeepPos <- getWinOverlapEachIRange(
        ranges(x)[strand(x) != "-", ], winWidth, winStep,
        readProp = nbOverlapBases/width(x)
        )
    mustKeepPos <- coverage(mustKeepPos) > 0

    # Calculate the windows that overlap the '-' ranges of x
    mustKeepNeg <- getWinOverlapEachIRange(
        ranges(x)[strand(x) != "+", ], winWidth, winStep,
        readProp = nbOverlapBases/width(x)
        )
    mustKeepNeg <- coverage(mustKeepNeg) > 0

    list(Positive = mustKeepPos, Negative = mustKeepNeg)
}

