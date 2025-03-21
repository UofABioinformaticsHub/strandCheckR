#' @title Get the window ranges that overlap each read fragment
#'
#' @description Calculate the window ranges that overlap each read fragment
#'
#' @param readInfo a list contains the read information
#' @param strand the considering strand
#' @param winWidth the width of the sliding window, 1000 by default.
#' @param winStep the step length to sliding the window, 100 by default.
#' @param readProp a read fragment is considered to be included in a window if
#' and only if at least \code{readProp} percent of it is in the window.
#' @param subset if we consider only a subset of the input reads
#' @param useCoverage either base on coverage or number of reads
#'
#' @return If \code{useCoverage=FALSE}: an IRanges object which contains the
#' range of sliding windows that overlap each read fragment.
#' If \code{useCoverage=TRUE}: a list of two objects, the first one is the
#' later IRanges object, the second one is an integer-Rle object which contains
#' the coverage of the input \code{readInfo}
#' @examples
#' library(Rsamtools)
#' file <- system.file('extdata','s2.sorted.bam',package = 'strandCheckR')
#' readInfo <- scanBam(file, param =
#' ScanBamParam(what = c("pos","cigar","strand")))
#' getWinOverlapEachReadFragment(readInfo[[1]],"+",1000,100,0.5)
#'
#' @importFrom GenomicAlignments extractAlignmentRangesOnReference
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom IRanges IRanges coverage
#' @export
getWinOverlapEachReadFragment <- function(
    readInfo, strand, winWidth, winStep, readProp, useCoverage = FALSE,
    subset = NULL
    )
{
    stopifnot(all(c("strand","pos","cigar") %in% names(readInfo)))
    if (is.null(subset)) {
        index <- which(readInfo$strand == strand)
    } else{
        index <- subset[which(readInfo$strand[subset] == strand)]
    }
    if (length(index) > 0) {
        # Read CIGAR strings to get fragment positions relative to windows
        position <- extractAlignmentRangesOnReference(
            readInfo$cigar[index], pos = readInfo$pos[index]
            )
        position <- as.data.frame(position)
        # Get the final window number
        maxWin <- ceiling((max(position$end) - winWidth) / winStep) + 1
        range <- IRanges(position$start, position$end, position$width)
        winrange <- getWinOverlapEachIRange(
            range, winWidth, winStep, readProp, maxWin
            )
        mcols(winrange) <- data.frame(alignment = index[position$group])
        if (useCoverage == TRUE) {
            return(list(Win = winrange, Coverage = coverage(range)))
        } else {
            return(list(Win = winrange))
        }
    } else{
        return(NULL)
    }
}
