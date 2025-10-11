#' @title Intersect the windows data frame with an annotation data frame
#'
#' @description Intersect the windows with an annotation data frame to get
#' features that overlap with each window
#'
#' @param windows data frame containing the strand information of the sliding
#' windows. Windows can be obtained using the function \code{getStrandFromBamFile}.
#' @param annotation a Grange object that you want to intersect with your
#' windows. It can have mcols which contains the information or features that
#' could be able to integrate to the input windows
#' @param getFeatureInfo whether to get the information of features in the
#' mcols of annotation data or not.
#' If FALSE the return windows will have an additional column indicating
#' whether a window overlaps with any range of the annotion data.
#' If TRUE the return windows will contain the information of features that
#' overlap each window
#' @param overlapCol the columnn name of the return windows indicating whether
#' a window overlaps with any range of the annotion data.
#' @param mcolsAnnot the column names of the mcols of the annotation data that
#' you want to get information
#' @param collapse character which is used collapse multiple features that
#' overlap with a same window into a string. If missing then we don't
#' collapse them.
#' @param ... used to pass parameters to GenomicRanges::findOverlaps
#'
#' @return the input windows DataFrame with some additional columns
#'
#' @seealso \code{\link{getStrandFromBamFile}}, \code{\link{plotHist}},
#' \code{\link{plotWin}}
#'
#' @examples
#' bamfilein = system.file('extdata','s2.sorted.bam',package = 'strandCheckR')
#' windows <- getStrandFromBamFile(file = bamfilein)
#' #add chr before chromosome names to be consistent with the annotation
#' windows$Seq <- paste0('chr',windows$Seq)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' annot <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' # get the transcript names that overlap with each window
#' windows <- intersectWithFeature(windows,annot,mcolsAnnot='tx_name')
#' # just want to know whether there's any transcript that
#' # overlaps with each window
#' windows <- intersectWithFeature(windows,annot,overlapCol='OverlapTranscript')
#' plotHist(windows,facets = 'OverlapTranscript')
#' plotWin(windows,facets = 'OverlapTranscript')
#'
#' @importFrom IRanges IRanges
#' @importFrom Seqinfo seqlevels
#' @importFrom GenomicRanges GRanges mcols findOverlaps
#' @importFrom S4Vectors from to
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom methods is
#' @export
intersectWithFeature <- function(
    windows, annotation, getFeatureInfo = FALSE, overlapCol = "OverlapFeature",
    mcolsAnnot, collapse, ...
    )
{
    # check annotation is a GRanges object
    stopifnot(is(annotation, "GRanges"))
    # check if windows contains required column
    reqWinCols <- c("Seq", "Start", "End")
    stopifnot(all(reqWinCols %in% colnames(windows)))
    w <- GRanges(
        seqnames = windows$Seq,
        ranges = IRanges(start = windows$Start, end = windows$End)
        )
    if (length(intersect(seqlevels(w), seqlevels(annotation))) == 0) {
        message("Windows do not have any sequence in annotation data.")
        return(windows)
    }

    # check mcolsAnnot parameter
    if (!missing(mcolsAnnot)) {
        if (length(mcolsAnnot) > 0) {
            mcolsAnnot <- intersect(mcolsAnnot, colnames(mcols(annotation)))
            if (length(mcolsAnnot) == 0) {
                message(
                    "mcols of annotation does not contain any column ",
                    paste0(mcolsAnnot, collapse = ",")
                    )
                getFeatureInfo <- FALSE
            } else {
                getFeatureInfo <- TRUE
            }
        } else {
            message("No column specified to get from mcols of annotation.")
            getFeatureInfo <- FALSE
        }
    }

    # Get overlap between windows and annotation
    ol <- findOverlaps(w, annotation, select = "all", ...)

    if (!getFeatureInfo) {
        windows[[overlapCol]] <- "NotOverlap"
        windows[[overlapCol]][as.integer(unique(from(ol)))] <- "Overlap"
    } else {
        ## check collapse
        stopifnot(missing(collapse) || is.character(collapse))

        ## Check mcolsAnnot
        if (missing(mcolsAnnot)) {
            mcolsAnnot <- mcols(annotation)
        }

        featureTo <- mcols(annotation)[to(ol), mcolsAnnot[1]]
        splitTo <- split(featureTo, f = from(ol), drop = TRUE)
        if (!missing(collapse)) {
            splitTo <- vapply(splitTo, function(s) {
                paste0(unique(s), collapse = collapse)
            }, character(1))
        }
        indexFrom <- as.integer(names(splitTo))
        windows[[mcolsAnnot[1]]] <- "unknown"
        windows[[mcolsAnnot[1]]][indexFrom] <- splitTo

        for (i in seq_along(mcolsAnnot)[-1]) {
            featureTo <- mcols(annotation)[to(ol), mcolsAnnot[i]]
            splitTo <- split(featureTo, f = from(ol))
            if (!missing(collapse)) {
                splitTo <- vapply(
                    splitTo,
                    function(s) {paste0(unique(s), collapse = collapse)},
                    character(1)
                    )
            }
            windows[[mcolsAnnot[i]]] <- "unknown"
            windows[[mcolsAnnot[i]]][indexFrom] <- splitTo
        }
    }
    return(windows)
}
