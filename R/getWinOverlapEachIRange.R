#' @title Get the ranges of sliding windows that overlap each range of an
#' IRanges object.
#'
#' @description Get the ranges of sliding windows that overlap each range of an
#' IRanges object.
#'
#' @details
#' This finds the windows that overlap each range of the input IRanges object.
#' Each range corresponds to a read fragment. This allows the total number of
#' read fragments within a window to be calculated simply using [IRanges::coverage()].
#'
#' @param x an IRanges object containing the start and end position of
#' each read fragment.
#' @param winWidth the width of the sliding window, 1000 by default.
#' @param winStep the step length to sliding the window, 100 by default.
#' @param readProp A read is considered to be included in a window if at least
#' \code{readProp} of it is in the window. Specified as a proportion.
#' @param maxWin The maximum window ID
#'
#' @return An IRanges object containing the index of the windows overlapping
#' each read fragment
#'
#' @examples
#' library(IRanges)
#' x <- IRanges(start=round(runif(100,1000,10000)),width=100)
#' getWinOverlapEachIRange(x)
#'
#' @importFrom IRanges start width IRanges
#' @export
getWinOverlapEachIRange <- function(
    x, winWidth = 1000L, winStep = 100L, readProp = 0.5, maxWin = Inf
    )
{
    if (length(readProp)>1 && length(readProp)!=length(x)){
        readProp <- readProp[1]
    }
    if (any(readProp <= 0) || any(readProp > 1)) {
        stop("readProp must be within the range (0, 1]")
    }

    stopifnot(is.numeric(maxWin))

    # Calculate the index of the first window that overlaps a fragment
    startWin <- ceiling(
        (start(x) - winWidth + (readProp) * width(x)) / winStep
    ) + 1
    startWin[startWin < 1] <- 1
    startWin[startWin > maxWin] <- maxWin

    # Calculate the index of the final window that overlaps a fragment
    endWin <- floor((end(x) - 1 - (readProp) * width(x)) / winStep) + 1
    endWin[endWin > maxWin] <- maxWin

    # Throw warning if there're some
    notSatisfy <- endWin < startWin
    if (sum(notSatisfy) > 0) {
        cat("Some reads do not have any window that overlap by the given read
proportion. Either read length is too big or window width is too small.\n")
    }
    # Return the complete set of windows for each fragment
    IRanges(start = startWin[!notSatisfy], end = endWin[!notSatisfy])
}
