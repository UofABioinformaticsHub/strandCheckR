#' @title Calculate the strand information based on coverage
#'
#' @description Calculate the coverage coming from '+'/'-' reads 
#' in all sliding windows
#' 
#' @param winPosAlignments a list that has a `Coverage` field containing 
#' coverage coming from positive reads
#' @param winNegAlignments a list that has a `Coverage` field containing 
#' coverage coming from negative reads
#' @param winWidth the length of the sliding window, 1000 by default.
#' @param winStep the step length to sliding the window, 100 by default.
#' 
#' @return a list of two vectors, containing a positive/negative coverage of 
#' the input positive/negative windows
#' 
#' @importFrom IRanges Views
#' @importFrom S4Vectors Rle
#' 
#' @keywords internal
.calculateStrandCoverage <- function(
    winPosAlignments, winNegAlignments, winWidth = 1000L, winStep = 100L
    ) 
{
    # make sure winPosAlignments$Coverage and winNegAlignments$Coverage
    # have the same length to avoid some warnings afterward
    if (is.null(winPosAlignments)) lastBasePos <- 0
    else lastBasePos <- length(winPosAlignments$Coverage)
    if (is.null(winNegAlignments)) lastBaseNeg <- 0
    else lastBaseNeg <- length(winNegAlignments$Coverage)
    lastBase <- max(lastBasePos, lastBaseNeg)
    if (lastBase == 0) {
        covPos <- Rle(0,0)
        covNeg <- Rle(0,0)
        maxCoverage <- Rle(0,0)
    } else{
        if (lastBaseNeg < lastBase) {
            winNegAlignments$Coverage <- c(
                winNegAlignments$Coverage,rep(0, lastBase - lastBaseNeg)
            )
        }
        if (lastBasePos < lastBase) {
            winPosAlignments$Coverage <- c(
                winPosAlignments$Coverage, rep(0, lastBase - lastBasePos)
            )
        }
        
        nbWin <- max(1,ceiling((lastBase - winWidth)/winStep) + 1)
        st <- seq(1, (nbWin - 1) * winStep + 1, by = winStep)
        end <- seq(winWidth, (nbWin - 1) * winStep + winWidth, by = winStep)
        covPos <- Views(winPosAlignments$Coverage, start = st, end = end)
        covPos <- Rle(sum(covPos))
        covNeg <- Views(winNegAlignments$Coverage, start = st, end = end)
        covNeg <- Rle(sum(covNeg))
        
        # calculate the max coverage in each window
        maxCoverage <- Rle(max(Views(
            winPosAlignments$Coverage + winNegAlignments$Coverage, 
            start = st, end = end
            )))
    }
    return(list(CovPos = covPos, CovNeg = covNeg, MaxCoverage = maxCoverage))
}
