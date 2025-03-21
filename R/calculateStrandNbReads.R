#' @title Calculate the strand information based on the number of reads
#'
#' @description Calculate the number of reads coming from '+'/'-' strands in 
#' all sliding wndows
#' @param winPosAlignments a list that has a `Win` field that contains 
#' information of sliding windows overalapping positive reads 
#' @param winNegAlignments a a list that has a `Win` field that contains 
#' information of sliding windows overalapping negative reads
#'  
#' @return a list of two vectors, containing a positive/negative number of 
#' reads of the input positive/negative windows
#' @importFrom IRanges coverage end
#' 
#' @keywords internal
.calculateStrandNbReads <- function(winPosAlignments, winNegAlignments) 
{   
    # Calculate strand information based on number of reads have the same 
    # length to avoid some warnings afterward
    if (is.null(winPosAlignments)) {
        NbPos <- NULL
        lastWinPos <- 0
    }
    else {
        NbPos <- coverage(winPosAlignments$Win)
        lastWinPos <- max(end(winPosAlignments$Win))
    }
    if (is.null(winNegAlignments)) {
        NbNeg <- NULL
        lastWinNeg <- 0
    }
    else {
        NbNeg <- coverage(winNegAlignments$Win)
        lastWinNeg <- max(end(winNegAlignments$Win))
    }
    
    # Find the last window in both sets of windows
    lastWin <- max(lastWinPos,lastWinNeg)
    if (lastWin == 0) {
        NbPos <- Rle(0,0)
        NbNeg <- Rle(0,0)
    } else{
        # Fill with zeroes if required make sure NbPos and NbNeg have 
        # the same length
        if (lastWinNeg < lastWin) 
            NbNeg <- c(NbNeg, rep(0, lastWin - lastWinNeg))
        if (lastWinPos < lastWin) 
            NbPos <- c(NbPos, rep(0, lastWin - lastWinPos))
    }
    return(list(NbPos = NbPos, NbNeg = NbNeg))
}
