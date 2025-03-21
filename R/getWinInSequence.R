#' @title Get window data frame with the correct sequence name and position 
#' @description Get the correct sequence name and position for each window 
#' @param Win a data frame contains the strand information of every window 
#' @param seqInfo a data frame that contains some key information of the 
#' sequences
#' @param winWidth the length of sliding window 
#' @param winStep the step length to sliding the window 
#' @keywords internal
.getWinInSequence <- function(Win, seqInfo, winWidth = 1000L, winStep = 100L) 
{
    # Check the correct columns are in the seqInfo df
    reqCols <- c("FirstBaseInPart", "LastBaseInPart")
    if (!all(reqCols %in% names(seqInfo))) {
        stop("seqInfo must contain the columns ", reqCols)
    }
    stopifnot(is.numeric(winWidth) || is.numeric(winStep))
    
    for (id in seq_along(seqInfo$Sequence)) {
        currentChr <- seqInfo$Sequence[id]
        # get id of the first window of the sequence
        idFirst <- ceiling(seqInfo$FirstBaseInPart[id]/winStep)
        # get id of the last window of the sequence
        idLast <- ceiling((seqInfo$LastBaseInPart[id] - winWidth + 1)/winStep)
        # get the windows of the sequence
        idRows <- which(Win$Start >= idFirst & Win$Start <= idLast)
        Win$Seq[idRows] <- Rle(currentChr)
        Win$Start[idRows] <- (Win$Start[idRows] - idFirst) + 1
    }
    return(Win[Win$Seq!="",])
}
