#' @title Calculate the first/last base/read of each sequence within each part
#' of the partition.
#' @param seqInfo a data frame that contains some key information of the 
#' sequences
#' @param winWidth the length of sliding window 
#' @param winStep the step length to sliding the window 
#' @keywords internal

.sequenceInfoInPartition <- function(seqInfo, winWidth, winStep) 
{
    reqNames <- c("Sequence", "Length", "NbReads")
    stopifnot(all(reqNames %in% names(seqInfo)))
    lengthInPart <- pmax(seqInfo$Length,winWidth)
    lengthOfPart <- 
        winStep * ceiling(c(0, cumsum(as.numeric(lengthInPart)))/winStep)
    nbReadsOfPart <- c(0, cumsum(seqInfo$NbReads))
    # update seqInfo
    seqInfo$FirstBaseInPart <- lengthOfPart[-length(lengthOfPart)] + 1
    seqInfo$LastBaseInPart <- lengthOfPart[-1]
    seqInfo$FirstReadInPart <- nbReadsOfPart[-length(nbReadsOfPart)] + 1
    seqInfo$LastReadInPart <- nbReadsOfPart[-1]
    return(seqInfo)
}

