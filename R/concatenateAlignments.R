#' @title Concatenate a list of Alignments into One
#' 
#' @description Concatenate a list of Alignments from multiple sequences into a 
#' single object
#' 
#' @details This method take a list of alignments across one or more sequences 
#' as output by \code{scanBam} 
#' and concatenates them into a single set of alignments which may include 
#' multiple sequences
#' 
#' @param readInfo a list returned by \code{scanBam} function, each element 
#' correspond to a sequence, containing the information of strand, starting 
#' position, cigar string, and eventually flag, qname
#' @param seqInfo a data frame that contains some key information of the 
#' sequences
#' 
#' @return the concatenated alignments of the input list
#'
#' @keywords internal
.concatenateAlignments <- function(readInfo, seqInfo)
{   
    nameScanWhat <- names(readInfo[[1]])
    nField <- length(nameScanWhat)
    nSeq <- length(readInfo)
    chk <- vapply(readInfo, function(x){
        all(names(x) %in% nameScanWhat)}, logical(1))
    stopifnot(all(chk))
    stopifnot(nrow(seqInfo) == nSeq)
    
    readInfo <- do.call(c, readInfo)
    concatAlignments <- lapply(seq_along(nameScanWhat),function(n){
        unlist(readInfo[seq(n,nSeq*nField,nField)], use.names = FALSE)})
    names(concatAlignments) <- nameScanWhat
    
    # shift the position of each record after concatening the sequences
    iPos <- which(nameScanWhat == "pos")
    for (i in seq_along(seqInfo$Sequence)) {
        range <- seqInfo$FirstReadInPart[i]:seqInfo$LastReadInPart[i]
        concatAlignments$pos[range] <- 
            readInfo[[iPos + (i - 1) * nField]] + seqInfo$FirstBaseInPart[i] - 1
    }
    return(concatAlignments)
}

