#' @title Test whether a bam file is single-end or paired-end
#' @description Check the first 100000 first reads of the bam file to see 
#' whether it is single-end or paired-end 
#' @param file the input bam file. Your bamfile should be sorted and have an 
#' index file located at the same path as well.
#' @param yieldSize the number of reads to be checked, 100000 by default.
#' 
#' @return return TRUE if the input file is paired end, and FALSE if it is 
#' single end
#' @export
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBam
#' @importFrom Rsamtools ScanBamParam
#' @examples
#' file <- system.file('extdata','s1.sorted.bam',package = 'strandCheckR')
#' checkPairedEnd(file)
#' @export

checkPairedEnd <- function(file, yieldSize = 1e+05) 
{
    message("Testing paired end by checking the first ", yieldSize, 
            " reads of file ", file)
    checkFile <- BamFile(file, yieldSize = yieldSize)
    flag <- scanBam(checkFile, param = ScanBamParam(what = "flag"))[[1]]$flag
    paired <- any(flag%%2 == 1)
    if (paired) {
        message("Your bam file is paired end") 
    } else {
        message("Your bam file is single end")
    }
    return(paired)
}

