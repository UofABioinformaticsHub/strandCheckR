#' @title Calculate the read fragments to be kept
#' @description Calculate the keeping probability of each read fragment based
#' on the keeping probability of the windows containing it. Then get the list
#' of read fragments to be kept.
#' @param fragments an IRanges object defining the starting, ending position of
#' each fragment
#' @param keptProbaW an Rle object define the kept probability of each
#' sliding window
#'
#' @return an integer vector of read fragment indices to be kept
#' @importFrom stats rbinom
#' @importFrom IRanges width start
#'
#' @keywords internal
.keptReadFragment <- function(fragments, keptProbaW) {
    if (length(fragments) > 0) {
        minL <- min(width(fragments))
        maxL <- max(width(fragments))
        proba <- rep(0, length(fragments))
        for (m in minL:maxL) {
            id <- which(width(fragments) == m)
            x <- keptProbaW[start(fragments)[id]]
            j <- 1
            while (j <= (m - 1)) {
                y <- keptProbaW[j + start(fragments)[id]]
                x <- (x >= y) * x + (y > x) * y # x <- pairwise max of x and y
                j <- j + 1
            }
            suppressWarnings(proba[id] <- rbinom(length(id), 1, as.vector(x)))
            rm(id)
            rm(x)
        }
        return(which(proba == 1))
    } else {
        return(c())
    }
}

