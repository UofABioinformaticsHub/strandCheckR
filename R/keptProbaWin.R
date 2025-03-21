#' @title Get the probability of being kept for each window
#' @description Calculate the keeping probability of each window based on its
#' positive/negative proportion
#' @param winPosAlignments an object returned by getWinOverlapEachReadFragment for
#' positive reads
#' @param winNegAlignments an object returned by getWinOverlapEachReadFragment for
#' negative reads
#' @param winWidth the width of the sliding window, 1000 by default.
#' @param winStep the winStep length to sliding the window, 100 by default.
#' @param threshold the strand proportion threshold to test whether to keep a
#' window or not.
#' @param pvalueThreshold threshold of p-value
#' @param minCov In the case that \code{useCoverage=FALSE}, if a window has
#' less than \code{minCov} reads, then it will be rejected regardless of the
#' strand proportion.
#' For the case that \code{useCoverage=TRUE}, if a window has max coverage
#' least than \code{minCov}, then it will be rejected. 0 by default
#' @param maxCov In the case that \code{useCoverage=FALSE}, if a window has
#' more than \code{maxCov} reads, then it will be kept regardless of the strand
#' proportion.
#' For the case that \code{useCoverage=TRUE}, if a window has max coverage
#' more than \code{maxCov}, then it will be kept. If 0 then it doesn't have
#' effect on selecting window. 0 by default.
#' @param errorRate the probability that an RNA read takes the false strand.
#' 0.01 by default
#' @param useCoverage if TRUE, then the strand information in each window
#' corresponds to the sum of coverage coming from positive/negative reads;
#' and not the number of positive/negative reads as default.
#' @param mustKeepWin the windows that must be kept regardless their strand
#' proportion
#' @param getWin if TRUE, the function will return a data frame containing the
#' information of all windows. It's FALSE by default.
#'
#' @return A list of 2 numeric-Rle objects containing keeping probability of
#' each +/- alignments.
#' If \code{getWin=TRUE} then the list contains an additional DataFrame for the
#' number of reads and coverage of the input window +/- alignments
#'
#' @importFrom IRanges coverage Views
#' @importFrom S4Vectors DataFrame Rle runValue runLength
#' @importFrom stats pnorm
#'
#' @keywords internal
.keptProbaWin <- function(
    winPosAlignments, winNegAlignments, winWidth, winStep, threshold,
    pvalueThreshold, errorRate, mustKeepWin, minCov, maxCov, getWin,
    useCoverage = FALSE
    )
{
    if (getWin) {
        stopifnot(
            "Coverage" %in% names(winPosAlignments) &&
                "Coverage" %in% names(winNegAlignments)
            )
        fromCoverage <- .calculateStrandCoverage(
            winPosAlignments = winPosAlignments,
            winNegAlignments = winNegAlignments, winWidth = winWidth,
            winStep = winStep
            )
        fromNbReads <- .calculateStrandNbReads(
            winPosAlignments = winPosAlignments,
            winNegAlignments = winNegAlignments
            )
        presentWin <- which(
            as.vector(fromCoverage$CovPos > 0 | fromCoverage$CovNeg > 0) == TRUE
            )
        Win <- DataFrame(
            Type = Rle("",length(presentWin)), Seq = Rle("",length(presentWin)),
            Start = presentWin, End = Rle(0,length(presentWin)),
            NbPos = fromNbReads$NbPos[presentWin],
            NbNeg = fromNbReads$NbNeg[presentWin],
            CovPos = fromCoverage$CovPos[presentWin],
            CovNeg = fromCoverage$CovNeg[presentWin],
            MaxCoverage = fromCoverage$MaxCoverage[presentWin],
            File = Rle("",length(presentWin))
            )
    } else if (useCoverage) {
        fromCoverage <- .calculateStrandCoverage(
            winPosAlignments, winNegAlignments, winWidth, winStep
            )
    } else {
        fromNbReads <- .calculateStrandNbReads(
            winPosAlignments, winNegAlignments
            )
    }
    if (useCoverage) {
        pos <- fromCoverage$CovPos
        neg <- fromCoverage$CovNeg
    } else {
        pos <- fromNbReads$NbPos
        neg <- fromNbReads$NbNeg
    }
    # Logit of the given threshold
    logitThreshold <- log(threshold/(1 - threshold))
    toTest <- (logitThreshold - abs(log(pos/neg)))/sqrt((pos + neg)/pos/neg)
    pvalue <- Rle(pnorm(runValue(toTest)), runLength(toTest))
    rm(toTest)
    pvalue[(pos == 0 | neg == 0)] <- 0
    if (minCov > 0) {
        if (useCoverage) {
            pvalue[fromCoverage$MaxCoverage < minCov] <- 1
        } else {
            pvalue[(pos + neg) < minCov] <- 1
        }
    }
    keptWin <- rep(TRUE, length(runValue(pvalue)))
    keptWin[runValue(pvalue) > pvalueThreshold] <- FALSE
    keptWin <- Rle(keptWin, runLength(pvalue))

    keptProbaPosWin <- keptWin * ((pos > neg) * ((pos - neg)/pos) +
        (pos < neg) * errorRate)
    keptProbaNegWin <- keptWin * ((pos < neg) * ((neg - pos)/neg) +
        (pos > neg) * errorRate)
    keptProbaPosWin[pos==0] <- 0
    keptProbaNegWin[neg==0] <- 0

    if (length(mustKeepWin) > 0) {
        if (length(mustKeepWin$Pos) > length(keptProbaPosWin)) {
            mustKeepWin$Pos <- mustKeepWin$Pos[length(keptProbaPosWin)]
        } else {
            mustKeepWin$Pos <- c(
                mustKeepWin$Pos,
                rep(0, length(keptProbaPosWin) - length(mustKeepWin$Pos))
                )
        }
        keptProbaPosWin <- mustKeepWin$Pos +
            (!mustKeepWin$Pos) * keptProbaPosWin
        if (length(mustKeepWin$Neg) > length(keptProbaNegWin)) {
            mustKeepWin$Neg <- mustKeepWin$Neg[length(keptProbaNegWin)]
        } else {
            mustKeepWin$Neg <- c(
                mustKeepWin$Neg,
                rep(0, length(keptProbaNegWin) - length(mustKeepWin$Neg))
                )
        }
        keptProbaNegWin <- mustKeepWin$Neg +
            (!mustKeepWin$Neg) * keptProbaNegWin
    }
    if (maxCov > 0) {
        if (useCoverage) {
            keepMorePos <- (pos > neg) * (fromCoverage$MaxCoverage >= maxCov)
            keepMoreNeg <- (pos < neg) * (fromCoverage$MaxCoverage >= maxCov)
        } else {
            keepMorePos <- (pos > neg) * (pos >= maxCov)
            keepMoreNeg <- (pos < neg) * (neg >= maxCov)
        }
        keptProbaPosWin <- keepMorePos + (!keepMorePos) * keptProbaPosWin
        keptProbaNegWin <- keepMoreNeg + (!keepMoreNeg) * keptProbaNegWin
    }
    if (getWin) {
        return(list(Pos = keptProbaPosWin, Neg = keptProbaNegWin, Win = Win))
    } else {
        return(list(Pos = keptProbaPosWin, Neg = keptProbaNegWin))
    }
}


