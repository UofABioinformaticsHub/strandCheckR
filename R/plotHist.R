#' @title Plot the histogram of positive proportions
#'
#' @description Plot the histogram of positive proportions of the input
#' data frame coming from \code{getStrandFromBamFile}
#'
#' @param windows data frame containing the strand information of the sliding
#' windows. Windows can be obtained using the function
#' \code{getStrandFromBamFile}.
#' @param save if TRUE, then the plot will be save into the file given by
#' \code{file} parameter
#' @param file the file name to save to plot
#' @param split an integer vector that specifies how you want to partition the
#' windows based on the coverage. By default \code{split} = c(10,100,1000),
#' which means that your windows will be partitionned into 4 groups, those have
#' coverage < 10, from 10 to 100, from 100 to 1000, and > 1000
#' @param breaks an integer giving the number of bins for the histogram
#' @param useCoverage if TRUE then plot the coverage strand information,
#' otherwise plot the number of reads strand information. FALSE by default
#' @param groupBy the columns that will be used to split the data.
#' @param normalizeBy instead of using the raw read count/coverage, we will
#' normalize it to a proportion by dividing it to the total number of read
#' count/coverage of windows that have the same value in the \code{normalizeBy}
#' columns.
#' @param heatmap if TRUE, then use heat map to plot the histogram, otherwise
#' use barplot. FALSE by default.
#' @param ... used to pass parameters to facet_wrap
#'
#' @return If \code{heatmap=FALSE}: a ggplot object
#'
#' @seealso \code{\link{getStrandFromBamFile}}, \code{\link{plotWin}}
#'
#' @examples
#' bamfilein = system.file('extdata','s1.sorted.bam',package = 'strandCheckR')
#' win  <- getStrandFromBamFile(file = bamfilein,sequences='10')
#' plotHist(win)
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom grid unit
#' @importFrom dplyr filter
#' @importFrom rlang !! sym
#' @export
plotHist <- function(
        windows, save = FALSE, file = "hist.pdf", groupBy = NULL,
        normalizeBy = NULL, split = c(10, 100, 1000), breaks = 100,
        useCoverage = FALSE, heatmap = FALSE, ...
)
{
    histWin <- .summarizeHist(
        windows, split = split, breaks = breaks, useCoverage = useCoverage,
        groupBy = groupBy, normalizeBy = normalizeBy
    )
    # The initial checks for appropriate input
    reqWinCols <- c("PosProp", "ReadCountProp")
    stopifnot(all(reqWinCols %in% colnames(histWin)))
    stopifnot(is.logical(save))
    allows_facet_wrap <- setdiff(colnames(histWin), c(reqWinCols, "Coverage"))
    groupBy <- intersect(groupBy, allows_facet_wrap)

    # Make the plot
    if (heatmap == FALSE) {
        g <- ggplot(
            histWin,
            aes(
                x = !!sym("PosProp"), y = !!sym("ReadCountProp"), fill = !!sym("Coverage")
            )
        ) +
            geom_bar(stat = "identity") +
            labs(
                x = "Proportion of Reads on '+' Strand",
                y = "Proportion of Windows"
            ) +
            theme_bw() +
            theme(plot.margin = unit(c(0.02, 0.04, 0.03, 0.02), "npc"))
        if (length(groupBy) > 0) {
            # Get any facet arguments from dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(facet_wrap))
            keepArgs <- names(dotArgs) %in% setdiff(allowed, "facets")
            argList <- c(list(facets = groupBy), dotArgs[keepArgs])
            myFacets <- do.call(facet_wrap, argList)
            g <- g + myFacets
        }
    } else {
        cov <- unique(histWin$Coverage)
        g <- list()
        if (!("File" %in% colnames(histWin))) {
            histWin$File <- "File"
        }
        groupBy <- groupBy[groupBy != "File"]
        if (length(groupBy) > 0) {
            # Get any facet arguments from dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(facet_wrap))
            keepArgs <- names(dotArgs) %in% setdiff(allowed, "facets")
            argList <- c(list(facets = groupBy), dotArgs[keepArgs])
            myFacets <- do.call(facet_wrap, argList)
        }
        for (i in seq_along(cov)) {
            l <- filter(histWin, Coverage == cov[i])
            g[[i]] <- ggplot(
                l,
                aes(x = !!sym("PosProp"), y = !!sym("File"), fill = !!sym("ReadCountProp"))
            ) +
                geom_tile() + theme_bw() +
                theme(
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()
                ) +
                scale_fill_gradient(
                    low = "white", high = "red", na.value = "white"
                )
            if (length(groupBy) > 0) {
                g[[i]] <- g[[i]] + myFacets
            }
            if (length(cov) > 1) {
                g[[i]] <- g[[i]] + ggtitle(paste0("Coverage ", cov[i]))
            }
        }
        g <- grid.arrange(grobs = g, nrow = length(cov))
    }

    if (save == TRUE) {
        message("The plot will be saved to the file ", file)
        dotArgs <- list(...)
        allowed <- names(formals(ggsave))
        keepArgs <- names(dotArgs) %in% setdiff(allowed, c("filename", "plot"))
        argList <- c(list(filename = file, plot = g), dotArgs[keepArgs])
        do.call(ggsave, argList)
    }
    if (heatmap == FALSE) return(g)
}

