\name{strandCheckR-package}
\alias{strandCheckR-package}
\alias{strandCheckR}
\docType{package}
\title{
    Quantify and Filter putative double strand DNA from strand-specific RNA 
    bam file
}
\description{
    This package aims to quantify and remove putative double strand DNA from a 
    strand-specific RNA sample. There are also options and methods to plot the 
    positive/negative proportions of all sliding windows, which allow users to 
    have an idea of how much the sample was contaminated and the appropriate 
    threshold to be used for filtering.
}

\details{
The package has some following main functions:

- \code{\link{getStrandFromBamFile}}: calculate positive/negative proprortion and 
sum of reads over all sliding windows from a bam file

- \code{\link{plotHist}}: plot histogram of positive proportion of windows 
calculated from \code{\link{getStrandFromBamFile}} method

- \code{\link{plotWin}}: plot positive proportion vs number of reads of 
windows calculated from \code{\link{getStrandFromBamFile}} method

- \code{\link{filterDNA}}: filter a bam file
}
\author{
Thu-Hien To & Steve Pederson

Maintainer: Thu-Hien To <tothuhien@gmail.com>
}
%\references{
%  This optional section can contain literature or other references for
%  background information.
%}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
%\keyword{ package }
%\seealso{
%  Optional links to other man pages
%}
\examples{
bamfilein <- system.file("extdata","s1.sorted.bam",package = "strandCheckR")
windows <- getStrandFromBamFile(bamfilein)
plotWin(windows)
plotHist(windows)
filterDNA(file = bamfilein,destination = "filter.bam")
}
