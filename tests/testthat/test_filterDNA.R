file <- system.file("extdata", "ex1.bam", package = "Rsamtools")
destination <- tempfile(fileext = ".bam")
stat <- tempfile(fileext = ".stat")

test_that("filterDNA returns proper DataFrame & output files", {

    win <- filterDNA(file,destination = destination, statFile = stat, getWin = TRUE)
    expect_true(is(win, "DataFrame"))

    require_cols <- c("Seq", "Start", "End", "NbPos", "NbNeg", "CovPos", "CovNeg", "MaxCoverage", "File")
    expect_true(all(require_cols %in% colnames(win)))

    expect_true(file.exists(destination))
    expect_true(file.exists(stat))
    ln <- readLines(stat)
    expect_true(all(grepl("(Sequences|Summary|Number|Removal|Total).+", ln)))

    unlink(destination)
    filterDNA(file,destination = destination, statFile = NULL) # Check the NULL works
    expect_true(file.exists(destination))
})


if (file.exists(destination)) unlink(destination)
if (file.exists(stat)) unlink(stat)
