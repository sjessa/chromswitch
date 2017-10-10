context("Helper functions for reading peaks")


test_that("MACS2 files are read into peaks correctly", {
    skip_on_os("windows")

    samples <- c("E068", "E071", "E074", "E101", "E102", "E110")
    conditions <- c(rep("Brain", 3), rep("Other", 3))
    bedfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
                            package = "chromswitch")

    metadata <- data.frame(Sample = samples, Condition = conditions)

    expect_equal(readNarrowPeak(bedfiles, metadata), H3K4me3, tolerance = 1e-2)

})
