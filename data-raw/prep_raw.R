# Code to prep raw data for an exported dataset to use for usage examples
# in the `chromswitch` package documentation
#
# Executed from the top-level directory

library(tidyverse)

raw <- list.files(pattern = "*chr19.nPk.bed",
                  full.names = TRUE) %>%
    map(read_tsv,
        col_names = c("chr", "start", "end", "name", "score","strand",
                      "signalValue", "pValue", "qValue", "peak"))

# Filter chr19 peaks for a few Roadmap brain & other samples to a smaller region
tidy <- raw %>%
    map(filter, start >= 54358955) %>%
    map(filter, end <= 55074918)

# Write to BED
samples <- c("brain1", "brain2", "brain3", "other1", "other2", "other3")
outfiles <- paste0("inst/extdata/", samples, ".H3K4me3.bed")

lapply(seq_along(tidy), function(i)
    write_tsv(x = tidy[[i]],
              path = outfiles[i],
              col_names = FALSE))

# Load it up as we would peaks in an actual analysis, and save the dataset
metadata <- data.frame(Sample = samples,
                       H3K4me3 = outfiles,
                       stringsAsFactors = FALSE)

H3K4me3 <- chromswitch::loadBed(metadata, "H3K4me3",
                     metadata_cols = c("name", "score","strand",
                     "signalValue", "pValue", "qValue", "peak"))

devtools::use_data(H3K4me3, overwrite = TRUE)


