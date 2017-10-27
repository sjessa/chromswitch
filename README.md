[![Build Status](https://travis-ci.org/sjessa/chromswitch.svg?branch=master)](https://travis-ci.org/sjessa/chromswitch)
[![Coverage Status](https://coveralls.io/repos/github/sjessa/chromswitch/badge.svg)](https://coveralls.io/github/sjessa/chromswitch)

# `chromswitch`: An R/Bioconductor package for detecting chromatin state switches

`chromswitch` implements a flexible method to detect chromatin state 
switches between samples in two biological conditions in a specific genomic
region of interest given peaks or chromatin state calls from ChIP-seq data.

# Quickstart

`chromswitch` is currently available from the development version of Bioconductor (3.6). To use Bioconductor 3.6, you will need  [BiocInstaller](https://bioconductor.org/packages/release/bioc/html/BiocInstaller.html):

```
source("https://bioconductor.org/biocLite.R")
biocLite("BiocInstaller")
````

Then, install `chromswitch` (dependencies should be installed automatically):

```
library(BiocInstaller)
useDevel() # This will install Bioconductor 3.6
biocLite("chromswitch")
library(chromswitch)
```

For an introduction to the package, check out the [vignette](https://bioconductor.org/packages/devel/bioc/vignettes/chromswitch/inst/doc/chromswitch_intro.html), or open it from the R console:

```
browseVignettes("chromswitch")
```

# Issues

Bug reports and issues are welcome, please [open an issue](https://github.com/sjessa/chromswitch/issues/new) in the GitHub
repository.
