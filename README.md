[![Build Status](https://travis-ci.org/sjessa/chromswitch.svg?branch=master)](https://travis-ci.org/sjessa/chromswitch)
[![Coverage Status](https://coveralls.io/repos/github/sjessa/chromswitch/badge.svg)](https://coveralls.io/github/sjessa/chromswitch)

# `chromswitch`: An R package for detecting chromatin state switches

`chromswitch` implements a flexible method to detect chromatin state 
switches between samples in two biological conditions in a specific genomic
region of interest given peaks called from ChIP-seq data.

# Quickstart

Install `chromswitch` using `devtools`:

```
devtools::install_github("sjessa/chromswitch", build_vignettes = TRUE)
library(chromswitch)
```

Open the vignette for an introduction to the package:

```
browseVignettes("chromswitch")
```
