[![Build Status](https://travis-ci.org/sjessa/chromswitch.svg?branch=master)](https://travis-ci.org/sjessa/chromswitch)
[![Coverage Status](https://coveralls.io/repos/github/sjessa/chromswitch/badge.svg)](https://coveralls.io/github/sjessa/chromswitch)

# `chromswitch`: An R/Bioconductor package for detecting chromatin state switches

`chromswitch` implements a flexible method to detect chromatin state 
switches between samples in two biological conditions in a specific genomic
region of interest given peaks or chromatin state calls from ChIP-seq data.

# Method

The `chromswitch` paper is now published _Bioinformatics_:

> Selin Jessa, Claudia L Kleinman (2018). chromswitch: a flexible method to detect chromatin state switches, Bioinformatics, [https://doi.org/10.1093/bioinformatics/bty075](https://doi.org/10.1093/bioinformatics/bty075)

For details on the method, see the [Supplementary Methods](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/PAP/10.1093_bioinformatics_bty075/2/bty075_supp.pdf?Expires=1523314027&Signature=MOHlKYzBdR7bNc64a2G-k5eg2TofADoXEVwFUF94ymNUy0U1CuhnvuY1G2PatkwY7toY2svdhLdnivFQZajRfVrEZrMTVaaFV9qEowwIV3OJUmyUbBX2sFR9Hm17yvkFBwhJYjr-Hjek~Q3wXbGVZ7shO8kTJwgiSCaFFI0gbO3qz3PB80tq1NQllArBuZ01qTjxh1eKHUldljsjAOG2jDYHOVHMHQOZwyyiNkhhwvtE~RbdmpmV4y-AqkHiJB4sy~hIT3PyqD67LyEJMndVoNgCa76MaqUVjhH-l6I53f4kQCNtc7tfD2Q8vshU0whtXZUdNAvHIMPTxVNazsiEjA__&Key-Pair-Id=APKAIUCZBIA4LVPAVW3Q).


# Installation

`chromswitch` is available from Bioconductor at [bioconductor.org/packages/chromswitch](https://bioconductor.org/packages/release/bioc/html/chromswitch.html). To install the package from the R console:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("chromswitch")
````


For an introduction to the package, check out the [vignette](https://bioconductor.org/packages/release/bioc/vignettes/chromswitch/inst/doc/chromswitch_intro.html), or open it from the R console:

```
browseVignettes("chromswitch")
```

# Issues

Bug reports and issues are welcome, please [open an issue](https://github.com/sjessa/chromswitch/issues/new) in the GitHub
repository.
