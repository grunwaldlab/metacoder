[![Build Status](https://travis-ci.org/grunwaldlab/metacoder.png?branch=master)](https://travis-ci.org/grunwaldlab/metacoder?branch=master) [![codecov.io](https://codecov.io/github/grunwaldlab/metacoder/coverage.svg?branch=master)](https://codecov.io/github/grunwaldlab/metacoder?branch=master)

An R package for metabarcoding research planning and analysis
-------------------------------------------------------------

Metabarcoding is revolutionizing microbial ecology and presenting new challenges:

-   Numerous database formats make taxonomic data difficult to parse, combine, and subset.
-   Stacked bar charts, commonly used to depict community diversity, lack taxonomic context.
-   Barcode loci and primers are a source of under-explored bias.

MetacodeR is an R package that attempts to addresses these issues:

-   Sources of taxonomic data can be extracted from any file format and manipulated.
-   Community diversity can be visualized by color and size in a tree plot.
-   Primer specificity can be estimated with *in silico* PCR.

### Documentation

Documentation is under construction at <http://grunwaldlab.github.io/metacoder>.

### Download the current version

While this project is in development it can be installed through Github:

    devtools::install_github(repo="grunwaldlab/metacoder", build_vignettes = TRUE)
    library(metacoder)

If you've built the vignettes, you can browse them with:

    browseVignettes(package="metacoder")

### Dependencies

The function that runs *in silico* PCR requires `primersearch` from the EMBOSS tool kit to be installed. This is not an R package, so it is not automatically installed. Type `?primersearch` after installing and loading MetcodeR for installation instructions.