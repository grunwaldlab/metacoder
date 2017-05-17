[![Build Status](https://travis-ci.org/grunwaldlab/metacoder.png?branch=master)](https://travis-ci.org/grunwaldlab/metacoder?branch=master) [![codecov.io](https://codecov.io/github/grunwaldlab/metacoder/coverage.svg?branch=master)](https://codecov.io/github/grunwaldlab/metacoder?branch=master)
[![Downloads from Rstudio mirror per month](http://cranlogs.r-pkg.org/badges/metacoder)](http://www.r-pkg.org/pkg/metacoder)
[![Downloads from Rstudio mirror](http://cranlogs.r-pkg.org/badges/grand-total/metacoder)](http://www.r-pkg.org/pkg/metacoder)
[![CRAN version](http://www.r-pkg.org/badges/version/metacoder)](http://cran.r-project.org/package=metacoder)

![](readme_figure.png)

An R package for metabarcoding research planning and analysis
-------------------------------------------------------------

Metabarcoding is revolutionizing microbial ecology and presenting new challenges:

-   Numerous database formats make taxonomic data difficult to parse, combine, and subset.
-   Stacked bar charts, commonly used to depict community diversity, lack taxonomic context.
-   Barcode loci and primers are a source of under-explored bias.

Metacoder is an R package that attempts to addresses these issues:

-   Sources of taxonomic data can be extracted from most file formats and manipulated.
-   Community diversity can be visualized by color and size in a tree plot.
-   Primer specificity can be estimated with *in silico* PCR.

### Documentation

Documentation is available at <http://grunwaldlab.github.io/metacoder_documentation>.

### Download the current version

Stable releases are available on CRAN and can be installed in the standard way:

    install.packages("metacoder")

The most recent version can be installed from Github:

    devtools::install_github("grunwaldlab/metacoder@dev")
    library(metacoder)
    

### Dependencies

The function that runs *in silico* PCR requires `primersearch` from the EMBOSS tool kit to be installed. This is not an R package, so it is not automatically installed. Type `?primersearch` after installing and loading MetcodeR for installation instructions.

### Citation

If you use metcoder in a publication, please cite our [article in PLOS Computational Biology](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005404):

Foster ZSL, Sharpton TJ, Grünwald NJ (2017) Metacoder: An R package for visualization and manipulation of community taxonomic diversity data. PLOS Computational Biology 13(2): e1005404. https://doi.org/10.1371/journal.pcbi.1005404

### License

This work is subject to the [MIT License](https://github.com/grunwaldlab/metacoder/blob/master/LICENSE).
