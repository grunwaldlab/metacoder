# MetacodeR

`MetacodeR` is an R package that provides a set of tools for:

- Evaluating potential metabarcoding primers/loci for taxonomic specificity and discrimination.
- Standardized parsing of taxonomic information.
- Visualization of statistics distributed over taxonomic classifications.

To accomplish these goals, metacoder leverages resources from other R packages, interfaces with external programs, and provides novel functions where needed to allow for entire analyses within R. Using R means that metacoder is entirely open-source, works on the Linux, Windows, and Apple operating systems, and is seamlessly integrated with the best free tools for statistical analysis.

## Download

While this project is in development it can be installed through github:

    devtools::install_github(repo="grunwaldlab/metacoder", build_vignettes=TRUE)
    library(metacoder)
