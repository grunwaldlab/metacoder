## Test environments and check results

Some packages seem to not be avaialble for R 4.0.0. I am not sure what to do about this. It works on my computer and winbuilder.

### local: ubuntu 18.04, R 3.6.3

0 errors | 0 warnings | 0 notes

### travis-ci: ubuntu 14.04.05, R development

Warning message:

package ‘BiocManager, zlibbioc’ is not available (for R Under development) 

The command "Rscript -e 'install.packages(c("BiocManager, zlibbioc"));if (!all(c("BiocManager, zlibbioc") %in% installed.packages())) { q(status = 1, save = "no")}'" failed and exited with 1 during .

### travis-ci: ubuntu 14.04.05, bioc-release

Warning message:

package ‘BiocManager, zlibbioc’ is not available (for R version 4.0.0) 

The command "Rscript -e 'install.packages(c("BiocManager, zlibbioc"));if (!all(c("BiocManager, zlibbioc") %in% installed.packages())) { q(status = 1, save = "no")}'" failed and exited with 1 during .


### Rhub: Windows Server 2008 R2 SP1, R-devel, 32/64 bit

* checking package dependencies ... ERROR
Packages required but not available:
  'taxa', 'stringr', 'ggplot2', 'igraph', 'scales', 'taxize', 'seqinr',
  'reshape2', 'zoo', 'traits', 'RCurl', 'ape', 'reshape', 'lazyeval',
  'dplyr', 'readr', 'rlang', 'biomformat', 'ggfittext', 'vegan',
  'ggrepel', 'cowplot', 'GA', 'Rcpp', 'svglite', 'tibble'

### Win builder: R devel

0 errors | 0 warnings | 0 notes

### Win builder: R version 4.0.0 (2020-04-24)

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are currently no downstream dependencies for this package.

