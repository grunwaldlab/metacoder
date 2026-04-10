An identical version was submitted a few months ago, but `traits` was not on CRAN at the time, which complicated things.
Now, `traits` is on CRAN again, so all issues broght up by the previous review should be addressed.

## Test environments and check results

### Local computer: Pop!_OS 22.04 LTS, R version 4.5.2

0 errors | 0 warnings | 1 notes

❯ checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: ‘R6’
    All declared Imports should be used.
    
I am not sure why this note exists. I am using `R6::` in about 10 places in the package.
All of the classes in the package are R6 classes, so it is definitely used.

### R-hub, and CI results:

https://github.com/grunwaldlab/metacoder/actions

### Winbuilder

0 errors | 0 warnings | 1 notes

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Zachary Foster <zacharyfoster1989@gmail.com>'

New submission

Package was archived on CRAN

Possibly misspelled words in DESCRIPTION:
  al (32:18)
  bioinformatics (30:16)
  et (32:15)
  metabarcoding (26:25)
  metagenomics (26:59, 27:5)
  microbiome (30:5)

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2025-10-18 as issues were not corrected
    despite reminders.
