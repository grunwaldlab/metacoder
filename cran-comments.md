## Test environments and check results

### Local computer: Pop!_OS 22.04 LTS, R version 4.4.2

0 errors | 0 warnings | 1 notes

❯ checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: ‘R6’
    All declared Imports should be used.
    
I am not sure why this note exists. I am using `R6::` in about 10 places in the package.
All of the classes in the package are R6 classes, so it is definitely used.



### Winbuilder

0 errors | 0 warnings | 1 notes

```
* checking CRAN incoming feasibility ... [16s] NOTE
Maintainer: 'Zachary Foster <zacharyfoster1989@gmail.com>'

New submission

Package was archived on CRAN

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2024-10-01 as issues were not corrected
    in time.

Suggests or Enhances not in mainstream repositories:
  traits
```

`metacoder` was removed from CRAN because it dependency `taxize` was removed from CRAN.
`taxize` is back on CRAN now, so `metacoder` is being resubmitted.
The package `traits` was also removed from CRAN because it depends on `taxize` as well.
It is only a suggested package for `metacoder` used for one function and I expect `traits` will be returned to CRAN as well soon.
I am not sure if it is required for all packages in "Suggests" to be on CRAN or not.


### Response to CRAN review

I added Scott Chamberlain and Kamil Slowikowski to authors in the DESCRIPTION.