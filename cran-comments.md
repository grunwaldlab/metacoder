## Test environments

* local desktop: ubuntu 12.04 install, R 3.2.3
* on travis-ci: ubuntu 12.04, devel [2016-08-21], R 3.3.0
* win-builder (devel [2016-08-21])

## R CMD check results

There were no ERRORs, WARNINGs.

only NOTE: This is a new submission.

## Downstream dependencies

There are currently no downstream dependencies for this package

## Resubmissions

### Downloading test failures for `extract_taxonomy`

This happended last time on one CRAN test computer:

 > test_check("metacoder")
  1. Failure: Exracting by obs_id works (@test--extract_taxonomy.R#27) -----------
  "Eukaryota" %in% result$taxon_data$name isn't true.

I think this is a rare error caused by the NCBI API. 
I have made these types of tests be skipped on CRAN and have made the tests print more information when they do fail, as suggested by Duncan Murdoch.