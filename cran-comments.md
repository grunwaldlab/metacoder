## Test environments

* local desktop: ubuntu 16.04 install, R 3.4.0
* on travis-ci: ubuntu 12.04.5 + R 3.4.0, devel [2017-05-22 r72718]
* win-builder

## R CMD check results

### local desktop:

There were no ERRORs, WARNINGs or NOTEs.

### travis-ci:

There were no ERRORs, WARNINGs or NOTEs.

### On winbuilder:

```
* checking whether package 'metacoder' can be installed ... WARNING
Found the following significant warnings:
  Warning: Installed Rcpp (0.12.11) different from Rcpp used to build dplyr (0.12.10).
  Warning: Installed R (R Under development (unstable) (2017-05-20 r72713)) different from R used to build dplyr (R Under development (unstable) (2017-05-20 r72708)).
See 'd:/RCompile/CRANguest/R-devel/metacoder.Rcheck/00install.out' for details.
```

I think this is because dplyr 0.6 will be released to CRAN soon, but not yet.
This version of metacoder is compatible with both the old and new versions of dplyr.
One of the major reasons for this update is to be compatible with dplyr 0.6.


## Downstream dependencies

There are currently no downstream dependencies for this package.

