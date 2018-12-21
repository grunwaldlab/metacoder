## Test environments

* local desktop: ubuntu 16.04 install, R 3.5.1
* on travis-ci: ubuntu 14.04.5 + R 3.5.1, devel, bioc-release
* win-builder

## R CMD check results

### local desktop:

(using `taxa` version 3.2 that is being submitted to CRAN at the same time as this)

There were no ERRORs, WARNINGs or NOTEs.

### travis-ci bioc-release:

There were no ERRORs, WARNINGs or NOTEs.

### On winbuilder:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Zachary Foster <zacharyfoster1989@gmail.com>'

New maintainer:
  Zachary Foster <zacharyfoster1989@gmail.com>
Old maintainer(s):
  ORPHANED
```

I somehow missed the email from Brian Ripley until I searched for it after seeing this NOTE. Sorry about that!

## Downstream dependencies

There are currently no downstream dependencies for this package.

