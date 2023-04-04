## Test environments and check results


### local: Pop!_OS 22.04 LTS, R 4.2.2 

0 errors | 0 warnings | 1 notes

❯ checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: ‘R6’
    All declared Imports should be used.
    
I am not sure why this note exists. I am using `R6::` in about 10 places in the package.
All of the classes in the package are R6 classes, so it is definitely used. 


### Rhub: Windows

0 errors | 0 warnings | 2 notes

* checking dependencies in R code ... NOTE
Namespace in Imports field not imported from: 'R6'
  All declared Imports should be used.
  
I am not sure why this note exists. I am using `R6::` in about 10 places in the package.
All of the classes in the package are R6 classes, so it is definitely used. 

* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'

Seems to be a bug with Rhub:

https://github.com/r-hub/rhub/issues/503



## Downstream dependencies

There are currently no downstream dependencies for this package.

