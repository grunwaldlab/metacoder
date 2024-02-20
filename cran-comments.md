## Test environments and check results


### local: Pop!_OS 22.04 LTS, R 4.3.2 

0 errors | 0 warnings | 1 notes

❯ checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: ‘R6’
    All declared Imports should be used.
    
I am not sure why this note exists. I am using `R6::` in about 10 places in the package.
All of the classes in the package are R6 classes, so it is definitely used.


### WinBuilder:

0 errors | 0 warnings | 0 notes
