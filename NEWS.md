# News 

## metacoder 0.1.0

### Breaking changes

* Many options and functions have been renamed (#115)

### New features

* dplyr functions for `taxmap` objects!
* Added a `print` method for `taxmap` objects
* new SILVA example dataset
* `extract_taxonomy` works on `SeqFastadna` class from `seqinr`
* `parse_mothur_summary` function: parses the mothur summary table
* `remove_redundant_names` function: removes components of names of taxa in subtaxa

### Changes

* Core functions are much faster
* More tests 
* Updated vignettes
* Many bug fixs and minor upgrades
* Legend now moves into plot if there is room (#118)
