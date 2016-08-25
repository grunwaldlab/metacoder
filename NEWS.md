# News 

## metacoder 0.1.2

### Breaking changes

* `plot_taxonomy` and the `plot` method have been renamed `heat_tree`.

### New features

* New introduction vignette
* Various minor bug fixes


## metacoder 0.1.1

### Breaking changes

* `taxon_levels` have been replaced with `n_supertaxa` to make names conceptually consistent. Note that this means what was `1` as `taxon_levels` is now `0` as `n_supertaxa`.

### New features

* Added `n_subtaxa` and `n_subtaxa_1` functions
* Added taxonomy parsing examples to vignettes


## metacoder 0.1.0

### Breaking changes

* Many options and functions have been renamed (#115)

### New features

* dplyr functions for `taxmap` objects!
* Added a `print` method for `taxmap` objects
* new SILVA example data set
* `extract_taxonomy` works on `SeqFastadna` class from `seqinr`
* `parse_mothur_summary` function: parses the mothur summary table
* `remove_redundant_names` function: removes components of names of taxa in subtaxa

### Changes

* Core functions are much faster
* More tests 
* Updated vignettes
* Many bug fixes and minor upgrades
* Legend now moves into plot if there is room (#118)
