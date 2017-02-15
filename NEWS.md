# News 

## Current 

### Improvements

* `heat_tree`: improved how the predicted bondries of text is calcuated, so text with any rotation, justification, or newlines influences margins correctly (i.e. does not get cut off).


### Minor changes

* Increased lengend text size and reduced number of labels
* `extract_taxonomy`: added `batch_size` option to help deal with invalid IDs better

### Breaking changes

* The `heat_tree` option `margin_size` funcion now takes four values instead of 2.

### Bug fixes

* `extract_taxonomy`: Fixed an error that occured when not all inputs could be classified and sequences were supplied
* Fixed bug in `primersearch` that cased the wrong primer sequence to be returned when primers match in the reverse direction
* Fixed a bug in `parse_mothur_summary` where "unclassified" had got changed to "untaxmap" during a search and replace
* Fixed outdated example code for `extract_taxonomy`
* Fixed a bug in `mutate_taxa` and `mutate_obs` that made replacing columns result in new columns with duplicate names. 


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
