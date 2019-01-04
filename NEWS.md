# News 

##  metacoder 0.3.1

### Bug fixes

* `calc_taxon_abund` no longer errors when a taxon has no observations associated with it.
* New `heat_tree_matrix` options to change size and color of row and column labels. 
* Fixed a bug causing the size legend not to be shown (issue [#249](https://github.com/grunwaldlab/metacoder/issues/249).
* Now when a node_color_interval is set but a edge_color_interval is not and edge_color is not used, the edge_color_interval is the same as the node_color_interval.

### New features

* Added parser for `dada2` results called `parse_dada2` and writers to convert back called `make_dada2_asv_table` and `make_dada2_tax_table`.

### Changes

* Started using the `viridis` colors for `heat_tree` by default (issue [#133](https://github.com/grunwaldlab/metacoder/issues/133).

##  metacoder 0.3.0

### Bug fixes

* Fixed bug in `calc_n_samples` where the message reported the number of taxa instead of the number of rows in the table.
* Fixed bug in `heat_tree_matrix` that happened when factors were used for treatments (issue [#240](https://github.com/grunwaldlab/metacoder/issues/240).
* `zero_low_counts` now ignores `NA`s instead of odd error.
* `compare_groups` now ignores `NA`s instaed of returning `NaN`

### Improvements

* Added `more_than` option to `calc_n_samples` so that users can set the minimum threshold for whether a sample is counted or not instead of it always 1.
* Added `calc_prop_samples` function for calculating the proportion of samples with a value greater than 0 (issues [#233](https://github.com/grunwaldlab/metacoder/issues/233).
* `primersearch` is faster and takes less memory by using `ape::DNAbin` objects internally.
* Made `calc_taxon_abund` about 5x faster.

### New features

* `taxmap` objects can be converted to `phyloseq` objects using `as_phyloseq`.
* Added parser for uBiome data.

### Changes

* `primersearch` now takes and returns a `taxmap` object with results added as tables. `primersearch_raw` is a new function that behaves like the old `primersearch` did, returning a table.
* The `dataset` option of many functions has been renamed to `data` to match the option name in the `taxa` package.
* Numerous spelling fixes.

##  metacoder 0.2.1

### Bug fixes

* Fixes numerous bugs in `heat_tree_matrix` that happen when the input data is not exactly like that produced by `compare_groups` (issues [#195](https://github.com/grunwaldlab/metacoder/issues/195), [#196](https://github.com/grunwaldlab/metacoder/issues/196), [#197](https://github.com/grunwaldlab/metacoder/issues/197)). 
* Fixed how `output_file` was used with `heat_tree_matrix`. Now whole plot is saved instead of last subplot.  (issue [#203](https://github.com/grunwaldlab/metacoder/issues/203))
* Fixed "unused argument" bug in `parse_mothur_tax_summary` when reading from a file path (issue [#211](https://github.com/grunwaldlab/metacoder/issues/211)).
* Fixed bug when in `zero_low_counts` when using `use_total = TRUE` (issue [#227](https://github.com/grunwaldlab/metacoder/issues/227)).
* Numerous other small fixes.
* Fixed `parse_phyloseq` error when arbitrary rank names were used.


### Improvements

* Node and edge legends can now be excluded individually (Thanks [\@grabear](https://github.com/grabear)!) (issue [#202](https://github.com/grunwaldlab/metacoder/issues/202)).
* The output of `heat_tree_matrix` always has a 1:1 aspect ratio. (issue [#205](https://github.com/grunwaldlab/metacoder/issues/205))
* Numerous calculation functions added, with more consistent behavior.

## metacoder 0.2.0

### Bug fixes

* Fixed bug in `subtaxa` that caused an error when all of `subset` is `FALSE`. (issue [#143](https://github.com/grunwaldlab/metacoder/issues/143))
* Fixed bug in `filter_taxa` that caused an error when all taxa are filtered out. (issue [#144](https://github.com/grunwaldlab/metacoder/issues/144))

### Breaking changes

* All taxmap-related manipulation functions have been moved to the [taxa](https://github.com/ropensci/taxa) package.
* `heat_tree` now uses the `taxmap` class defined in the [taxa](https://github.com/ropensci/taxa) package.
* Numerous changes (i.e. upgrades) to `primersearch`

### Improvements

* Upgraded `primersearch` output to be cleaner and have info like the amplicon sequence and primer binding sites.
* Added functions to identift and remove taxa with ambiguous names like "unknown"
* code from [ggrepel](https://github.com/slowkow/ggrepel) package now used to avoid overlapping labels. Thanks Kamil Slowikowski!
* New function `heat_tree_matrix` to make plotting a pairwise matrix of heat trees for comparing treatments.
* New parser named `parse_mothur_tax_summary` for mothur *.tax.summary file made by [classify.seqs](https://www.mothur.org/wiki/Classify.seqs).
* New parser named `parse_mothur_taxonomy` for mothur *.taxonomy file made by [classify.seqs](https://www.mothur.org/wiki/Classify.seqs).
* New parser named `parse_qiime_biom` for the QIIME BIOM output.
* New parser named `parse_phyloseq` to convert phyloseq objects.
* New parser named `parse_newick` to parse newick files.
* New parser named `parse_unite_general` for unite general FASTA release. (issue [#154](https://github.com/grunwaldlab/metacoder/issues/154))
* New parser named `parse_rdp` for RDP FASTA release. (issue [#160](https://github.com/grunwaldlab/metacoder/issues/160))
* New parser named `parse_silva_fasta` for SILVA FASTA release. (issue [#162](https://github.com/grunwaldlab/metacoder/issues/162))
* New function `calc_obs_props` to calculate proportions from observation counts (issue [#167](https://github.com/grunwaldlab/metacoder/issues/167)
* New parser named `parse_greengenes` for the Greengenes database. (issue [#?](https://github.com/grunwaldlab/metacoder/issues/?))
* New writer named `write_greengenes` to create an imitation of the Greengenes database format. 
* New writer named `write_rdp` to create an imitation of the RDP database format. 
* New writer named `write_mothur_taxonomy` to create an imitation of the mothur taxonomy format. 
* New writer named `write_unite_general` to create an imitation of the UNITE general FASTA release. 
* New writer named `write_silva_fasta` to create an imitation of the SILVA FASTA release. 
* New function named `compare_treatments` to compare multiple samples in multiple treatments, applying a user-defined function.
* New function named `calc_taxon_abund` to sum observation values for each taxon.
* Added `col_names` option to `calc_taxon_abund` to set names of output columns.

##  metacoder 0.1.3 

### Improvements

* Provided helpful error message when the `evaluation nested too deeply: infinite recursion / options(expressions=)?` occurs due to too many labels being printed.
* `heat_tree`: improved how the predicted bondries of text is calcuated, so text with any rotation, justification, or newlines influences margins correctly (i.e. does not get cut off).
* `heat_tree`: Can now save multiple file outputs in different formats at once

### Minor changes

* `heat_tree` now gives a warning if infinite values are given to it
* `extract_taxonomy`: There is now a warning message if class regex does not match ([issue #123](https://github.com/grunwaldlab/metacoder/issues/123))
* `heat_tree`: Increased lengend text size and reduced number of labels
* `extract_taxonomy`: added `batch_size` option to help deal with invalid IDs better
* Added CITATION file


### Breaking changes

* The `heat_tree` option `margin_size` funcion now takes four values instead of 2.

### Bug fixes

* `heat_tree`: Fixed bug when color is set explicitly (e.g. "grey") instead of raw numbers and the legend is not removed. Now a mixure of raw numbers and color names can be used. 
* Fixed bugs caused by dplyr version update
* Fixed bug in `heat_tree` that made values not in the input taxmap object not associate with the right taxa. See [this post](https://groups.google.com/d/msgid/metacoder-discussions/c9d8ecc2-1efa-4baf-946e-0f105575da2e%40googlegroups.com).
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
