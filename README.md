


[![Build Status](https://travis-ci.org/grunwaldlab/metacoder.png?branch=master)](https://travis-ci.org/grunwaldlab/metacoder?branch=master)





## MetacodeR: an R package for metabarcoding research planning and analysis

*Zachary S. L. Foster and Niklaus J. Grünwald*

### Introduction 

Metabarcoding is revolutionizing microbial ecology and presenting new challenges:

* Numerous database formats make taxonomic data difficult to parse, combine, and subset.
* Stacked bar charts, commonly used to depict community diversity, lack taxonomic context.
* Barcode loci and primers are a source of under-explored bias.

MetacodeR is an R package that attempts to addresses these issues:

* Sources of taxonomic data can be extracted from any file format and manipulated. 
* Community diversity can be visualized by color and size in a tree plot we call a **Metadiversity plot**.
* Primer specificity can be estimated with *In silico* PCR.

### Extracting taxonomic data

Most databases have a unique file format and taxonomic hierarchy/nomenclature.
Taxonomic data can be extracted from any file format using the **extract_taxonomy** function.
Classifications can be parsed offline or retrieved from online databases if a taxon name, taxon ID, or sequence ID is present.
A regular expression with capture groups and a corresponding key is used to define how to parse the file.
The example code below parses the 16s Ribosome Database Project training set for Mothur.
R can be used to download files from the internet and decompress them.
The code below downloads the compressed data to a temporary directory:


```r
rdp_fasta_url <- "http://mothur.org/w/images/b/b5/Trainset10_082014.rdp.tgz"
temp_dir_path <- tempdir()
local_file_path <- file.path(temp_dir_path, basename(rdp_fasta_url))
download.file(url = rdp_fasta_url, destfile = local_file_path, quiet = TRUE)
```

Next we will uncompress the archive and identify the fasta file.


```r
unpacked_file_paths <- untar(local_file_path, list = TRUE) # get contents
untar(local_file_path, exdir = temp_dir_path)
unpacked_fasta_path <- file.path(temp_dir_path, 
                                  unpacked_file_paths[grepl("fasta$", unpacked_file_paths)])
```

The file can then be parsed using the `ape` package and the taxonomy data in the headers can be extracted by `extract_taxonomy`:


```r
library(metacoder)
seqs <- ape::read.FASTA(unpacked_fasta_path)
cat(names(seqs)[1]) # print an example of the sequence headers
```

```
## AB294171_S001198039	Root;Bacteria;Firmicutes;Bacilli;Lactobacillales;Carnobacteriaceae;Alkalibacterium
```

```r
data <- extract_taxonomy(seqs, regex = "^(.*)\\t(.*)",
                         key = c("item_info", "class_name"),
                         class_tax_sep = ";", database = "none")
```

Not, that this command will take a while to process. The resulting object contains sequence information associated with an inferred taxonomic hierarchy:


```r
head(taxon_data(data)) # print a few rows of the taxon data
```

```
##   taxon_id parent_id item_count rank              name
## 1        1      <NA>      10650    1              Root
## 2        2         1      10240    2          Bacteria
## 3        3         2       1956    3        Firmicutes
## 4        4         3       1214    4           Bacilli
## 5        5         4        411    5   Lactobacillales
## 6        6         5         39    6 Carnobacteriaceae
```


### Metadiversity plots

The hierarchical nature of taxonomic data makes it difficult to plot effectively.
Most often, bar charts, stacked bar charts, or pie graphs are used, but these are ineffective when plotting many taxa or multiple ranks.
MetacodeR maps taxonomic data (e.g. sequence abundance) to color/size of tree components in what we call a **Metadiversity plot**:




```r
plot(data, vertex_size = item_count, vertex_label = name, vertex_color = item_count)
```

![](README_files/figure-html/unnamed-chunk-9-1.png)




The default size range displayed is optimized for each plot.
Only a few options are needed to make effective plots, yet many are available for customization of publication-ready graphics:




```r
plot(data, vertex_size = item_count, edge_color = rank,
     vertex_label = name, vertex_color = item_count,
     vertex_color_range = c("cyan", "magenta", "green"),
     edge_color_range   = c("#555555", "#EEEEEE"),
     layout = "davidson-harel", overlap_avoidance = 0.5)
```

![](README_files/figure-html/unnamed-chunk-12-1.png)




### Subsetting

Taxonomic data can be easily subset using the **subset** function.
The user can choose preserve or remove the subtaxa, supertaxa, and sequence data of the subset.
For example, **subset** can be used to look at just the Archaea:




```r
plot(subset(data, name == "Archaea"), vertex_size = item_count, 
     vertex_label = name, vertex_color = item_count, layout = "fruchterman-reingold")
```

![](README_files/figure-html/unnamed-chunk-15-1.png)




To make the Archaea-Bacteria division more clear, the "Root" taxon can be removed, causing two trees to be plotted:




```r
subsetted <- subset(data, rank > 1)
plot(subsetted, vertex_size = item_count, vertex_label = name,
     vertex_color = item_count, tree_label = name, layout = "davidson-harel")
```

![](README_files/figure-html/unnamed-chunk-18-1.png)



### Taxonomically balanced sampling

When calculating statistics for taxa, the amount of data should be balanced among taxa and there should be enough data per taxon to make unbiased estimates.
Random samples from large reference databases are biased toward overrepresented taxa.
The function **taxonomic_sample** is used to create taxonomically balanced random samples.
The acceptable range of sequence or subtaxa counts can be defined for each taxonomic rank; taxa with too few are excluded and taxa with too many are randomly sampled.
The code below samples the data such that rank 6 taxa will have 5 sequences and rank 3 taxa (phyla) will have less than 100:




```r
sampled <- taxonomic_sample(subsetted, max_counts = c("3" = 100, "6" = 5), min_counts = c("6" = 5))
sampled <- subset(sampled, item_count > 0, itemless = FALSE) 
```




```r
plot(sampled, vertex_size = item_count, vertex_label = item_count, overlap_avoidance = 0.5,
     vertex_color = item_count, layout = "davidson-harel")
```

![](README_files/figure-html/unnamed-chunk-23-1.png)




### In silico PCR

The function **primersearch** is a wrapper for an EMBOSS tool that implements *in silico* PCR.
The code below estimates the coverage of the universal bacterial primer pair 357F/519F: 


```r
pcr <- primersearch(sampled, forward = "CTCCTACGGGAGGCAGCAG", reverse = "GWATTACCGCGGCKGCTG",
                    pair_name = "357F_519R",  mismatch = 10)
head(taxon_data(pcr))
```

```
##   taxon_id parent_id item_count rank count_amplified prop_amplified              name
## 2        2      <NA>       1881    1            1658      0.8814460          Bacteria
## 3        3         2        230    2             221      0.9608696        Firmicutes
## 4        4         3        100    3              97      0.9700000           Bacilli
## 5        5         4         38    4              37      0.9736842   Lactobacillales
## 6        6         5          3    5               3      1.0000000 Carnobacteriaceae
## 7        7         6          1    6               1      1.0000000   Alkalibacterium
```

The proportion of sequences amplified can be represented by color in a metadiversity plot:




```r
plot(pcr, vertex_size = item_count, vertex_label = name, vertex_color = prop_amplified,
     vertex_color_range =  c("red", "cyan"), vertex_color_trans = "radius", tree_label = name)
```

![](README_files/figure-html/unnamed-chunk-27-1.png)



This plot makes it apparent that no Archaea were amplified and most Bacteria were amplified, but not all.
The data can also be subset to better see what did not get amplified:




```r
library(magrittr) # Adds optional %>% operator for chaining commands
pcr %>%
  subset(name == "Bacteria") %>%
  subset(count_amplified < item_count, subtaxa = FALSE) %>% 
  plot(vertex_size = item_count, vertex_label = name, vertex_color = prop_amplified,
       vertex_color_range =  c("red", "cyan"),
       vertex_color_interval = c(0, 1), vertex_color_trans = "radius")
```

![](README_files/figure-html/unnamed-chunk-30-1.png)




### Differential metadiversity plots

We can compare the effectiveness of two primer pairs, 357F/519F and 515F/1100R, by plotting the difference in proportions amplified by each.
First, the same sequences are amplified with 515F/1100R and results for the two primer pairs combined:




```r
pcr_2 <- primersearch(sampled, forward = "GTGCCAGCMGCCGCGGTAA", reverse = "AGGGTTGCGCTCGTTG",
                      pair_name = "515F_1100R", mismatch = 10)
pcr$taxon_data$count_amplified_2 <- taxon_data(pcr_2, "count_amplified")
pcr$taxon_data$prop_diff <- taxon_data(pcr, "prop_amplified") - taxon_data(pcr_2, "prop_amplified")
```

Then, taxa that are not amplified by both pairs can be subset and the difference in amplification plotted.
In the plot below, blue corresponds to taxa amplified by 357F/519F but not 515F/1100R and brown is the opposite:





```r
pcr %>%
  subset(name == "Bacteria") %>%
  subset(count_amplified < item_count | count_amplified_2 < item_count, subtaxa = FALSE) %>%
  plot(vertex_size = item_count, vertex_label = name,
       vertex_color = prop_diff, vertex_color_range = diverging_palette(),
       vertex_color_interval = c(-1, 1), vertex_color_trans = "radius")
```

![](README_files/figure-html/unnamed-chunk-34-1.png)




### Conclusions

Intuitive plotting will reveal subtle patterns in complex data and the ability to estimate the effectiveness of new primers/barcodes will accelerate the adoption of metabarcoding to understudied groups of organisms.

### Availability

This package is currently being developed and can be installed from GitHub using the following code:


```r
devtools::install_github("grunwaldlab/metacoder")
```

### Aknowledgements

We thank Tom Sharpton for sharing his metagenomics expertise and advising us.
MetacodeR's major dependencies are taxize, igraph, and ggplot2.

### Documentation

Documentation is under construction at http://grunwaldlab.github.io/metacoder.

### Download the current version

While this project is in development it can be installed through github:

    devtools::install_github(repo="grunwaldlab/metacoder")
    library(metacoder)

If you've built the vignettes, you can browse them with:

    browseVignettes(package="metacoder")
