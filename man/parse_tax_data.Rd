% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/old_taxa--taxmap--parsers.R
\name{parse_tax_data}
\alias{parse_tax_data}
\title{Convert one or more data sets to taxmap}
\usage{
parse_tax_data(
  tax_data,
  datasets = list(),
  class_cols = 1,
  class_sep = ";",
  sep_is_regex = FALSE,
  class_key = "taxon_name",
  class_regex = "(.*)",
  class_reversed = FALSE,
  include_match = TRUE,
  mappings = c(),
  include_tax_data = TRUE,
  named_by_rank = FALSE
)
}
\arguments{
\item{tax_data}{A table, list, or vector that contains the names of taxa that represent
[taxonomic
classifications](https://en.wikipedia.org/wiki/Taxonomy_(biology)#Classifying_organisms).
Accepted representations of classifications include: * A list/vector or table with column(s) of
taxon names: Something like `"Animalia;Chordata;Mammalia;Primates;Hominidae;Homo"`. What
separator(s) is used (";" in this example) can be changed with the `class_sep` option. For
tables, the classification can be spread over multiple columns and the separator(s) will be
applied to each column, although each column could just be single taxon names with no
separator. Use the `class_cols` option to specify which columns have taxon names. * A list in
which each entry is a classifications. For example, `list(c("Animalia", "Chordata", "Mammalia",
"Primates", "Hominidae", "Homo"), ...)`. * A list of data.frames where each represents a
classification with one taxon per row. The column that contains taxon names is specified using
the `class_cols` option. In this instance, it only makes sense to specify a single column.}

\item{datasets}{Additional lists/vectors/tables that should be included in the resulting `taxmap`
object. The `mappings` option is use to specify how these data sets relate to the `tax_data`
and, by inference, what taxa apply to each item.}

\item{class_cols}{(`character` or `integer`) The names or indexes of columns that contain
classifications if the first input is a table. If multiple columns are specified, they will be
combined in the order given. Negative column indexes mean "every column besides these columns".}

\item{class_sep}{(`character`) One or more separators that delineate taxon names in a
classification. For example, if one column had `"Homo sapiens"` and another had
`"Animalia;Chordata;Mammalia;Primates;Hominidae"`, then `class_sep = c(" ", ";")`. All
separators are applied to each column so order does not matter.}

\item{sep_is_regex}{(`TRUE`/`FALSE`) Whether or not `class_sep` should be used as a [regular
expression](https://en.wikipedia.org/wiki/Regular_expression).}

\item{class_key}{(`character` of length 1) The identity of the capturing groups defined using
`class_regex`. The length of `class_key` must be equal to the number of capturing groups
specified in `class_regex`. Any names added to the terms will be used as column names in the
output. At least one `"taxon_name"` must be specified. Only `"info"` can be used multiple
times. Each term must be one of those described below: * `taxon_name`: The name of a taxon. Not
necessarily unique, but are interpretable by a particular `database`. Requires an internet
connection. * `taxon_rank`: The rank of the taxon. This will be used to add rank info into the
output object that can be accessed by `out$taxon_ranks()`. * `info`: Arbitrary taxon info you
want included in the output. Can be used more than once.}

\item{class_regex}{(`character` of length 1) A regular expression with capturing groups
indicating the locations of data for each taxon in the `class` term in the `key` argument. The
identity of the information must be specified using the `class_key` argument. The `class_sep`
option can be used to split the classification into data for each taxon before matching. If
`class_sep` is `NULL`, each match of `class_regex` defines a taxon in the classification.}

\item{class_reversed}{If `TRUE`, then classifications go from specific to general. For example:
`Abditomys latidens : Muridae : Rodentia : Mammalia : Chordata`.}

\item{include_match}{(`logical` of length 1) If `TRUE`, include the part of the input matched by
`class_regex` in the output object.}

\item{mappings}{(named `character`) This defines how the taxonomic information in `tax_data`
applies to data set in `datasets`. This option should have the same number of inputs as
`datasets`, with values corresponding to each data set. The names of the character vector
specify what information in `tax_data` is shared with info in each `dataset`, which is
specified by the corresponding values of the character vector. If there are no shared
variables, you can add `NA` as a placeholder, but you could just leave that data out since it
is not benefiting from being in the taxmap object. The names/values can be one of the
following: * For tables, the names of columns can be used. * `"\{\{index\}\}"` : This means to use
the index of rows/items * `"\{\{name\}\}"`  : This means to use row/item names. * `"\{\{value\}\}"` :
This means to use the values in vectors or lists. Lists will be converted to vectors using
[unlist()].}

\item{include_tax_data}{(`TRUE`/`FALSE`) Whether or not to include `tax_data` as a dataset, like
those in `datasets`.}

\item{named_by_rank}{(`TRUE`/`FALSE`) If  `TRUE` and the input is a table with columns named by
ranks or a list of vectors with each vector named by ranks, include that rank info in the
output object, so it can be accessed by `out$taxon_ranks()`. If `TRUE`, taxa with different
ranks, but the same name and location in the taxonomy, will be considered different taxa.
Cannot be used with the `sep`, `class_regex`, or `class_key` options.}
}
\description{
Reads taxonomic information and associated data in tables, lists, and vectors and stores it in a
[taxmap()] object. [Taxonomic
classifications](https://en.wikipedia.org/wiki/Taxonomy_(biology)#Classifying_organisms) must be
present.
}
\examples{
 # Read a vector of classifications
 my_taxa <- c("Mammalia;Carnivora;Felidae",
              "Mammalia;Carnivora;Felidae",
              "Mammalia;Carnivora;Ursidae")
 parse_tax_data(my_taxa, class_sep = ";")

 # Read a list of classifications
 my_taxa <- list("Mammalia;Carnivora;Felidae",
                "Mammalia;Carnivora;Felidae",
                "Mammalia;Carnivora;Ursidae")
 parse_tax_data(my_taxa, class_sep = ";")

 # Read classifications in a table in a single column
 species_data <- data.frame(tax = c("Mammalia;Carnivora;Felidae",
                                    "Mammalia;Carnivora;Felidae",
                                    "Mammalia;Carnivora;Ursidae"),
                           species_id = c("A", "B", "C"))
 parse_tax_data(species_data, class_sep = ";", class_cols = "tax")

 # Read classifications in a table in multiple columns
 species_data <- data.frame(lineage = c("Mammalia;Carnivora;Felidae",
                                        "Mammalia;Carnivora;Felidae",
                                        "Mammalia;Carnivora;Ursidae"),
                            species = c("Panthera leo",
                                        "Panthera tigris",
                                        "Ursus americanus"),
                            species_id = c("A", "B", "C"))
 parse_tax_data(species_data, class_sep = c(" ", ";"),
                class_cols = c("lineage", "species"))

 # Read classification tables with one column per rank
 species_data <- data.frame(class = c("Mammalia", "Mammalia", "Mammalia"),
                            order = c("Carnivora", "Carnivora", "Carnivora"),
                            family = c("Felidae", "Felidae", "Ursidae"),
                            genus = c("Panthera", "Panthera", "Ursus"),
                            species = c("leo", "tigris", "americanus"),
                            species_id = c("A", "B", "C"))
  parse_tax_data(species_data, class_cols = 1:5)
  parse_tax_data(species_data, class_cols = 1:5,
                 named_by_rank = TRUE) # makes `taxon_ranks()` work

 # Classifications with extra information
 my_taxa <- c("Mammalia_class_1;Carnivora_order_2;Felidae_genus_3",
              "Mammalia_class_1;Carnivora_order_2;Felidae_genus_3",
              "Mammalia_class_1;Carnivora_order_2;Ursidae_genus_3")
 parse_tax_data(my_taxa, class_sep = ";",
                class_regex = "(.+)_(.+)_([0-9]+)",
                class_key = c(my_name = "taxon_name",
                              a_rank = "taxon_rank",
                              some_num = "info"))


  # --- Parsing multiple datasets at once (advanced) ---
  # The rest is one example for how to classify multiple datasets at once.

  # Make example data with taxonomic classifications
  species_data <- data.frame(tax = c("Mammalia;Carnivora;Felidae",
                                     "Mammalia;Carnivora;Felidae",
                                     "Mammalia;Carnivora;Ursidae"),
                             species = c("Panthera leo",
                                         "Panthera tigris",
                                         "Ursus americanus"),
                             species_id = c("A", "B", "C"))

  # Make example data associated with the taxonomic data
  # Note how this does not contain classifications, but
  # does have a varaible in common with "species_data" ("id" = "species_id")
  abundance <- data.frame(id = c("A", "B", "C", "A", "B", "C"),
                          sample_id = c(1, 1, 1, 2, 2, 2),
                          counts = c(23, 4, 3, 34, 5, 13))

  # Make another related data set named by species id
  common_names <- c(A = "Lion", B = "Tiger", C = "Bear", "Oh my!")

  # Make another related data set with no names
  foods <- list(c("ungulates", "boar"),
                c("ungulates", "boar"),
                c("salmon", "fruit", "nuts"))

  # Make a taxmap object with these three datasets
  x = parse_tax_data(species_data,
                     datasets = list(counts = abundance,
                                     my_names = common_names,
                                     foods = foods),
                     mappings = c("species_id" = "id",
                                  "species_id" = "{{name}}",
                                  "{{index}}" = "{{index}}"),
                     class_cols = c("tax", "species"),
                     class_sep = c(" ", ";"))

  # Note how all the datasets have taxon ids now
  x$data

  # This allows for complex mappings between variables that other functions use
  map_data(x, my_names, foods)
  map_data(x, counts, my_names)

}
\seealso{
Other parsers: 
\code{\link{extract_tax_data}()},
\code{\link{lookup_tax_data}()},
\code{\link{parse_dada2}()},
\code{\link{parse_edge_list}()},
\code{\link{parse_greengenes}()},
\code{\link{parse_mothur_tax_summary}()},
\code{\link{parse_mothur_taxonomy}()},
\code{\link{parse_newick}()},
\code{\link{parse_phylo}()},
\code{\link{parse_phyloseq}()},
\code{\link{parse_qiime_biom}()},
\code{\link{parse_rdp}()},
\code{\link{parse_silva_fasta}()},
\code{\link{parse_ubiome}()},
\code{\link{parse_unite_general}()}
}
\concept{parsers}
