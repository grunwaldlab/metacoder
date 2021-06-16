#' Convert a table with an edge list to taxmap
#'
#' Converts a table containing an edge list into a [taxa::taxmap()] object.
#' An "edge list" is two columns in a table, where each row defines a taxon-supertaxon relationship.
#' The contents of the edge list will be used as taxon IDs.
#' The whole table will be included as a data set in the output object.
#'
#' @param input A table containing an edge list encoded by two columns.
#' @param taxon_id The name/index of the column containing the taxon IDs.
#' @param supertaxon_id The name/index of the column containing the taxon IDs for the supertaxon of the IDs in `taxon_col`.
#' @param taxon_name xxx
#' @param taxon_rank xxx
#'
#' @family parsers
#'
#' @export
parse_edge_list <- function(input, taxon_id, supertaxon_id, taxon_name, taxon_rank = NULL) {

  # Create empty taxmap object
  output <- taxmap()

  # Make taxon ID characters
  input[taxon_id] <- as.character(input[[taxon_id]])
  input[supertaxon_id] <- as.character(input[[supertaxon_id]])

  # Add edge list
  output$edge_list <- data.frame(from = input[[supertaxon_id]],
                                 to = input[[taxon_id]],
                                 stringsAsFactors = FALSE)

  # Add taxa
  output$taxa <- lapply(seq_len(nrow(input)), function(i) {
    my_name <- input[[taxon_name]][i]
    if (is.null(taxon_rank)) {
      my_rank <- NULL
    } else {
      my_rank <- input[[taxon_rank]][i]
    }
    my_id <- input[[taxon_id]][i]
    taxon(name = my_name, rank = my_rank, id = my_id)
  })
  names(output$taxa) <- input[[taxon_id]]

  # Add data
  input <- dplyr::mutate(input, taxon_id = taxon_ids(output))
  input <- dplyr::select(input, taxon_id, everything())
  output$data <- list(input = input)

  return(output)
}


#' Convert one or more data sets to taxmap
#'
#' Reads taxonomic information and associated data in tables, lists, and vectors and stores it in a
#' [taxa::taxmap()] object. [Taxonomic
#' classifications](https://en.wikipedia.org/wiki/Taxonomy_(biology)#Classifying_organisms) must be
#' present.
#'
#' @param tax_data A table, list, or vector that contains the names of taxa that represent
#'   [taxonomic
#'   classifications](https://en.wikipedia.org/wiki/Taxonomy_(biology)#Classifying_organisms).
#'   Accepted representations of classifications include: * A list/vector or table with column(s) of
#'   taxon names: Something like `"Animalia;Chordata;Mammalia;Primates;Hominidae;Homo"`. What
#'   separator(s) is used (";" in this example) can be changed with the `class_sep` option. For
#'   tables, the classification can be spread over multiple columns and the separator(s) will be
#'   applied to each column, although each column could just be single taxon names with no
#'   separator. Use the `class_cols` option to specify which columns have taxon names. * A list in
#'   which each entry is a classifications. For example, `list(c("Animalia", "Chordata", "Mammalia",
#'   "Primates", "Hominidae", "Homo"), ...)`. * A list of data.frames where each represents a
#'   classification with one taxon per row. The column that contains taxon names is specified using
#'   the `class_cols` option. In this instance, it only makes sense to specify a single column.
#' @param datasets Additional lists/vectors/tables that should be included in the resulting `taxmap`
#'   object. The `mappings` option is use to specify how these data sets relate to the `tax_data`
#'   and, by inference, what taxa apply to each item.
#' @param class_cols (`character` or `integer`) The names or indexes of columns that contain
#'   classifications if the first input is a table. If multiple columns are specified, they will be
#'   combined in the order given. Negative column indexes mean "every column besides these columns".
#' @param class_sep (`character`) One or more separators that delineate taxon names in a
#'   classification. For example, if one column had `"Homo sapiens"` and another had
#'   `"Animalia;Chordata;Mammalia;Primates;Hominidae"`, then `class_sep = c(" ", ";")`. All
#'   separators are applied to each column so order does not matter.
#' @param sep_is_regex (`TRUE`/`FALSE`) Whether or not `class_sep` should be used as a [regular
#'   expression](https://en.wikipedia.org/wiki/Regular_expression).
#' @param class_key (`character` of length 1) The identity of the capturing groups defined using
#'   `class_regex`. The length of `class_key` must be equal to the number of capturing groups
#'   specified in `class_regex`. Any names added to the terms will be used as column names in the
#'   output. At least one `"taxon_name"` must be specified. Only `"info"` can be used multiple
#'   times. Each term must be one of those described below: * `taxon_name`: The name of a taxon. Not
#'   necessarily unique, but are interpretable by a particular `database`. Requires an internet
#'   connection. * `taxon_rank`: The rank of the taxon. This will be used to add rank info into the
#'   output object that can be accessed by `out$taxon_ranks()`. * `info`: Arbitrary taxon info you
#'   want included in the output. Can be used more than once.
#' @param class_regex (`character` of length 1) A regular expression with capturing groups
#'   indicating the locations of data for each taxon in the `class` term in the `key` argument. The
#'   identity of the information must be specified using the `class_key` argument. The `class_sep`
#'   option can be used to split the classification into data for each taxon before matching. If
#'   `class_sep` is `NULL`, each match of `class_regex` defines a taxon in the classification.
#' @param class_reversed If `TRUE`, then classifications go from specific to general. For example:
#'   `Abditomys latidens : Muridae : Rodentia : Mammalia : Chordata`.
#' @param include_match (`logical` of length 1) If `TRUE`, include the part of the input matched by
#'   `class_regex` in the output object.
#' @param mappings (named `character`) This defines how the taxonomic information in `tax_data`
#'   applies to data set in `datasets`. This option should have the same number of inputs as
#'   `datasets`, with values corresponding to each data set. The names of the character vector
#'   specify what information in `tax_data` is shared with info in each `dataset`, which is
#'   specified by the corresponding values of the character vector. If there are no shared
#'   variables, you can add `NA` as a placeholder, but you could just leave that data out since it
#'   is not benefiting from being in the taxmap object. The names/values can be one of the
#'   following: * For tables, the names of columns can be used. * `"{{index}}"` : This means to use
#'   the index of rows/items * `"{{name}}"`  : This means to use row/item names. * `"{{value}}"` :
#'   This means to use the values in vectors or lists. Lists will be converted to vectors using
#'   [unlist()].
#' @param include_tax_data (`TRUE`/`FALSE`) Whether or not to include `tax_data` as a dataset, like
#'   those in `datasets`.
#' @param named_by_rank (`TRUE`/`FALSE`) If  `TRUE` and the input is a table with columns named by
#'   ranks or a list of vectors with each vector named by ranks, include that rank info in the
#'   output object, so it can be accessed by `out$taxon_ranks()`. If `TRUE`, taxa with different
#'   ranks, but the same name and location in the taxonomy, will be considered different taxa.
#'   Cannot be used with the `sep`, `class_regex`, or `class_key` options.
#'
#' @family parsers
#'
#' @examples
#'  # Read a vector of classifications
#'  my_taxa <- c("Mammalia;Carnivora;Felidae",
#'               "Mammalia;Carnivora;Felidae",
#'               "Mammalia;Carnivora;Ursidae")
#'  parse_tax_data(my_taxa, class_sep = ";")
#'
#'  # Read a list of classifications
#'  my_taxa <- list("Mammalia;Carnivora;Felidae",
#'                 "Mammalia;Carnivora;Felidae",
#'                 "Mammalia;Carnivora;Ursidae")
#'  parse_tax_data(my_taxa, class_sep = ";")
#'
#'  # Read classifications in a table in a single column
#'  species_data <- data.frame(tax = c("Mammalia;Carnivora;Felidae",
#'                                     "Mammalia;Carnivora;Felidae",
#'                                     "Mammalia;Carnivora;Ursidae"),
#'                            species_id = c("A", "B", "C"))
#'  parse_tax_data(species_data, class_sep = ";", class_cols = "tax")
#'
#'  # Read classifications in a table in multiple columns
#'  species_data <- data.frame(lineage = c("Mammalia;Carnivora;Felidae",
#'                                         "Mammalia;Carnivora;Felidae",
#'                                         "Mammalia;Carnivora;Ursidae"),
#'                             species = c("Panthera leo",
#'                                         "Panthera tigris",
#'                                         "Ursus americanus"),
#'                             species_id = c("A", "B", "C"))
#'  parse_tax_data(species_data, class_sep = c(" ", ";"),
#'                 class_cols = c("lineage", "species"))
#'
#'  # Read classification tables with one column per rank
#'  species_data <- data.frame(class = c("Mammalia", "Mammalia", "Mammalia"),
#'                             order = c("Carnivora", "Carnivora", "Carnivora"),
#'                             family = c("Felidae", "Felidae", "Ursidae"),
#'                             genus = c("Panthera", "Panthera", "Ursus"),
#'                             species = c("leo", "tigris", "americanus"),
#'                             species_id = c("A", "B", "C"))
#'   parse_tax_data(species_data, class_cols = 1:5)
#'   parse_tax_data(species_data, class_cols = 1:5,
#'                  named_by_rank = TRUE) # makes `taxon_ranks()` work
#'
#'  # Classifications with extra information
#'  my_taxa <- c("Mammalia_class_1;Carnivora_order_2;Felidae_genus_3",
#'               "Mammalia_class_1;Carnivora_order_2;Felidae_genus_3",
#'               "Mammalia_class_1;Carnivora_order_2;Ursidae_genus_3")
#'  parse_tax_data(my_taxa, class_sep = ";",
#'                 class_regex = "(.+)_(.+)_([0-9]+)",
#'                 class_key = c(my_name = "taxon_name",
#'                               a_rank = "taxon_rank",
#'                               some_num = "info"))
#'
#'
#'   # --- Parsing multiple datasets at once (advanced) ---
#'   # The rest is one example for how to classify multiple datasets at once.
#'
#'   # Make example data with taxonomic classifications
#'   species_data <- data.frame(tax = c("Mammalia;Carnivora;Felidae",
#'                                      "Mammalia;Carnivora;Felidae",
#'                                      "Mammalia;Carnivora;Ursidae"),
#'                              species = c("Panthera leo",
#'                                          "Panthera tigris",
#'                                          "Ursus americanus"),
#'                              species_id = c("A", "B", "C"))
#'
#'   # Make example data associated with the taxonomic data
#'   # Note how this does not contain classifications, but
#'   # does have a varaible in common with "species_data" ("id" = "species_id")
#'   abundance <- data.frame(id = c("A", "B", "C", "A", "B", "C"),
#'                           sample_id = c(1, 1, 1, 2, 2, 2),
#'                           counts = c(23, 4, 3, 34, 5, 13))
#'
#'   # Make another related data set named by species id
#'   common_names <- c(A = "Lion", B = "Tiger", C = "Bear", "Oh my!")
#'
#'   # Make another related data set with no names
#'   foods <- list(c("ungulates", "boar"),
#'                 c("ungulates", "boar"),
#'                 c("salmon", "fruit", "nuts"))
#'
#'   # Make a taxmap object with these three datasets
#'   x = parse_tax_data(species_data,
#'                      datasets = list(counts = abundance,
#'                                      my_names = common_names,
#'                                      foods = foods),
#'                      mappings = c("species_id" = "id",
#'                                   "species_id" = "{{name}}",
#'                                   "{{index}}" = "{{index}}"),
#'                      class_cols = c("tax", "species"),
#'                      class_sep = c(" ", ";"))
#'
#'   # Note how all the datasets have taxon ids now
#'   x$data
#'
#'   # This allows for complex mappings between variables that other functions use
#'   map_data(x, my_names, foods)
#'   map_data(x, counts, my_names)
#'
#' @export
parse_tax_data <- function(tax_data, datasets = list(), class_cols = 1,
                           class_sep = ";", sep_is_regex = FALSE,
                           class_key = "taxon_name", class_regex = "(.*)",
                           class_reversed = FALSE,
                           include_match = TRUE,
                           mappings = c(), include_tax_data = TRUE,
                           named_by_rank = FALSE) {

  # Check for nonsensical options
  if (length(datasets) != length(mappings)) {
    stop(paste0('The `mappings` option must have the same number of values (',
                length(datasets), ') as the `datasets` option.'))
  }
  if (length(mappings) > 0 && is.null(names(mappings))) {
    stop(paste0('The mapping options must be named.'))
  }
  for (i in seq_len(length(datasets))) {
    valid_mappings <- c("{{index}}", "{{name}}", "{{value}}", NA)
    if (is.data.frame(datasets[[i]])) {
      valid_mappings <- c(valid_mappings, colnames(datasets[[i]]))
    }
    if (is.data.frame(tax_data)) {
      valid_mappings <- c(valid_mappings, colnames(tax_data))
    }
    mappings_used <- c(mappings[[i]], names(mappings[[i]]))
    invalids <- mappings_used[! mappings_used %in% valid_mappings]
    if (length(invalids) > 0) {
      stop(paste0('Invalid inputs to the `mappings` found for input "', invalids[1], '". ',
                  'The names and values of `mappings` must be one of the following: ',
                  paste0(valid_mappings, collapse = ", ")))
    }
  }
  if (! "taxon_name" %in% class_key) {
    stop('At least one value in "class_key" must be "taxon_name".')
  }
  if ("taxon_rank" %in% class_key && named_by_rank) {
    stop('"named_by_rank = TRUE" does not make sense when "taxon_rank" is in "class_key".')
  }

  # Check for zero-length inputs
  if (length_of_thing(tax_data) <= 0) {
    return(taxmap(data = list(tax_data = tax_data)))
  }

  # Check that column exists
  for (a_col in class_cols) {
    check_class_col(tax_data, a_col)
  }

  # Translate negative column indexes to positive indexes
  if (is.numeric(class_cols)) {
    if (any(class_cols < 0) && any(class_cols > 0)) {
      stop("Cannot use both negative and positive column indexes at once.")
    } else if (any(class_cols < 0)) {
      class_cols <- setdiff(seq_len(ncol(tax_data)), -class_cols)
    }
  }

  # Get classificatons
  is_list_of_frames <- FALSE
  if (is.character(tax_data)) { # is a character vector
    parsed_tax <- multi_sep_split(tax_data, fixed = !sep_is_regex, split = class_sep)
  } else if (is.data.frame(tax_data)) { # is a data.frame
    parsed_tax <- lapply(seq_len(nrow(tax_data)),
                         function(i) {
                           class_source <- unlist(lapply(tax_data[i, class_cols], as.character))
                           unlist(multi_sep_split(class_source,
                                                         fixed = !sep_is_regex,
                                                         split = class_sep))
                         })
  } else if (is.list(tax_data) && is.data.frame(tax_data[[1]])) { # is a list of data.frames
    if (length(class_cols) > 1) {
      stop("When the taxonomy source is a list of data.frames, it does not make sense to specify multiple columns.")
    }
    parsed_tax <- unlist(lapply(tax_data, function(x) {
      lapply(seq_len(nrow(x)), function(i) {
        as.character(x[1:i, class_cols])
      })
    }), recursive = FALSE)
    tax_data <- do.call(rbind, tax_data)
    is_list_of_frames <- TRUE
  } else if (is.list(tax_data) && is.character(tax_data[[1]])) { # is a list of characters
    parsed_tax <- lapply(tax_data, function(x) unlist(multi_sep_split(x, split = class_sep)))
  } else {
    stop("Unknown format for first input. Cannot parse taxonomic information.")
  }

  # Remove white space
  parsed_tax <- lapply(parsed_tax, trimws)

  # Reverse order of taxa in classifications
  if (class_reversed) {
    parsed_tax <- lapply(parsed_tax, rev)
  }

  # Check for NAs in input
  na_indexes <- which(vapply(parsed_tax, function(x) any(is.na(x)), logical(1)))
  if (length(na_indexes) > 0) {
    message('The following ', length(na_indexes),' of ', length(parsed_tax),
            ' (', to_percent(length(na_indexes) / length(parsed_tax)), ')',
            ' input indexes have `NA` in their classifications:\n',
            limited_print(na_indexes, prefix = "  ", type = "silent"))
  }

  # Extract out any taxon info
  expected_col_count <- count_capture_groups(class_regex) + 1
  if (is.null(class_sep)) { # Use mutliple matches of the class regex instead of sep
    taxon_info <- lapply(parsed_tax, function(x)
      data.frame(stringr::str_match_all(x, class_regex), stringsAsFactors = FALSE))
  } else { # only use first match if sep is applied
    taxon_info <- lapply(parsed_tax, function(x) {
      validate_regex_match(x, class_regex)
      output <- data.frame(stringr::str_match(x, class_regex),
                           stringsAsFactors = FALSE)
      if (ncol(output) != expected_col_count) { # stringr::str_match found no matches (all NA)
        output <- as.data.frame(matrix(rep(NA_character_, nrow(output) * expected_col_count),
                                       nrow = nrow(output), ncol = expected_col_count),
                                stringsAsFactors = FALSE)
        colnames(output) <- paste0("X", seq_len(expected_col_count))
      }
      rownames(output) <- names(x)
      return(output)
    })
  }

  # Clean up classifications to only include taxon name and rank
  if ("taxon_rank" %in% class_key) {
    parsed_tax <- lapply(taxon_info, function(x) {
      stats::setNames(x[[which(class_key == "taxon_name") + 1]],
                      x[[which(class_key == "taxon_rank") + 1]])
    })
    named_by_rank <- TRUE
  } else if (named_by_rank)  {
    parsed_tax <- lapply(taxon_info, function(x) {
      stats::setNames(x[[which(class_key == "taxon_name") + 1]],
                      rownames(x))
    })
  } else {
    parsed_tax <- lapply(taxon_info, function(x) {
      x[[which(class_key == "taxon_name") + 1]]
    })
  }

  # combine taxon info into a single table
  taxon_info <- do.call(rbind, c(taxon_info, make.row.names = FALSE))
  if (is.null(names(class_key))) { # If all capture groups are not named
    taxon_info_colnames <- make.unique(paste0(class_key, "_match"),
                                       sep = "_")
  } else { # Some capture groups might be unnamed still
    taxon_info_colnames <- ifelse(names(class_key) == "",
                                  paste0(class_key, "_match"),
                                  names(class_key))
    taxon_info_colnames <- make.unique(taxon_info_colnames, sep = "_")
  }
  names(taxon_info) <- c("match", taxon_info_colnames)

  # Create taxmap object
  output <- taxmap(.list = parsed_tax, named_by_rank = named_by_rank)

  # Add taxon ids to extracted info and add to data
  if (ncol(taxon_info) > 2) {
    my_supertaxa <- output$supertaxa(output$input_ids,
                                     value = "taxon_ids",
                                     include_input = TRUE)
    input_index <- rep(seq_len(length(my_supertaxa)),
                                  vapply(my_supertaxa, length, numeric(1)))
    output$data$class_data <- cbind(list(taxon_id = unlist(lapply(my_supertaxa, rev)),
                                         input_index = rep(seq_len(length(my_supertaxa)),
                                                           vapply(my_supertaxa, length, numeric(1)))),
                                    taxon_info,
                                    stringsAsFactors = FALSE)
    if (include_match) {
      match_col_index <- which(colnames(output$data$class_data) == "match")
      output$data$class_data <- cbind(output$data$class_data[-match_col_index],
                                      list(regex_match = output$data$class_data$match),
                                      stringsAsFactors = FALSE)
    } else {
      output$data$class_data$match  <- NULL
    }
    output$data$class_data <- tibble::as_tibble(output$data$class_data)
  }

  # Add taxonomic source to datasets
  if (include_tax_data) {
    datasets <- c(list(tax_data = tax_data), datasets)
    mappings <- c("{{index}}" = "{{index}}", mappings)
  }


  # Add additional data sets
  name_datset <- function(dataset, sort_var) {
    if (!is.na(sort_var)) {
      target_value <- get_sort_var(dataset, sort_var)
      source_value <- get_sort_var(tax_data, names(sort_var))
      obs_taxon_ids <- unname(output$input_ids[match(target_value, source_value)])
      if (is.data.frame(dataset)) {
        dataset <- dplyr::bind_cols(dplyr::tibble(taxon_id = obs_taxon_ids),
                                    dataset)
      } else {
        names(dataset) <- obs_taxon_ids
      }
    }
    return(dataset)
  }

  output$data <- c(output$data,
                   stats::setNames(lapply(seq_len(length(datasets)),
                                          function(i) name_datset(datasets[[i]],
                                                                  mappings[i])),
                                   names(datasets)))

  # Convert additional tables to tibbles
  are_tables <- vapply(output$data, is.data.frame, logical(1))
  output$data[are_tables] <- lapply(output$data[are_tables], tibble::as_tibble)

  # Fix incorrect taxon ids in data if a list of data.frames is given
  if (is_list_of_frames && include_tax_data) {
    output$data$tax_data <- output$data$tax_data[!duplicated(output$data[[1]]), ]
  }

  # Check format of data sets
  check_taxmap_data(output)

  # Put tax_data first
  if ("tax_data" %in% names(output$data)) {
    tax_data_index <- which(names(output$data) == "tax_data")
    output$data <- c(output$data[tax_data_index], output$data[-tax_data_index])
  }

  return(output)
}



#' Convert one or more data sets to taxmap
#'
#' Looks up taxonomic data from NCBI sequence IDs, taxon IDs, or taxon names
#' that are present in a table, list, or vector. Also can incorporate additional
#' associated datasets.
#'
#' @param tax_data A table, list, or vector that contain sequence IDs, taxon
#'   IDs, or taxon names.
#'   * tables: The `column` option must be used to specify which column
#'   contains the sequence IDs, taxon IDs, or taxon names.
#'   * lists: There must be only one item per list entry unless the `column`
#'   option is used to specify what item to use in each list entry.
#'   * vectors: simply a vector of sequence IDs, taxon IDs, or taxon names.
#' @param type What type of information can be used to look up the
#'   classifications. Takes one of the following values:
#'   * `"seq_id"`: A database sequence ID with an associated classification
#'   (e.g. NCBI accession numbers).
#'   * `"taxon_id"`: A reference database taxon ID (e.g. a NCBI taxon ID)
#'   * `"taxon_name"`: A single taxon name (e.g. "Homo sapiens" or "Primates")
#'   * `"fuzzy_name"`: A single taxon name, but check for misspellings first.
#'   Only use if you think there are misspellings. Using `"taxon_name"` is
#'   faster.
#' @param column (`character` or `integer`) The name or index of the column that
#'   contains information used to lookup classifications. This only applies when
#'   a table or list is supplied to `tax_data`.
#' @param datasets Additional lists/vectors/tables that should be included in
#'   the resulting `taxmap` object. The `mappings` option is use to specify how
#'   these data sets relate to the `tax_data` and, by inference, what taxa apply
#'   to each item.
#' @param mappings (named `character`) This defines how the taxonomic
#'   information in `tax_data` applies to data in `datasets`. This option
#'   should have the same number of inputs as `datasets`, with values
#'   corresponding to each dataset. The names of the character vector specify
#'   what information in `tax_data` is shared with info in each `dataset`, which
#'   is specified by the corresponding values of the character vector. If there
#'   are no shared variables, you can add `NA` as a placeholder, but you could
#'   just leave that data out since it is not benefiting from being in the
#'   taxmap object. The names/values can be one of the following:
#'   * For tables, the names of columns can be used.
#'   * `"{{index}}"` : This means to use the index of rows/items
#'   * `"{{name}}"`  : This means to use row/item names.
#'   * `"{{value}}"` : This means to use the values in vectors or lists. Lists
#'   will be converted to vectors using [unlist()].
#' @param database (`character`) The name of a database to use to look up
#'   classifications. Options include "ncbi", "itis", "eol", "col", "tropicos",
#'  and "nbn".
#' @param include_tax_data (`TRUE`/`FALSE`) Whether or not to include `tax_data`
#'   as a dataset, like those in `datasets`.
#' @param use_database_ids (`TRUE`/`FALSE`) Whether or not to use downloaded
#'   database taxon ids instead of arbitrary, automatically-generated taxon ids.
#' @param ask  (`TRUE`/`FALSE`) Whether or not to prompt the user for input.
#'   Currently, this would only happen when looking up the taxonomy of a taxon
#'   name with multiple matches. If `FALSE`, taxa with multiple hits are treated
#'   as if they do not exist in the database. This might change in the future if
#'   we can find an elegant way of handling this.
#'
#' @section Failed Downloads: If you have invalid inputs or a download fails for
#'   another reason, then there will be a "unknown" taxon ID as a placeholder
#'   and failed inputs will be assigned to this ID. You can remove these using
#'   [filter_taxa()] like so: `filter_taxa(result, taxon_ids != "unknown")`. Add
#'   `drop_obs = FALSE` if you want the input data, but want to remove the
#'   taxon.
#'
#' @family parsers
#'
#' @examples \dontrun{
#'
#'   # Look up taxon names in vector from NCBI
#'   lookup_tax_data(c("homo sapiens", "felis catus", "Solanaceae"),
#'                   type = "taxon_name")
#'
#'   # Look up taxon names in list from NCBI
#'   lookup_tax_data(list("homo sapiens", "felis catus", "Solanaceae"),
#'                   type = "taxon_name")
#'
#'   # Look up taxon names in table from NCBI
#'   my_table <- data.frame(name = c("homo sapiens", "felis catus"),
#'                          decency = c("meh", "good"))
#'   lookup_tax_data(my_table, type = "taxon_name", column = "name")
#'
#'   # Look up taxon names from NCBI with fuzzy matching
#'   lookup_tax_data(c("homo sapienss", "feles catus", "Solanacese"),
#'                   type = "fuzzy_name")
#'
#'   # Look up taxon names from a different database
#'   lookup_tax_data(c("homo sapiens", "felis catus", "Solanaceae"),
#'                   type = "taxon_name", database = "ITIS")
#'
#'   # Prevent asking questions for ambiguous taxon names
#'   lookup_tax_data(c("homo sapiens", "felis catus", "Solanaceae"),
#'                   type = "taxon_name", database = "ITIS", ask = FALSE)
#'
#'   # Look up taxon IDs from NCBI
#'   lookup_tax_data(c("9689", "9694", "9643"), type = "taxon_id")
#'
#'   # Look up sequence IDs from NCBI
#'   lookup_tax_data(c("AB548412", "FJ358423", "DQ334818"),
#'                   type = "seq_id")
#'
#'   # Make up new taxon IDs instead of using the downloaded ones
#'   lookup_tax_data(c("AB548412", "FJ358423", "DQ334818"),
#'                   type = "seq_id", use_database_ids = FALSE)
#'
#'
#'   # --- Parsing multiple datasets at once (advanced) ---
#'   # The rest is one example for how to classify multiple datasets at once.
#'
#'   # Make example data with taxonomic classifications
#'   species_data <- data.frame(tax = c("Mammalia;Carnivora;Felidae",
#'                                      "Mammalia;Carnivora;Felidae",
#'                                      "Mammalia;Carnivora;Ursidae"),
#'                              species = c("Panthera leo",
#'                                          "Panthera tigris",
#'                                          "Ursus americanus"),
#'                              species_id = c("A", "B", "C"))
#'
#'   # Make example data associated with the taxonomic data
#'   # Note how this does not contain classifications, but
#'   # does have a varaible in common with "species_data" ("id" = "species_id")
#'   abundance <- data.frame(id = c("A", "B", "C", "A", "B", "C"),
#'                           sample_id = c(1, 1, 1, 2, 2, 2),
#'                           counts = c(23, 4, 3, 34, 5, 13))
#'
#'   # Make another related data set named by species id
#'   common_names <- c(A = "Lion", B = "Tiger", C = "Bear", "Oh my!")
#'
#'   # Make another related data set with no names
#'   foods <- list(c("ungulates", "boar"),
#'                 c("ungulates", "boar"),
#'                 c("salmon", "fruit", "nuts"))
#'
#'   # Make a taxmap object with these three datasets
#'   x = lookup_tax_data(species_data,
#'                       type = "taxon_name",
#'                       datasets = list(counts = abundance,
#'                                       my_names = common_names,
#'                                       foods = foods),
#'                       mappings = c("species_id" = "id",
#'                                    "species_id" = "{{name}}",
#'                                    "{{index}}" = "{{index}}"),
#'                       column = "species")
#'
#'   # Note how all the datasets have taxon ids now
#'   x$data
#'
#'   # This allows for complex mappings between variables that other functions use
#'   map_data(x, my_names, foods)
#'   map_data(x, counts, my_names)
#' }
#' @export
lookup_tax_data <- function(tax_data, type, column = 1, datasets = list(),
                            mappings = c(), database = "ncbi",
                            include_tax_data = TRUE, use_database_ids = TRUE,
                            ask = TRUE) {

  # Check for zero-length inputs
  if (length_of_thing(tax_data) <= 0) {
    return(taxmap(data = list(tax_data = tax_data)))
  }

  # Make sure taxize is installed
  check_for_pkg("taxize")

  # Check that a supported database is being used
  supported_databases <- names(database_list)
  database <- tolower(database)
  if (! database %in% supported_databases) {
    stop(paste0('The database "', database,
                '" is not a valid database for looking up that taxonomy of ',
                'sequnece ids. Valid choices include:\n',
                limited_print(supported_databases, type = "silent")))
  }

  # Check that column exists
  check_class_col(tax_data, column)

  # Hidden parameters
  batch_size <- 100
  max_print <- 10
  internal_class_sep <- "||||"
  internal_class_name <- "___class_string___"

  # Define lookup functions
  format_class_table <- function(class_table) {
    # Complain about failed queries
    failed_queries <- is.na(class_table)
    if (sum(failed_queries) > 0) {
      failed_names <- unique(names(failed_queries[failed_queries]))
      error_msg <- paste0('The following ', length(failed_names),
                          ' unique inputs could not be looked up:\n  ',
                          limited_print(failed_names, type = "silent"))
      if (ask) {
        error_msg <- paste0(error_msg, "\n",
                            'This probably means they dont exist in the database "', database, '".')
      } else {
        error_msg <- paste0(error_msg, "\n",
                            'This probably means they dont exist in the database "', database,
                            '" or have multiple matches. ',
                            'Use "ask = TRUE" to specify which is the correct match when multiple matches occur.')
      }
      warning(error_msg, call. = FALSE)
    }

    # Rename columns of result
    col_names <- c(paste0(database, "_name"),
                   paste0(database, "_rank"),
                   paste0(database, "_id"))
    class_table[! failed_queries] <- lapply(class_table[! failed_queries],
                                       function(x) stats::setNames(x, col_names))

    # Add placeholders for failed requests
    class_table[failed_queries] <- lapply(seq_len(sum(failed_queries)),
                                     function(i) {
                                       out <- dplyr::tibble(a = "unknown taxon",
                                                            b = "unknown rank",
                                                            c = "unknown")
                                       colnames(out) <- col_names
                                       return(out)
                                     })
    return(class_table)
  }


  use_taxon_id <- function(ids) {
    message("Looking up classifications for ", length(unique(ids)),
            ' unique taxon IDs from database "', database, '"...')

    # Look up classifications
    lookup_all <- function(ids) {
      lookup_one <- function(id) {
        taxize::classification(id, ask = FALSE, rows = 1, db = database,
                               message = FALSE)
      }
      output <- progress_lapply(ids, lookup_one)
      return(output)
    }
    tryCatch(msgs <- utils::capture.output(raw_result <- map_unique(ids, lookup_all),
                                           type = "message"),
             error = function(e) stop(e))


    # Remove repeated messages (e.g. no NCBI API key)
    on.exit(message(paste0(unique(msgs), collapse = "\n")))

    # Reformat result
    result <- stats::setNames(unlist(raw_result, recursive = FALSE), ids)
    format_class_table(result)
  }

  use_seq_id <- function(ids) {
    # Check that a supported database is being used
    supported_databases <- c("ncbi")
    if (! database %in% supported_databases) {
      stop(call. = FALSE,
           paste0('The database "', database,
                  '" is not a valid database for looking up that taxonomy of ',
                  'sequnece ids. Valid choices include:\n',
                  limited_print(supported_databases, type = "silent")))
    }

    # Look up classifications
    message("Looking up classifications for ", length(unique(ids)), " unique sequence IDs from NCBI...")
    lookup_all <- function(ids) {
      lookup_one <- function(id) {
        taxize::classification(taxize::genbank2uid(id)[1], db = database)
      }
      output <- progress_lapply(ids, lookup_one)
      return(output)
    }
    tryCatch(msgs <- utils::capture.output(raw_result <- map_unique(ids, lookup_all),
                                  type = "message"),
             error = function(e) stop(e))

    # Remove repeated messages (e.g. no NCBI API key)
    on.exit(message(paste0(unique(msgs), collapse = "\n")))

    # Reformat result
    result <- stats::setNames(unlist(raw_result, recursive = FALSE), ids)
    format_class_table(result)
  }

  use_taxon_name <- function(my_names) {
    message("Looking up classifications for ", length(unique(my_names)),
            ' unique taxon names from database "', database, '"...')

    # Look up classifications
    lookup_all <- function(my_names) {
      lookup_one <- function(my_name) {
        taxize::classification(my_name, ask = ask, db = database,
                               messages = FALSE)
       }
      output <- progress_lapply(my_names, lookup_one)
      return(output)
    }
    tryCatch(msgs <- utils::capture.output(raw_result <- map_unique(my_names, lookup_all),
                                  type = "message"),
             error = function(e) stop(e))

    # Remove repeated messages (e.g. no NCBI API key)
    on.exit(message(paste0(unique(msgs), collapse = "\n")))

    # Reformat result
    result <- stats::setNames(unlist(raw_result, recursive = FALSE), my_names)
    format_class_table(result)
  }

  use_taxon_name_fuzzy <- function(my_names) {
    message("Looking up classifications for ", length(unique(my_names)),
            ' unique taxon names from database "', database, '" using fuzzy name matching...')

    # Look up similar taxon names
    corrected <- map_unique(my_names, correct_taxon_names, database = database)

    # Check for not found names
    not_found <- unique(names(corrected[is.na(corrected)]))
    if (length(not_found) > 0) {
      warning(call. = FALSE,
              "No taxon name was found that was similar to the following ",
              length(not_found), " inputs:\n  ",
              limited_print(type = "silent", not_found))

    }

    # Replace not found values with original input
    corrected[is.na(corrected)] <- names(corrected[is.na(corrected)])

    # Check for changed names
    changed <- tolower(names(corrected)) != tolower(corrected) & ! is.na(corrected) & ! is.na(names(corrected))
    if (any(changed)) {
      before <- names(corrected[changed])
      after <- corrected[changed]
      change_char <- unique(paste0('"', before, '" -> "', after, '"'))
      message("The following ", length(change_char), " names were corrected before looking up classifications:\n  ",
              limited_print(type = "silent", change_char))
    }

    # Run standard taxon name lookup
    use_taxon_name(corrected)
  }

  lookup_funcs <- list("seq_id" = use_seq_id,
                       "taxon_id" = use_taxon_id,
                       "taxon_name" = use_taxon_name,
                       "fuzzy_name" = use_taxon_name_fuzzy)

  # Get query information
  if (is.data.frame(tax_data)) { # is table
    query <- as.character(tax_data[[column]])
  } else if (is.list(tax_data)) { # is list
    query <- vapply(tax_data,
                    function(x) as.character(x[[column]]),
                    character(1))
  } else if (is.vector(tax_data)) { # is vector
    query <- as.character(tax_data)
  }

  # Look up taxonomic classifications
  if (! type %in% names(lookup_funcs)) {
    stop(paste0('Invalid "type" option. It must be one of the following:\n  ',
                paste0(names(lookup_funcs), collapse = ", ")))
  }
  classifications <- lookup_funcs[[type]](query)
  class_strings <- unlist(lapply(classifications, function(x) {
    lapply(seq_len(nrow(x)), function(i) {
      paste0(as.character(unlist(x[1:i, 1])), "___", as.character(unlist(x[1:i, 2])), collapse = internal_class_sep)
    })
  }))
  combined_class <- do.call(rbind, unname(classifications))
  internal_class_frame <- stats::setNames(data.frame(class_strings,
                                                     stringsAsFactors = FALSE),
                                          internal_class_name)
  combined_class <- cbind(internal_class_frame,
                          combined_class,
                          stringsAsFactors = FALSE)

  # Add mapping columns to classfication data
  tax_data_indexes <- cumsum(vapply(classifications, nrow, numeric(1)))
  mappping_cols <- unique(c(names(mappings), "{{index}}", "{{name}}"))
  if (!is.data.frame(tax_data)) {
    mappping_cols <- c(mappping_cols, "{{value}}")
  }
  for (col in mappping_cols) {
    combined_class[tax_data_indexes, col] <- get_sort_var(tax_data, col)
  }

  # Add input data to datasets included in the resulting object
  if (include_tax_data) {
    datasets <- c(list(query_data = tax_data), datasets)
    mappings <- c("{{index}}" = "{{index}}", mappings)
  }

  # Make taxmap object
  output <- parse_tax_data(tax_data = combined_class,
                           datasets = datasets,
                           class_cols = 1,
                           class_sep = internal_class_sep,
                           class_key = c('taxon_name', 'taxon_rank'),
                           class_regex = '^(.+)___(.+)$',
                           mappings = mappings,
                           include_tax_data = include_tax_data)
  output$data$class_data <- NULL

  if (include_tax_data) {
    # Remove mapping columns from output
    output$data$tax_data[mappping_cols] <- NULL

    # Remove class column from output
    output$data$tax_data[internal_class_name] <- NULL

    # Fix incorrect taxon ids for tax_data
    #   This is due to the "{{index}}" being interpreted as a column name,
    #   which is needed for the user-defined data sets to be parsed right.
    output$data$tax_data$taxon_id <- output$input_ids

    # Remove duplicates from `tax_data`
    output$data$tax_data <- output$data$tax_data[!duplicated(output$data$tax_data), ]
  }

  # Replace standard taxon ids with database taxon ids
  if (use_database_ids) {
    taxon_id_col <- paste0(database, "_id")
    # I am not sure why the following line works...
    new_ids <- unique(combined_class[[taxon_id_col]])[match(output$taxon_ids(),
                                                            unique(output$input_ids))]
    output$replace_taxon_ids(new_ids)
  }

  return(output)
}


#' Get a vector from a vector/list/table to be used in mapping
#'
#' @param data A vector/list/table
#' @param var What to get.
#'   * For tables, the names of columns can be used.
#'   * `"{{index}}"` : This means to use the index of rows/items
#'   * `"{{name}}"`  : This means to use row/item names.
#'   * `"{{value}}"` : This means to use the values in vectors or lists. Lists
#'
#' @keywords internal
get_sort_var <- function(data, var) {
  if (is.data.frame(data) && var %in% colnames(data)) {
    return(data[[var]])
  } else if (var == "{{index}}") {
    if (is.data.frame(data)) {
      return(seq_len(nrow(data)))
    } else {
      return(seq_len(length(data)))
    }
  } else if (var == "{{name}}") {
    if (is.data.frame(data)) {
      return(rownames(data))
    } else {
      if (is.null(names(data))) {
        return(rep(NA_character_, length(data)))
      } else {
        return(names(data))
      }
    }
  }  else if (var == "{{value}}") {
    if (is.data.frame(data)) {
      stop("The `{{value}}` setting of the `mappings` option cannot be used with data.frames.")
    } else {
      return(unlist(data))
    }
  }  else {
    stop(paste0('No column named "', var, '"."'))
  }
}



#' Extracts taxonomy info from vectors with regex
#'
#' Convert taxonomic information in a character vector into a [taxmap()] object.
#' The location and identity of important information in the input is specified
#' using a [regular expression](https://en.wikipedia.org/wiki/Regular_expression)
#' with capture groups and a corresponding key. An object of type [taxmap()] is
#' returned containing the specified information. See the `key` option for
#' accepted sources of taxonomic information.
#'
#' @param tax_data A vector from which to extract taxonomy information.
#' @param key (`character`) The identity of the capturing groups defined using
#'   `regex`. The length of `key` must be equal to the number of capturing
#'   groups specified in `regex`. Any names added to the terms will be used as
#'   column names in the output. Only `"info"` can be used multiple times. Each
#'   term must be one of those described below:
#'   * `taxon_id`: A unique numeric id for a taxon for a particular `database`
#'   (e.g. ncbi accession number). Requires an internet connection.
#'   * `taxon_name`: The name of a taxon (e.g. "Mammalia" or "Homo sapiens").
#'   Not necessarily unique, but interpretable by a particular `database`.
#'   Requires an internet connection.
#'   * `fuzzy_name`: The name of a taxon, but check for misspellings first.
#'   Only use if you think there are misspellings. Using `"taxon_name"` is
#'   faster.
#'   * `class`: A list of taxon information that constitutes the full taxonomic
#'   classification (e.g. "K_Mammalia;P_Carnivora;C_Felidae"). Individual
#'   taxa are separated by the `class_sep` argument and the information is
#'   parsed by the `class_regex` and `class_key` arguments.
#'   * `seq_id`: Sequence ID for a particular database that is associated with a
#'   taxonomic classification. Currently only works with the "ncbi" database.
#'   * `info`: Arbitrary taxon info you want included in the output. Can be used
#'   more than once.
#' @param regex (`character` of length 1) A regular expression with capturing
#'   groups indicating the locations of relevant information. The identity of
#'   the information must be specified using the `key` argument.
#' @param class_key (`character` of length 1) The identity of the capturing
#'   groups defined using `class_regex`. The length of `class_key` must be equal
#'   to the number of capturing groups specified in `class_regex`. Any names
#'   added to the terms will be used as column names in the output. Only
#'   `"info"` can be used multiple times. Each term must be one of those
#'   described below:
#'   * `taxon_name`: The name of a taxon. Not necessarily unique.
#'   * `taxon_rank`: The rank of the taxon. This will be used to add rank info
#'   into the output object that can be accessed by `out$taxon_ranks()`.
#'   * `info`: Arbitrary taxon info you want included in the output. Can be used
#'   more than once.
#' @param class_regex (`character` of length 1)
#'   A regular expression with capturing groups indicating the locations of data
#'   for each taxon in the `class` term in the `key` argument. The identity of
#'   the information must be specified using the `class_key` argument. The
#'   `class_sep` option can be used to split the classification into data for
#'   each taxon before matching. If `class_sep` is `NULL`, each match of
#'   `class_regex` defines a taxon in the classification.
#' @param class_sep (`character` of length 1)
#'   Used with the `class` term in the `key` argument. The character(s) used to
#'   separate individual taxa within a classification. After the string defined
#'   by the `class` capture group in `regex` is split by `class_sep`, its
#'   capture groups are extracted by `class_regex` and defined by `class_key`.
#'   If `NULL`, every match of `class_regex` is used instead with first
#'   splitting by `class_sep`.
#' @param sep_is_regex (`TRUE`/`FALSE`) Whether or not `class_sep` should be
#'   used as a [regular expression](https://en.wikipedia.org/wiki/Regular_expression).
##' @param class_rev (`logical` of length 1)
#'   Used with the `class` term in the `key` argument. If `TRUE`, the order of
#'   taxon data in a classification is reversed to be specific to broad.
#' @param database (`character` of length 1) The name of the database that
#'   patterns given in `parser` will apply to. Valid databases include "ncbi",
#'   "itis", "eol", "col", "tropicos", "nbn", and "none". `"none"` will cause no
#'   database to be queried; use this if you want to not use the internet. NOTE:
#'   Only `"ncbi"` has been tested extensively so far.
#' @param include_match (`logical` of length 1) If `TRUE`, include the part of
#'   the input matched by `regex` in the output object.
#' @param include_tax_data (`TRUE`/`FALSE`) Whether or not to include `tax_data`
#'   as a dataset.
#'
#' @section Failed Downloads: If you have invalid inputs or a download fails for
#'   another reason, then there will be a "unknown" taxon ID as a placeholder
#'   and failed inputs will be assigned to this ID. You can remove these using
#'   [filter_taxa()] like so: `filter_taxa(result, taxon_ids != "unknown")`. Add
#'   `drop_obs = FALSE` if you want the input data, but want to remove the
#'   taxon.
#'
#' @family parsers
#'
#' @return Returns an object of type [taxmap()]
#'
#' @examples
#'
#' \dontrun{
#'
#'   # For demonstration purposes, the following example dataset has all the
#'   # types of data that can be used, but any one of them alone would work.
#'   raw_data <- c(
#'   ">id:AB548412-tid:9689-Panthera leo-tax:K_Mammalia;P_Carnivora;C_Felidae;G_Panthera;S_leo",
#'   ">id:FJ358423-tid:9694-Panthera tigris-tax:K_Mammalia;P_Carnivora;C_Felidae;G_Panthera;S_tigris",
#'   ">id:DQ334818-tid:9643-Ursus americanus-tax:K_Mammalia;P_Carnivora;C_Felidae;G_Ursus;S_americanus"
#'   )
#'
#'   # Build a taxmap object from classifications
#'   extract_tax_data(raw_data,
#'                    key = c(my_seq = "info", my_tid = "info", org = "info", tax = "class"),
#'                    regex = "^>id:(.+)-tid:(.+)-(.+)-tax:(.+)$",
#'                    class_sep = ";", class_regex = "^(.+)_(.+)$",
#'                    class_key = c(my_rank = "info", tax_name = "taxon_name"))
#'
#'   # Build a taxmap object from taxon ids
#'   # Note: this requires an internet connection
#'   extract_tax_data(raw_data,
#'                    key = c(my_seq = "info", my_tid = "taxon_id", org = "info", tax = "info"),
#'                    regex = "^>id:(.+)-tid:(.+)-(.+)-tax:(.+)$")
#'
#'   # Build a taxmap object from ncbi sequence accession numbers
#'   # Note: this requires an internet connection
#'   extract_tax_data(raw_data,
#'                    key = c(my_seq = "seq_id", my_tid = "info", org = "info", tax = "info"),
#'                    regex = "^>id:(.+)-tid:(.+)-(.+)-tax:(.+)$")
#'
#'   # Build a taxmap object from taxon names
#'   # Note: this requires an internet connection
#'   extract_tax_data(raw_data,
#'                    key = c(my_seq = "info", my_tid = "info", org = "taxon_name", tax = "info"),
#'                    regex = "^>id:(.+)-tid:(.+)-(.+)-tax:(.+)$")
#' }
#' @export
extract_tax_data <- function(tax_data, key, regex, class_key = "taxon_name",
                             class_regex = "(.*)", class_sep = NULL,
                             sep_is_regex = FALSE,
                             class_rev = FALSE, database = "ncbi",
                             include_match = FALSE, include_tax_data = TRUE) {

  # Check for zero-length inputs
  if (length_of_thing(tax_data) <= 0) {
    return(taxmap(data = list(tax_data = tax_data)))
  }

  # Check regex/keys
  key <- validate_regex_key_pair(regex, key, multiple_allowed = "info")
  class_key <- validate_regex_key_pair(class_regex, class_key, multiple_allowed = "info")

  # classification sep
  if (!is.null(class_sep) && (class(class_sep) != "character" | length(class_sep) != 1)) {
    stop('"class_sep" must be a character vector of length 1 or NULL')
  }

  # Boolean options
  if (class(class_rev) != "logical" | length(class_rev) != 1) {
    stop('"class_rev" must be TRUE/FALSE')
  }
  if (class(include_match) != "logical" | length(include_match) != 1) {
    stop('"include_match" must be TRUE/FALSE')
  }
  if (class(include_tax_data) != "logical" | length(include_tax_data) != 1) {
    stop('"include_tax_data" must be TRUE/FALSE')
  }

  # database
  valid_databases <- c(names(database_list), "none")
  if (! database %in% valid_databases) {
    stop(paste0('Unknown database "', database,
                '". Accepted database names include:\n    "',
                paste0(valid_databases, collapse = ", ")))
  }

  # Extract capture groups
  parsed_input <- data.frame(stringr::str_match(tax_data, regex), stringsAsFactors = FALSE)
  colnames(parsed_input) <- c("input", names(key))
  parsed_input <- cbind(parsed_input[-1], list(input = parsed_input$input),
                        stringsAsFactors = FALSE)

  # Complain about failed matches
  failed <- which(apply(is.na(parsed_input), MARGIN = 1, FUN = all))
  if (length(failed) > 0) {
    warning(paste0("The following ", length(failed), " input indexes failed to match the regex supplied:\n",
                   limited_print(failed, type = "silent")), call. = FALSE)
    parsed_input <- parsed_input[-failed, ]
  }

  # Use parse_tax_data if the input is a classification
  if ("class" %in% key) {
    output <- parse_tax_data(tax_data = parsed_input,
                             class_cols = which(key == "class"),
                             class_sep = class_sep,
                             class_key = class_key,
                             sep_is_regex = sep_is_regex,
                             class_regex = class_regex,
                             include_match = include_match,
                             include_tax_data = include_tax_data)
    if (!include_match) {
      output$data$tax_data[which(key == "class") + 1] <- NULL
    }
  }


  # Use lookup_tax_data for taxon names, ids, and sequence ids
  if (any(key %in% c("taxon_name", "taxon_id", "seq_id", "fuzzy_name"))) {
    my_type <- key[key != "info"][1]
    output <- lookup_tax_data(tax_data = parsed_input, type = my_type,
                              column = names(my_type),
                              database = database,
                              include_tax_data = include_tax_data)

    if (!include_match) {
      output$data$query_data$match <- NULL
    }
  }

  return(output)
}



#' Check that all match input
#'
#' Ensure that all of a character vector matches a regex.
#' Inputs that do not match are excluded.
#'
#' @param input (\code{character})
#' @param regex (\code{character} of length 1)
#'
#' @return \code{character} Parts of \code{input} matching \code{regex}
#'
#' @keywords internal
validate_regex_match <- function(input, regex) {
  # check which input values match
  input <- as.character(input)
  not_matching <- (! grepl(pattern = regex, x = input)) & (! is.na(input))

  # complain about those that dont
  if (sum(not_matching) > 0) {
    stop(call. = FALSE,
         paste0(collapse = "",
                c("The following ", sum(not_matching), " of ", length(input),
                  " input(s) could not be matched by the regex supplied:\n",
                  limited_print(input[not_matching], prefix = "  ", type = "silent"))))
  }
}


#' Check a regex-key pair
#'
#' Checks that the number of capture groups in the regex matches the length of the key.
#' Checks that only certain values of \code{key} can appear more that once.
#' Adds names to keys that will be used for column names in the output of \code{extract_taxonomy}.
#' Uses non-standard evaluation to get the name of input variables.
#'
#' @param regex (\code{character})
#' A regex with capture groups
#' @param key (\code{character})
#' A key corresponding to \code{regex}
#' @param multiple_allowed (\code{character})
#' Values of \code{key_options} that can appear more than once.
#'
#' @return Returns the result of \code{\link{match.arg}} on the key.
#'
#' @keywords internal
#' @rdname validate_regex_key_pair
validate_regex_key_pair <- function(regex, key, multiple_allowed) {

  # Non-standard evaluation
  regex_var_name <- deparse(substitute(regex))
  key_var_name <- deparse(substitute(key))

  # Check that the keys used are valid
  allowed <- c("taxon_id", "taxon_name", "info", "class", "seq_id", "fuzzy_name", "taxon_rank")
  invalid_keys <- key[! key %in% allowed]
  if (length(invalid_keys) > 0) {
    stop(paste0('Invalid key value "', invalid_keys[1], '" given.\n',
                'Valid options are: ', paste0(collapse = ", ", allowed)))
  }

  # Check key length
  regex_capture_group_count <- count_capture_groups(regex)
  key_length <- length(key)
  if (key_length != regex_capture_group_count) {
    stop(paste0(collapse = "",
                'The number of capture groups in "', regex_var_name, '" and the length of "',
                key_var_name, '" do not match.\n',
                'The key has ', key_length, ' term(s) and the regex has ', regex_capture_group_count))
  }

  # Check that only keys in `multiple_allowed` appear more than once
  counts <- table(key)
  duplicated_keys <- names(counts[counts > 1])
  invalid_duplicates <- duplicated_keys[! duplicated_keys %in% multiple_allowed]
  if (length(invalid_duplicates) > 0) {
    stop(paste0(collapse = "",
                'The following values in `', key_var_name, '` have been used more than once: ',
                paste0(collapse =", ", invalid_duplicates), '\n',
                '  Only the following keys can be duplicated: ',
                paste0(collapse =", ", multiple_allowed)))
  }

  # Apply key name defaults
  key_names <- names(key)
  if (is.null(key_names)) { key_names <- rep("", length(key)) }
  key_names[key_names == ""] <- key[key_names == ""]

  names(key) <- key_names
  return(key)
}

#' Count capture groups
#'
#' Count the number of capture groups in a regular expression.
#'
#' @param regex (\code{character} of length 1)
#'
#' @return \code{numeric} of length 1
#'
#' @source http://stackoverflow.com/questions/16046620/regex-to-count-the-number-of-capturing-groups-in-a-regex
#'
#' @keywords internal
count_capture_groups <- function(regex) {
  new_regex <- paste0(regex, "|") # Allow the regex to match nothing
  ncol(stringr::str_match(string = "", pattern = new_regex)) - 1
}



#' Check for name/index in input data
#'
#' Used by parse_tax_data and lookup_tax_data to check that columm/class_col is valid for the input data
#'
#' @param tax_data A table, list, or vector that contain sequence IDs, taxon
#'   IDs, or taxon names.
#'   * tables: The `column` option must be used to specify which column
#'   contains the sequence IDs, taxon IDs, or taxon names.
#'   * lists: There must be only one item per list entry unless the `column`
#'   option is used to specify what item to use in each list entry.
#'   * vectors: simply a vector of sequence IDs, taxon IDs, or taxon names.
#' @param column (`character` or `integer`) The name or index of the column that
#'   contains information used to lookup classifications. This only applies when
#'   a table or list is supplied to `tax_data`.
#'
#' @keywords internal
check_class_col <- function(tax_data, column) {
  if (is.data.frame(tax_data)) {
    if (is.numeric(column)) {
      if (column == 0 || abs(column) > ncol(tax_data)) {
        stop(call. = FALSE,
             'Column index "', column, '" out of bounds. Must be between 1 and ',
             ncol(tax_data), '.')
      }
    } else if (! column %in% colnames(tax_data)) {
      stop(call. = FALSE,
           'No column "', column, '" in input table. Valid columns include:\n  ',
           limited_print(colnames(tax_data), type = "silent"))
    }
  } else if (is.list(tax_data) || is.vector(tax_data)) {
    my_lengths <- vapply(tax_data, length, numeric(1))
    had_col_name <- vapply(tax_data, function(x) column %in% names(x), logical(1))
    if (is.numeric(column)) {
      if (column < 1 || any(column > my_lengths)) {
        stop(call. = FALSE,
             'Column index "', column, '" out of bounds for inputs:\n',
             limited_print(which(column > my_lengths), type = "silent"))
      }
    } else if (! all(had_col_name)) {
      stop(call. = FALSE,
           'No item named "', column, '" in the following inputs:\n',
           limited_print(which(! had_col_name), type = "silent"))
    }
  } else {
    stop(call. = FALSE,
         'Cannot read input of class "', class(tax_data)[1], '". Input must be a table, list or vector.')
  }
}


#' Look up official names from potentially misspelled names
#'
#' Look up official names from potentially misspelled names using Global Names
#' Resolver (GNR). If a result from the chosen database is present, then it is
#' used, otherwise the NCBI result is used and if that does not exist, then the
#' first result is used. Names with no match will return NA.
#'
#' @param names Potentially misspelled taxon names
#' @param database The database the names are being looked up for. If `NULL`, do
#'   not consider database.
#'
#' @return vector of names
#'
#' @keywords internal
correct_taxon_names <- function(names, database = "ncbi") {
  # Look up what the database is called in GNR
  database_key <-  c(itis = "ITIS",
                     ncbi = "NCBI",
                     eol = "EOL",
                     col = "Catalogue of Life",
                     tropicos = "Tropicos - Missouri Botanical Garden",
                     gbif = "GBIF Backbone Taxonomy",
                     wiki = "Wikispecies")
  if (! is.null(database)) {
    if (! tolower(database) %in% names(database_key)) {
      warning(call. = FALSE,
              'Can not check taxon names for database "', database,
              '". Using NCBI instead.')
      database = "ncbi"
    }
    gnr_database <- database_key[tolower(database)]
  }

  # Query GRN
  result <- taxize::gnr_resolve(names)

  # Pick results
  pick_one <- function(n) {
    one_tax_data <- result[result$user_supplied_name == n, ]
    if (nrow(one_tax_data) == 0) {
      return(NA_character_)
    } else if (nrow(one_tax_data) == 1) {
      return(one_tax_data$matched_name)
    } else if (is.null(database) || ! gnr_database %in% one_tax_data$data_source_title) { # no database preference
      most_common <- names(sort(table(one_tax_data$matched_name), decreasing = TRUE)[1])
      if (is.null(most_common)) {
        return(NA_character_)
      } else {
        return(most_common)
      }
    } else {
      return(one_tax_data$matched_name[one_tax_data$data_source_title == gnr_database][1])
    }
  }

  vapply(names, pick_one, character(1))

}
