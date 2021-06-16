#' Taxonomy class
#'
#' Stores a taxonomy composed of [taxon()] objects organized in a tree
#' structure. This differs from the [hierarchies()] class in how the [taxon()]
#' objects are stored. Unlike [hierarchies()], each taxon is only stored once
#' and the relationships between taxa are stored in an [edge
#' list](https://en.wikipedia.org/wiki/Adjacency_list).
#'
#' @export
#' @param ... Any number of object of class [hierarchy()] or character
#'   vectors.
#' @param .list An alternate to the `...` input. Any number of object of class
#'   [hierarchy()] or character vectors in a list. Cannot be used with `...`.
#' @inheritParams parse_raw_heirarchies_to_taxonomy
#'
#' @return An `R6Class` object of class `Taxonomy`
#' @family classes
#'
#' @template taxonomyegs

taxonomy <- function(..., .list = NULL, named_by_rank = FALSE) {
  Taxonomy$new(..., .list = .list, named_by_rank = named_by_rank)
}

Taxonomy <- R6::R6Class(
  "Taxonomy",
  lock_class = TRUE,
  public = list(
    taxa = NULL,
    edge_list = NULL, # Note: this should be made of taxon ids, not indexes
    input_ids = NULL, # Only used by `Taxmap` right now

    # --------------------------------------------------------------------------
    # A simple wrapper to make future changes easier
    # it returns ids named by ids for consistency with other funcs
    taxon_ids = function() {
      stats::setNames(self$edge_list$to, self$edge_list$to)
    },

    # --------------------------------------------------------------------------
    # A simple wrapper to make future changes easier
    taxon_names = function() {
      vapply(self$taxa[self$taxon_ids()],
             function(x) {
               if (is.null(x$name$name)) {
                 return(NA_character_)
               } else {
                 return(x$name$name)
               }
             },
             character(1))
    },


    # --------------------------------------------------------------------------
    #  Set taxon names
    set_taxon_names = function(value) {
      if (length(value) !=  length(self$taxa))
      {
        stop(call. = FALSE, 'New taxon names must be of the same length as the number of taxa and in the same order as taxa.')
      }

      for (i in seq_len(length(self$taxa))) {
        self$taxa[[i]]$name$name <- as.character(value[[i]])
      }
    },


    # --------------------------------------------------------------------------
    # Return the taxon ranks in a taxonomy() or taxmap() object.
    # They are in the order taxa appear in the edge list.
    taxon_ranks = function() {
      vapply(self$taxa[self$taxon_ids()],
             function(x) {
               if (is.null(x$rank$name)) {
                 return(NA_character_)
               } else {
                 return(x$rank$name)
               }
             },
             character(1))
    },

    # --------------------------------------------------------------------------
    #  Set taxon names
    set_taxon_ranks = function(value) {
      if (length(value) !=  length(self$taxa))
      {
        stop(call. = FALSE, 'New taxon ranks must be of the same length as the number of taxa and in the same order as taxa.')
      }

      for (i in seq_len(length(self$taxa))) {
        self$taxa[[i]]$rank$name <- as.character(value[[i]])
      }
    },

    # --------------------------------------------------------------------------
    #  Set taxon authorities
    set_taxon_auths = function(value) {
      if (length(value) !=  length(self$taxa))
      {
        stop(call. = FALSE, 'New taxon authorities must be of the same length as the number of taxa and in the same order as taxa.')
      }

      for (i in seq_len(length(self$taxa))) {
        self$taxa[[i]]$authority <- as.character(value[[i]])
      }
    },

    # --------------------------------------------------------------------------
    # A simple wrapper to make future changes easier
    taxon_indexes = function() {
      stats::setNames(seq_len(nrow(self$edge_list)), self$taxon_ids())
    },

    # --------------------------------------------------------------------------
    # Constructor
    initialize = function(..., .list = NULL, named_by_rank = FALSE) {
      # Get intput
      input <- get_dots_or_list(..., .list = .list)

      # Parse input
      if (length(input) > 0 && "Hierarchy" %in% class(input[[1]])) {
        parsed_data <- parse_heirarchies_to_taxonomy(input)
      } else {
        parsed_data <- parse_raw_heirarchies_to_taxonomy(input, named_by_rank = named_by_rank)
      }
      self$taxa <- parsed_data$taxa
      self$edge_list <- parsed_data$edge_list
      self$input_ids <- parsed_data$input_ids

      # Convert numeric IDs to alpha
      if (length(self$taxa) > 0) {
        id_len <- as.integer(log(length(self$taxa), 25) + 1)
        self$replace_taxon_ids(convert_base(self$taxon_ids(),
                                            min_length = id_len))
      }
    },

    # --------------------------------------------------------------------------
    print = function(indent = "") {
      cat(paste0(indent, "<Taxonomy>\n"))
      taxon_names <- vapply(self$taxa, function(x) x$name$name, character(1))
      taxon_ids <- names(self$taxa)
      if (length(self$taxa) > 0) {
        limited_print(paste(tid_font(taxon_ids), taxon_names,
                            sep = punc_font(". ")),
                      sep = punc_font(", "),
                      mid = punc_font(" ... "),
                      trunc_char = punc_font("[truncated]"),
                      prefix = paste0(indent, "  ",
                                      length(self$taxa), " taxa:"),
                      type = "cat")
        limited_print(private$make_graph(),
                      sep = punc_font(", "),
                      mid = punc_font(" ... "),
                      trunc_char = punc_font("[truncated]"),
                      prefix = paste0(indent, "  ",
                                      nrow(self$edge_list), " edges:"),
                      type = "cat")
      } else {
        cat("Empty taxonomy")
      }
      invisible(self)
    },

    # --------------------------------------------------------------------------
    # Returns the names of things to be accessible using non-standard evaluation
    all_names = function() {
      output <- c()

      # Add functions included in the package
      output <- c(output, private$nse_accessible_funcs)

      # Add the name to the name of the name and return
      names(output) <- paste0(names(output),
                              ifelse(names(output) == "", "", "$"), output)
      return(output)
    },


    # --------------------------------------------------------------------------
    # Looks for names of data in a expression for use with NSE
    names_used = function(...) {
      decompose <- function(x) {
        if (any(class(x) %in% c("call", "(", "{")) && x[[1]] != "$") {
          return(lapply(1:length(x), function(i) decompose(x[[i]])))
        } else {
          return(deparse(x))
        }
      }

      expressions <- lapply(lazyeval::lazy_dots(...), function(x) x$expr)
      if (length(expressions) == 0) {
        return(character(0))
      } else {
        names_used <- unlist(lapply(1:length(expressions),
                                    function(i) decompose(expressions[[i]])))
        my_names <- self$all_names()
        return(my_names[my_names %in% names_used & my_names == make.names(my_names)])
      }
    },

    # --------------------------------------------------------------------------
    # Get data by name
    get_data = function(name = NULL, allow_external = TRUE, ...) {
      # Get default if name is NULL
      if (is.null(name)) {
        name = unique(self$all_names(...))
      }

      # Check that names provided are valid or an external object is used
      my_names <- self$all_names(...)
      if (any(unknown <- ! name %in% my_names)) {
        if (allow_external && (! is.null(obj_taxon_ids <- self$get_data_taxon_ids(name)))) {
          return(name)
        } else {
          stop(paste0("Cannot find the following data: ",
                      paste0(name[unknown], collapse = ", "), "\n ",
                      "Valid choices include: ",
                      paste0(my_names, collapse = ", "), "\n "))
        }
      }

      # Check for ambiguous terms
      name <- my_names[match(name, my_names)]
      name_counts <- table(my_names)
      ambiguous_names <- names(name_counts[name_counts > 1])
      ambiguous_names_used <- ambiguous_names[ambiguous_names %in% name]
      value_to_be_used <- unique(names(name)[name %in% ambiguous_names_used])
      if (length(ambiguous_names_used) > 0) {
        warning(call. = FALSE,
                'Ambiguous names used in non-standard evaluation:\n',
                limited_print(prefix = "  ", type = "silent", ambiguous_names_used),
                'These could refer to values in different datasets with the same name.',
                ' The following values will be used:\n',
                limited_print(prefix = "  ", type = "silent", value_to_be_used)
                )
      }

      # Format output
      output <- lapply(names(name),
                       function(x) eval(parse(text = paste0("self$", x))))
      names(output) <- name

      # Name each thing returned by taxon id if possibile
      #   This only applies to taxmap objects
      output[] <- lapply(seq_len(length(output)), function(index) {
        data_location <- names(name[index])
        if (startsWith(data_location, "data[['")) {
          data_name <- strsplit(data_location,
                                split =  "'", fixed = TRUE)[[1]][2]
          return(stats::setNames(output[[index]],
                                 self$get_data_taxon_ids(data_name)))
        } else {
          return(output[[index]])
        }
      })

      # Run any functions and return their results instead
      is_func <- vapply(output, is.function, logical(1))
      output[is_func] <- lapply(which(is_func), function(i) {
        if (length(formals(output[[i]])) > 0 && ! names(output[i]) %in% names(self)) {
          return(output[[i]](self))
        } else {
          return(output[[i]]())
        }
      })

      return(output)
    },

    # --------------------------------------------------------------------------
    # Get data in a taxonomy or taxmap object by name
    get_data_frame = function(...) {
      x <- self$get_data(...)
      if (length(unique(vapply(x, length, 1))) == 1) {
        tibble::as_tibble(data.frame(x, stringsAsFactors = FALSE))
      } else {
        stop("variables not of equal length")
      }
    },

    # --------------------------------------------------------------------------
    # Get a list of all data in an expression used with non-standard evaluation
    data_used = function(...) {
      my_names_used <- self$names_used(...)
      self$get_data(my_names_used)
    },


    # --------------------------------------------------------------------------
    # Return data for supertaxa (i.e. all taxa the target taxa are a part of)
    # of each taxon in a taxonomy() or taxmap() object.
    supertaxa = function(subset = NULL, recursive = TRUE, simplify = FALSE,
                         include_input = FALSE, value = "taxon_indexes", na = FALSE) {
      # non-standard argument evaluation
      data_used <- eval(substitute(self$data_used(subset)))
      subset <- rlang::eval_tidy(rlang::enquo(subset), data = data_used)
      subset <- private$parse_nse_taxon_subset(subset)

      # Get supertaxa
      parent_index <- match(self$edge_list$from, self$edge_list$to)
      recursive_part <- function(taxon, n_recursions) {
        supertaxon <- parent_index[taxon]
        if (n_recursions) {
          if (is.na(supertaxon)) {
            output <- c(taxon, supertaxon)
          } else {
            if (is.numeric(n_recursions)) {
              n_recursions <- n_recursions - 1
            }
            output <- c(taxon, recursive_part(supertaxon,
                                              n_recursions = n_recursions))
          }
        } else {
          output <- c(taxon, supertaxon)
        }
        return(unname(output))
      }

      if (is.numeric(recursive)) {
        n_recursions <- recursive - 1 # This makes 1 the same as FALSE
      } else {
        n_recursions <- recursive
      }

      if (is.numeric(recursive) && recursive == 0) {
        output <- setNames(lapply(subset, function(x) numeric(0)), subset)
      } else {
        output <- lapply(subset, recursive_part, n_recursions = n_recursions)
      }

      # Remove query taxa from output
      if (! include_input) {
        output <- lapply(output, `[`, -1)
      }

      # Remove NAs from output
      if (! na) {
        output <- lapply(output, function(x) x[!is.na(x)])
      }

      # Look up values
      if (!is.null(value)) {
        possible_values <- self$get_data(value)[[1]]
        if (is.null(names(possible_values))) {
          output <- lapply(output, function(i) possible_values[i])
        } else {
          output <- lapply(output, function(i) possible_values[self$taxon_ids()[i]])
        }
      }

      # Reduce dimensionality
      if (simplify) {
        output <- simplify(output)
      }

      return(output)
    },


    # --------------------------------------------------------------------------
    # Apply function to supertaxa of each taxon
    supertaxa_apply = function(func, subset = NULL, recursive = TRUE,
                               simplify = FALSE, include_input = FALSE,
                               value = "taxon_indexes", na = FALSE, ...) {
      my_sup <- eval(substitute(self$supertaxa(subset = subset,
                                               recursive = recursive,
                                               simplify = FALSE,
                                               include_input = include_input,
                                               value = value,
                                               na = na)))
      output <- lapply(my_sup, func, ...)
      if (simplify) {
        output <- simplify(output)
      }
      return(output)
    },


    # --------------------------------------------------------------------------
    # Return the root taxa for a taxonomy() or taxmap() object.
    # Can also be used to get the roots of a subset of taxa.
    roots = function(subset = NULL, value = "taxon_indexes") {
      # non-standard argument evaluation
      data_used <- eval(substitute(self$data_used(subset)))
      subset <- rlang::eval_tidy(rlang::enquo(subset), data = data_used)
      subset <- private$parse_nse_taxon_subset(subset)

      # Get roots
      parents <- self$supertaxa(subset = subset, recursive = TRUE,
                                include_input = TRUE, value = "taxon_indexes",
                                na = FALSE)
      is_global_root <- vapply(parents, length, numeric(1)) == 1
      if (missing(subset)) {
        is_root <- is_global_root
      } else {
        is_root <- is_global_root | vapply(parents,
                                           FUN.VALUE = logical(1),
                                           function(x) ! any(x[-1] %in% subset))
      }
      output <- unname(subset[is_root])

      # Look up values
      if (! is.null(value)) {
        possible_values <- self$get_data(value)[[1]]
        if (is.null(names(possible_values))) {
          output <- possible_values[output]
        } else {
          output <- possible_values[self$taxon_ids()[output]]
        }
      } else {
        names(output) <- self$taxon_ids()[output]
      }

      return(output)
    },


    # --------------------------------------------------------------------------
    # Return the stem taxa for a taxonomy() or a taxmap() object.
    # Stem taxa are all those from the roots to the first taxon with more than one subtaxon.
    stems = function(subset = NULL, value = "taxon_indexes", simplify = FALSE,
                     exclude_leaves = FALSE) {
      # non-standard argument evaluation
      data_used <- eval(substitute(self$data_used(subset)))
      subset <- rlang::eval_tidy(rlang::enquo(subset), data = data_used)
      subset <- private$parse_nse_taxon_subset(subset)

      # Get roots to start search
      my_roots <- self$roots(subset = subset, value = "taxon_indexes")

      # Search until taxa with multiple subtaxa are found
      parent_index <- match(self$edge_list$from, self$edge_list$to)
      recursive_part <- function(taxon) {
        children <- which(parent_index == taxon)
        if (length(children) == 0 && ! exclude_leaves) {
          output <- taxon
        } else if (length(children) == 1) {
          output <- c(taxon, recursive_part(children))
        } else {
          output <- numeric(0)
        }
        return(unname(output))
      }
      output <- lapply(my_roots, recursive_part)

      # Look up values
      if (!is.null(value)) {
        possible_values <- self$get_data(value)[[1]]
        if (is.null(names(possible_values))) {
          output <- lapply(output, function(i) possible_values[i])
        } else {
          output <- lapply(output, function(i) possible_values[self$taxon_ids()[i]])
        }
      }

      # Reduce dimensionality
      if (simplify) {
        output <- simplify(output)
      }

      return(output)
    },


    # --------------------------------------------------------------------------
    # Leaf taxa are taxa with no subtaxa.
    leaves = function(subset = NULL, recursive = TRUE, simplify = FALSE, value = "taxon_indexes") {
      # Find taxa without subtaxa (leaves)
      childless_taxa <- which(self$n_subtaxa_1() == 0)

      # Subset subtaxa results to just leaves
      eval(substitute(self$subtaxa_apply(func = function(x) x[names(x) %in% names(childless_taxa)],
                                         subset = subset,
                                         simplify = simplify,
                                         recursive = recursive,
                                         include_input = FALSE,
                                         value = value)))
    },


    # --------------------------------------------------------------------------
    # Apply a function to the leaves of each taxon.
    # This is similar to using leaves() with lapply() or sapply().
    leaves_apply = function(func, subset = NULL, recursive = TRUE, simplify = FALSE,
                            value = "taxon_indexes", ...) {
      my_sub <- eval(substitute(self$leaves(subset = subset,
                                            recursive = recursive,
                                            simplify = FALSE,
                                            value = value)))
      output <- lapply(my_sub, func, ...)
      if (simplify) {
        output <- simplify(output)
      }
      return(output)
    },


    # --------------------------------------------------------------------------
    # A branch is anything that is not a root, stem, or leaf.
    # Its the interior of the tree after the first split starting from the roots
    branches = function(subset = NULL, value = "taxon_indexes") {
      # non-standard argument evaluation
      data_used <- eval(substitute(self$data_used(subset)))
      subset <- rlang::eval_tidy(rlang::enquo(subset), data = data_used)
      subset <- private$parse_nse_taxon_subset(subset)

      # Return empty list if `subset` has no values
      if (length(subset) == 0) {
        return(vector(mode = class(self$get_data(value)[[1]])))
      }

      # Prefilter subset
      if (length(subset) == length(self$taxa) && all(subset == seq_len(length(self$taxon_ids())))) {
        filtered_self <- self
      } else {
        filtered_self <- filter_taxa(self, subset)
      }

      # Get branches
      output <- which(filtered_self$is_branch())

      # Look up values
      if (!is.null(value)) {
        possible_values <- self$get_data(value)[[1]]
        if (is.null(names(possible_values))) {
          output <- possible_values[output]
        } else {
          output <- possible_values[self$taxon_ids()[output]]
        }
      }

      return(output)
    },


    # --------------------------------------------------------------------------
    # An internode is any taxon with a single immediate supertaxon and a single immediate subtaxon.
    internodes = function(subset = NULL, value = "taxon_indexes") {
      # non-standard argument evaluation
      data_used <- eval(substitute(self$data_used(subset)))
      subset <- rlang::eval_tidy(rlang::enquo(subset), data = data_used)
      subset <- private$parse_nse_taxon_subset(subset)

      # Return empty list if `subset` has no values
      if (length(subset) == 0) {
        return(vector(mode = class(self$get_data(value)[[1]])))
      }

      # Prefilter subset
      if (length(subset) == length(self$taxa) && all(subset == seq_len(length(self$taxon_ids())))) {
        filtered_self <- self
      } else {
        filtered_self <- filter_taxa(self, subset)
      }

      # Get branches
      output <- which(filtered_self$is_internode())

      # Look up values
      if (!is.null(value)) {
        possible_values <- self$get_data(value)[[1]]
        if (is.null(names(possible_values))) {
          output <- possible_values[output]
        } else {
          output <- possible_values[self$taxon_ids()[output]]
        }
      }

      return(output)
    },


    # --------------------------------------------------------------------------
    # Return data for the subtaxa of each taxon in an taxonomy() or taxmap() object.
    subtaxa = function(subset = NULL, recursive = TRUE,
                       simplify = FALSE, include_input = FALSE,
                       value = "taxon_indexes") {
      # non-standard argument evaluation
      data_used <- eval(substitute(self$data_used(subset)))
      subset <- rlang::eval_tidy(rlang::enquo(subset), data = data_used)
      subset <- private$parse_nse_taxon_subset(subset)

      # Return empty list if `subset` has no values
      if (length(subset) == 0) {
        if (simplify) {
          return(vector(mode = class( self$get_data(value)[[1]])))
        } else {
          return(list())
        }
      }

      # Get subtaxa
      parent_index <- match(self$edge_list$from, self$edge_list$to)

      get_children <- function(taxon) {
        which(parent_index == taxon)
      }

      recursive_part <- function(taxon) {
        # Get immediate children of current taxon
        children <- get_children(taxon)
        # Run this function on them to get their output
        child_output <- lapply(children, recursive_part) # stops if no children
        child_output <- stats::setNames(unlist(child_output, recursive = FALSE),
                                        unlist(lapply(child_output, names)))
        # Get all subtaxa from the names of the child output
        child_taxa <- c(taxon, as.numeric(names(child_output)))
        # Combine the child output with the subtaxa for the current taxon
        output <- stats::setNames(c(list(child_taxa), child_output),
                                  c(taxon, names(child_output)))
        return(output)
      }

      if (recursive) {
        starting_taxa <- unname(self$roots(subset = subset,
                                           value = "taxon_indexes"))
        output <- stats::setNames(
          unlist(lapply(starting_taxa, recursive_part),
                 recursive = FALSE)[as.character(subset)],
          names(subset)
        )
      } else {
        output <- lapply(subset, function(x) c(x, get_children(x)))
      }

      # Remove query taxa from output
      if (! include_input) {
        output <- lapply(output, `[`, -1)
      }

      # Simulate limited recursion
      #
      # To increase speed, the recursive algorithm only searches starting at
      # root taxa, but this makes it hard to limit the number of rankes returned
      # below each taxon during recursion. Instead, a finite number of
      # recursions are simulated by filtering the results of tarversing the
      # entire tree and comparing rank depth between each taxon and its subtaxa.
      if (is.numeric(recursive) && recursive >= 0) {
        all_taxa <- unique(c(self$map_data(taxon_ids, taxon_indexes)[names(output)],
                             unlist(output)))
        rank_depth <- vapply(self$supertaxa(all_taxa), length, numeric(1))
        output_names <- names(output)
        output <- lapply(seq_along(output), function(i) {
          subtaxa_ids <- self$taxon_ids()[output[[i]]]
          subtaxa_depth <- rank_depth[subtaxa_ids]
          current_depth <- rank_depth[names(output[i])]
          passing_taxa <- subtaxa_depth - current_depth <= recursive
          return(output[[i]][passing_taxa])
        })
        names(output) <- output_names
      }

      # Look up values
      if (!is.null(value)) {
        possible_values <- self$get_data(value)[[1]]
        if (is.null(names(possible_values))) {
          output <- lapply(output, function(i) possible_values[i])
        } else {
          output <- lapply(output, function(i) possible_values[self$taxon_ids()[i]])
        }
      }

      # Reduce dimensionality
      if (simplify) {
        output <- simplify(output)
      }

      return(output)
    },


    # --------------------------------------------------------------------------
    # Apply a function to the subtaxa for each taxon.
    # This is similar to using subtaxa() with lapply() or sapply().
    subtaxa_apply = function(func, subset = NULL, recursive = TRUE,
                             simplify = FALSE, include_input = FALSE,
                             value = "taxon_indexes", ...) {
      my_sub <- eval(substitute(self$subtaxa(subset = subset,
                                             recursive = recursive,
                                             simplify = FALSE,
                                             include_input = include_input,
                                             value = value)))
      output <- lapply(my_sub, func, ...)
      if (simplify) {
        output <- simplify(output)
      }
      return(output)
    },


    # --------------------------------------------------------------------------
    # Get classifications of taxa
    classifications = function(value = "taxon_names", sep = ";") {
      vapply(self$supertaxa(recursive = TRUE, include_input = TRUE,
                            value = value, na = FALSE),
             function(x) paste0(rev(x), collapse = sep), character(1))
    },

    # --------------------------------------------------------------------------
    # Get ID classifications of taxa
    id_classifications = function(sep = ";") {
      self$classifications(value = "taxon_ids", sep = sep)
    },

    # --------------------------------------------------------------------------
    # Get number of supertaxa for each taxon
    n_supertaxa = function() {
      vapply(self$supertaxa(recursive = TRUE, include_input = FALSE,
                            value = "taxon_indexes", na = FALSE),
             length, numeric(1))
    },

    # --------------------------------------------------------------------------
    # Get number of immediate supertaxa (not supertaxa of supertaxa) for each taxon
    n_supertaxa_1 = function() {
      vapply(self$supertaxa(recursive = FALSE, include_input = FALSE,
                            value = "taxon_indexes", na = FALSE),
             length, numeric(1))
    },

    # --------------------------------------------------------------------------
    # Get number of subtaxa for each taxon
    n_subtaxa = function() {
      vapply(self$subtaxa(recursive = TRUE, include_input = FALSE,
                          value = "taxon_indexes"),
             length, numeric(1))
    },

    # --------------------------------------------------------------------------
    # Get number of subtaxa for each taxon, not including subtaxa of subtaxa
    n_subtaxa_1 = function() {
      vapply(self$subtaxa(recursive = FALSE, include_input = FALSE,
                          value = "taxon_indexes"),
             length, numeric(1))
    },

    # --------------------------------------------------------------------------
    # Get number of leaves for each taxon
    n_leaves = function() {
      vapply(self$leaves(recursive = TRUE, value = "taxon_indexes"), length, numeric(1))
    },

    # --------------------------------------------------------------------------
    # Get number of leaves for each taxon, not including leaves of subtaxa etc.
    n_leaves_1 = function() {
      vapply(self$leaves(recursive = FALSE, value = "taxon_indexes"), length, numeric(1))
    },

    # --------------------------------------------------------------------------
    # Test if taxa are roots
    is_root = function() {
      stats::setNames(is.na(self$edge_list$from), self$taxon_ids())
    },

    # --------------------------------------------------------------------------
    # Test if taxa are leaves
    is_leaf = function() {
      self$n_subtaxa() == 0
    },

    # --------------------------------------------------------------------------
    # Test if taxa are stems
    is_stem = function() {
      stats::setNames(self$taxon_ids() %in% self$stems(simplify = TRUE,
                                                       value = "taxon_ids"),
                      self$taxon_ids())
    },

    # --------------------------------------------------------------------------
    # Test if taxa are branches
    is_branch = function() {
      stats::setNames(! (self$is_root() | self$is_leaf() | self$is_stem()),
                      self$taxon_ids())
    },

    # --------------------------------------------------------------------------
    # Test if taxa are "internodes"
    is_internode = function() {
      stats::setNames(self$n_subtaxa_1() == 1 & self$n_supertaxa_1() == 1,
                      self$taxon_ids())
    },

    # --------------------------------------------------------------------------
    # Filter taxa in a taxonomy() or taxmap() object with a series of conditions.
    filter_taxa = function(..., subtaxa = FALSE, supertaxa = FALSE,
                           drop_obs = TRUE, reassign_obs = TRUE,
                           reassign_taxa = TRUE, invert = FALSE,
                           keep_order = TRUE) {

      # non-standard argument evaluation
      selection <- private$parse_nse_taxon_subset(...)

      # Get taxa of subset
      if (is.logical(subtaxa) && subtaxa == FALSE) {
        subtaxa = 0
      }
      if (is.logical(supertaxa) && supertaxa == FALSE) {
        supertaxa = 0
      }
      taxa_subset <- unique(c(selection,
                              self$subtaxa(subset = selection,
                                           recursive = subtaxa,
                                           value = "taxon_indexes",
                                           include_input = FALSE,
                                           simplify = TRUE),
                              self$supertaxa(subset = selection,
                                             recursive = supertaxa,
                                             value = "taxon_indexes",
                                             na = FALSE, simplify = TRUE,
                                             include_input = FALSE)
      ))

      # Preserve original order
      if (keep_order) {
        taxa_subset <- sort(taxa_subset)
      }

      # Invert selection
      if (invert) {
        taxa_subset <- (1:nrow(self$edge_list))[-taxa_subset]
      }

      # Reassign taxonless observations
      if ("Taxmap" %in% class(self)) {
        reassign_obs <- parse_possibly_named_logical(
          reassign_obs,
          self$data,
          default <- formals(self$filter_taxa)$reassign_obs
        )
        process_one <- function(data_index) {

          reassign_one <- function(parents) {
            included_parents <- parents[parents %in% taxa_subset]
            return(self$taxon_ids()[included_parents[1]])
          }

          # Get the taxon ids of the current object
          if (is.null((data_taxon_ids <-
                       self$get_data_taxon_ids(data_index, warn = TRUE)))) {
            return(NULL) # if there is no taxon id info, dont change anything
          }

          # Generate replacement taxon ids
          to_reassign <- ! data_taxon_ids %in% self$taxon_ids()[taxa_subset]
          supertaxa_key <- self$supertaxa(
            subset = unique(data_taxon_ids[to_reassign]),
            recursive = TRUE, simplify = FALSE, include_input = FALSE,
            value = "taxon_indexes", na = FALSE
          )
          reassign_key <- vapply(supertaxa_key, reassign_one, character(1))
          new_data_taxon_ids <- reassign_key[data_taxon_ids[to_reassign]]

          # Replace taxon ids
          if (is.data.frame(self$data[[data_index]])) {
            self$data[[data_index]][to_reassign, "taxon_id"] <- new_data_taxon_ids
          } else {
            names(self$data[[data_index]])[to_reassign] <- new_data_taxon_ids
          }
        }

        unused_output <- lapply(seq_along(self$data)[reassign_obs], process_one)
      }

      # Reassign subtaxa
      if (reassign_taxa) {
        reassign_one <- function(parents) {
          included_parents <- parents[parents %in% taxa_subset]
          return(self$taxon_ids()[included_parents[1]])
        }

        to_reassign <- ! self$edge_list$from %in% self$taxon_ids()[taxa_subset]
        supertaxa_key <- self$supertaxa(
          subset = unique(self$taxon_ids()[to_reassign]),
          recursive = TRUE, simplify = FALSE, include_input = FALSE,
          value = "taxon_indexes", na = FALSE)
        reassign_key <- vapply(supertaxa_key, reassign_one, character(1)
        )
        self$edge_list[to_reassign, "from"] <-
          reassign_key[self$taxon_ids()[to_reassign]]
      }

      # Remove taxonless observations
      if ("Taxmap" %in% class(self)) {
        drop_obs <- parse_possibly_named_logical(
          drop_obs,
          self$data,
          default = formals(self$filter_taxa)$drop_obs
        )
        process_one <- function(my_index) {

          # Get the taxon ids of the current object
          if (is.null((data_taxon_ids <-
                       self$get_data_taxon_ids(my_index)))) {
            return(NULL) # if there is no taxon id info, dont change anything
          }

          obs_subset <- data_taxon_ids %in% self$taxon_ids()[taxa_subset]
          private$remove_obs(data = my_index,
                             indexes = obs_subset,
                             unname_only = ! drop_obs[my_index])
        }
        unused_output <- lapply(seq_along(self$data), process_one)
      }

      # Remove filtered taxa
      private$remove_taxa(taxa_subset)

      return(self)
    },

    # --------------------------------------------------------------------------
    # Sort the edge list and taxon list
    arrange_taxa = function(...) {
      # Sort edge list
      data_used <- self$data_used(...)
      data_used <- data_used[! names(data_used) %in% names(self$edge_list)]
      if (length(data_used) == 0) {
        self$edge_list <- dplyr::arrange(self$edge_list, ...)
      } else {
        target_with_extra_cols <- dplyr::bind_cols(data_used, self$edge_list)
        self$edge_list <-
          dplyr::arrange(target_with_extra_cols, ...)[, -seq_along(data_used)]
      }

      # Reorder taxa list to be the same order
      self$taxa <- self$taxa[self$edge_list$to]

      return(self)
    },

    # --------------------------------------------------------------------------
    # Randomly sample some number of taxa
    sample_n_taxa = function(size, taxon_weight = NULL, obs_weight = NULL,
                             obs_target = NULL, use_subtaxa = TRUE,
                             collapse_func = mean, ...) {

      # non-standard argument evaluation
      data_used <- eval(substitute(self$data_used(taxon_weight, obs_weight)))
      taxon_weight <- rlang::eval_tidy(rlang::enquo(taxon_weight),
                                       data = data_used)
      obs_weight <- rlang::eval_tidy(rlang::enquo(obs_weight),
                                     data = data_used)

      # Calculate observation component of taxon weights
      if (is.null(obs_weight) || !"Taxmap" %in% class(self)) {
        taxon_obs_weight <- rep(1, nrow(self$edge_list))
      } else {
        if (is.null(obs_target)) {
          stop(paste("If the option `obs_weight` is used, then `obs_target`",
                     "must also be defined."))
        }
        my_obs <- self$obs(obs_target, recursive = use_subtaxa,
                           simplify = FALSE)
        taxon_obs_weight <- vapply(my_obs,
                                   function(x) collapse_func(obs_weight[x]),
                                   numeric(1))
      }
      taxon_obs_weight <- taxon_obs_weight / sum(taxon_obs_weight)

      # Calculate taxon component of taxon weights
      if (is.null(taxon_weight)) {
        taxon_weight <- rep(1, nrow(self$edge_list))
      }
      taxon_weight <- taxon_weight / sum(taxon_weight)

      # Combine observation and taxon weight components
      combine_func <- prod
      weight <- mapply(taxon_weight, taxon_obs_weight,
                       FUN = function(x, y) combine_func(c(x,y)))
      weight <- weight / sum(weight)

      # Sample
      sampled_rows <- sample.int(nrow(self$edge_list), size = size,
                                 replace = FALSE, prob = weight)
      self$filter_taxa(sampled_rows, ...)
    },

    # --------------------------------------------------------------------------
    # Randomly sample some proportion of taxa
    sample_frac_taxa = function(size = 1, taxon_weight = NULL,
                                obs_weight = NULL, obs_target = NULL,
                                use_subtaxa = TRUE, collapse_func = mean, ...) {
      self$sample_n_taxa(size = size * nrow(self$edge_list),
                         taxon_weight = taxon_weight,
                         obs_weight = obs_weight, obs_target = obs_target,
                         use_subtaxa = use_subtaxa,
                         collapse_func = collapse_func, ...)
    },


    # --------------------------------------------------------------------------
    # Creates a named vector that maps the values of two variables associated with taxa
    map_data = function(from, to, warn = TRUE) {
      # non-standard argument evaluation
      data_used <- eval(substitute(self$data_used(from, to)))
      # to_data <- rlang::eval_tidy(rlang::enquo(to), data = data_used)
      # from_data <- rlang::eval_tidy(rlang::enquo(from), data = data_used)

      # check that arguments have taxon ids and evaluate
      validate_and_eval <- function(unparsed) {
        parsed <- rlang::eval_tidy(rlang::enquo(unparsed),
                                   data = data_used)
        if (! private$valid_taxon_ids(names(parsed))) {
          stop(paste0("The value `", deparse(match.call()$unparsed),
                      "` is not named by taxon id or contains invalid ids. ",
                      "Use `taxon_ids()` to see the valid ids. ",
                      "Use `warn = FALSE` to ignore this."))
        }
        return(parsed)
      }
      to_data <- eval(substitute(validate_and_eval(to)))
      from_data <- eval(substitute(validate_and_eval(from)))

      # Check for multiple different values of `to` for each `from`
      is_one_to_one <- vapply(unique(names(from_data)),
                              function(n) {
                                matches <- to_data[names(to_data)==n]
                                matches <- matches[!is.na(matches)]
                                length(unique(matches)) <= 1
                              },
                              logical(1))
      if (warn && any(! is_one_to_one)) {
        warning(paste0('There are multiple unique values of "',
                       deparse(match.call()$to),
                       '" for at least one value of "',
                       deparse(match.call()$from),
                       '". Only the first instance will be returned. ',
                       'To get all instances, use the `obs` function.'))
      }

      # Map values using taxon ids
      self$map_data_(from = from_data, to = to_data)
    },

    # --------------------------------------------------------------------------
    # map_data without NSE
    map_data_ = function(from, to) {
      stats::setNames(to[match(names(from), names(to))], from)
    },

    # --------------------------------------------------------------------------
    # Replace taxon ids
    replace_taxon_ids = function(new_ids) {
      # Check that new ids are unique
      duplicate_ids <- unique(new_ids[duplicated(new_ids)])
      if (any(duplicated(new_ids))) {
        stop(paste0("New taxon IDs must be unique. ",
                    "The following ", length(duplicate_ids),
                    " taxon ids are not unique:\n",
                    limited_print(duplicate_ids, type = "silent")))
      }

      # Check that new ids are the same length as old ids
      if (length(new_ids) != length(self$taxon_ids())) {
        stop(paste0('The number of new taxon IDs (', length(new_ids),
                    ') is different than the current number of taxa (',
                    length(self$taxon_ids()), ').'))
      }

      # Replace taxon ids in datasets
      names(new_ids) <- self$taxon_ids()
      if (!is.null(self$data)) {
        self$data <- lapply(self$data, function(x) {
          if (is.data.frame(x) && "taxon_id" %in% colnames(x) && private$valid_taxon_ids(x$taxon_id)) {
            x$taxon_id <- new_ids[x$taxon_id]
          } else if (!is.null(names(x)) && private$valid_taxon_ids(names(x))) {
            names(x) <- new_ids[names(x)]
          }
          return(x)
        })
      }

      # Replace taxon ids
      if (length(self$taxa) > 0) {
        self$input_ids <- new_ids[self$input_ids]
        self$edge_list$to <- new_ids[self$edge_list$to]
        self$edge_list$from <- new_ids[self$edge_list$from]
        names(self$taxa) <- new_ids[names(self$taxa)]
      }

      # Return modified object
      return(self)
    },

    # --------------------------------------------------------------------------
    # Remove the names of parent taxa in the begining of their children's names
    remove_redundant_names = function() {
      new_names <- vapply(supertaxa(self, recursive = FALSE, include_input = TRUE),
                          function(x) gsub(self$taxon_names()[x[1]],
                                           pattern = paste0("^", self$taxon_names()[x[2]], "[_ ]+"),
                                           replacement = ""),
                          character(1))
      lapply(seq_along(new_names), function(i) {
        self$taxa[[i]]$name$name <- new_names[i]
      })
      return(self)
    },

    # pop = function(ranks = NULL, names = NULL, ids = NULL) {
    #   taxa_rks <- vapply(self$taxa, function(x) x$rank$name, "")
    #   taxa_nms <- vapply(self$taxa, function(x) x$name$name, "")
    #   taxa_ids <- vapply(self$taxa, function(x) x$id$id, numeric(1))
    #   todrop <- which(taxa_rks %in% ranks |
    #                     taxa_nms %in% names |
    #                     taxa_ids %in% ids)
    #   tmp_taxa <- self$taxa
    #   tmp_taxa[todrop] <- NULL
    #   parsed_data <- parse_heirarchies_to_taxonomy(hierarchy(tmp_taxa))
    #   self$taxa <- parsed_data$taxa
    #   self$edge_list <- parsed_data$edge_list
    #   self$input_ids <- parsed_data$input_ids
    #   return(self)
    # }

    # --------------------------------------------------------------------------
    # Return taxonomy information in a taxon x rank table
    taxonomy_table = function(subset = NULL, value = "taxon_names",
                              use_ranks = NULL, add_id_col = FALSE) {

      # non-standard argument evaluation of subset
      data_used <- eval(substitute(self$data_used(subset)))
      subset <- rlang::eval_tidy(rlang::enquo(subset), data = data_used)
      if (is.null(subset)) {
        subset <- self$leaves(simplify = TRUE)
      }
      subset <- private$parse_nse_taxon_subset(subset)

      # Get supertaxa of each taxon, named by ranks
      obj_has_ranks <- ! all(is.na(self$taxon_ranks()))
      class_list <- self$supertaxa(subset = subset, value = value,
                              include_input = TRUE)
      if (is.null(use_ranks)) {
        if (obj_has_ranks) {
          rank_names <- self$supertaxa(subset = subset, value = "taxon_ranks",
                                  include_input = TRUE)
          ranks_used <- rev(rank_names[[which.max(vapply(rank_names, length, numeric(1)))]])
        } else {
          rank_names <- self$supertaxa(subset = subset, value = "n_supertaxa",
                                  include_input = TRUE)
          rank_names <- lapply(rank_names, function(x) paste0("rank_", x + 1))
          ranks_used <- rev(rank_names[[which.max(vapply(rank_names, length, numeric(1)))]])
        }
      } else if (is.logical(use_ranks)) {
        if (use_ranks == TRUE) {
          if (obj_has_ranks) {
            rank_names <- self$supertaxa(subset = subset, value = "taxon_ranks",
                                    include_input = TRUE)
            ranks_used <- rev(rank_names[[which.max(vapply(rank_names, length, numeric(1)))]])
          } else {
            stop(call. = FALSE, "option `use_ranks is `TRUE`, but there is no rank information.")
          }
        } else {
          rank_names <- self$supertaxa(subset = subset, value = "n_supertaxa",
                                  include_input = TRUE)
          rank_names <- lapply(rank_names, function(x) paste0("rank_", x + 1))
          ranks_used <- rev(rank_names[[which.max(vapply(rank_names, length, numeric(1)))]])
        }
      } else if (is.character(use_ranks)) {
        rank_names <- self$supertaxa(subset = subset, value = "taxon_ranks",
                                include_input = TRUE)
        ranks_used <- use_ranks
      } else if (is.numeric(use_ranks)) {
        rank_names <- self$supertaxa(subset = subset, value = "n_supertaxa",
                                include_input = TRUE)
        rank_names <- lapply(rank_names, function(x) paste0("rank_", x + 1))
        ranks_used <-  rev(paste0("rank_", use_ranks))
      } else {
        stop(call. = FALSE, "Invalid `use_ranks` input. See ?taxonomy_table for valid inputs.")
      }
      class_list <- lapply(seq_len(length(class_list)),
                           function(i) stats::setNames(class_list[[i]], rank_names[[i]]))

      # Check for unused ranks
      all_ranks <- unique(unlist(rank_names))
      unused_ranks <- all_ranks[! all_ranks %in% ranks_used]
      if (length(unused_ranks) > 0) {
        message('The following ranks will not be included because the order cannot be determined:\n',
                limited_print(unused_ranks, prefix = "  ", type = "silent"),
                'See the section on the `use_ranks` option in ?taxonomy_table if you want to change which ranks are used.')
      }

      # Make table
      output <- do.call(rbind, lapply(class_list, function(x) x[ranks_used]))
      colnames(output) <- ranks_used

      # Add taxon ID column
      if (add_id_col) {
        output <- cbind(data.frame(stringsAsFactors = FALSE,
                                   taxon_id = self$taxon_ids()[subset]),
                        output)
      }

      return(dplyr::as_tibble(output))
    },

    # Make text print out of a tree
    print_tree = function(value = "taxon_names") {
        # Make tree with taxon IDs to avoid potential loop errors
        tree_data <- data.frame(stringsAsFactors = FALSE,
                                taxa = self$taxon_ids(),
                                subtaxa = I(self$subtaxa(recursive = FALSE,
                                                         value = "taxon_ids")))
        trees <- lapply(self$roots(value = "taxon_ids"), function(my_root) {
          cli::tree(tree_data, root = my_root)
        })
        trees <- unlist(trees)

        # Replace taxon ids with other value
        if (value != "taxon_ids") {
          value <- self$get_data(value)[[1]]
          ids_in_tree <- sub(trees, pattern = "^[\u2502 \u251C\u2500\u2514]*", replacement = "")
          trees <- vapply(seq_len(length(trees)), FUN.VALUE = character(1),
                          function(i) {
                            sub(trees[i], pattern = ids_in_tree[i],
                                replacement = value[ids_in_tree[i]])

                          })
        }

        # Return tree
        # cat(paste0(trees, collapse = "\n"))
        class(trees) <- unique(c("tree", "character"))
        return(trees)
    }

  ),

  private = list(
    # --------------------------------------------------------------------------
    nse_accessible_funcs = c("taxon_names",
                             "taxon_ids",
                             "taxon_indexes",
                             "classifications",
                             "n_supertaxa",
                             "n_supertaxa_1",
                             "n_subtaxa",
                             "n_subtaxa_1",
                             "n_leaves",
                             "n_leaves_1",
                             "taxon_ranks",
                             "is_root",
                             "is_stem",
                             "is_branch",
                             "is_leaf",
                             "is_internode"),

    # --------------------------------------------------------------------------
    make_graph = function() {
      apply(self$edge_list, 1, function(x) paste0(tid_font(x), collapse = punc_font("->")))
    },

    # --------------------------------------------------------------------------
    # Remove taxa NOT in "el_indexes"
    remove_taxa = function(el_indexes) {
      # Remove taxa objects
      self$taxa <- self$taxa[self$taxon_ids()[el_indexes]]

      # Remove corresponding rows in the edge list
      self$edge_list <- self$edge_list[el_indexes, , drop = FALSE]

      # Replace and edges going to removed taxa with NA
      to_replace <- ! self$edge_list$from %in% self$taxon_ids()
      if (sum(to_replace) > 0) {
        self$edge_list[to_replace, "from"] <- as.character(NA)
      }
    },

    # --------------------------------------------------------------------------
    # Takes one ore more NSE expressions and resolves them to indexes of edgelist rows
    # Each expression can resolve to taxon ids, edgelist indexes, or logical.
    parse_nse_taxon_subset = function(...) {
      # Non-standard argument evaluation
      selection <- lapply(rlang::quos(...), rlang::eval_tidy, data = self$data_used(...))

      # Remove any NULL selections
      is_blank <- vapply(selection, is.null, logical(1))
      selection <- selection[! is_blank]

      # Default to all taxa if no selection is provided
      if (length(selection) == 0) {
        return(self$taxon_indexes())
      }

      # Check that index input is valid
      is_num <- vapply(selection, is.numeric, logical(1))
      unused <- lapply(selection[is_num],
                       function(x) {
                         invalid_indexes <- x[x < 1 | x > length(self$taxon_ids())]
                         if (length(invalid_indexes) > 0) {
                           stop(call. = FALSE,
                                paste0("The following taxon indexes are invalid:\n",
                                       limited_print(invalid_indexes, type = "silent")))
                         }
                       })

      # Convert taxon_ids to indexes
      is_char <- vapply(selection, is.character, logical(1))
      selection[is_char] <- lapply(selection[is_char],
                                   function(x) {
                                     result <- match(x, self$taxon_ids())
                                     invalid_ids <- x[is.na(result)]
                                     if (length(invalid_ids) > 0) {
                                       stop(call. = FALSE,
                                            paste0("The following taxon IDs do not exist:\n",
                                                   limited_print(invalid_ids, type = "silent")))
                                     }
                                     return(result)
                                   })

      # Convert logical to indexes
      is_tf <- vapply(selection, is.logical, logical(1))
      selection[is_tf] <- lapply(selection[is_tf],
                                 function(x) {
                                   if (length(x) != length(self$taxon_ids())) {
                                     stop(call. = FALSE,
                                          paste0("TRUE/FALSE vector (length = ",
                                                 length(x),
                                                 ") must be the same length as the number of taxa (",
                                                 length(self$taxon_ids()), ")"))
                                   }
                                   which(x)
                                 })

      # Combine index lists.
      intersect_with_dups <- function(a, b) {
        #taken from http://r.789695.n4.nabble.com/intersect-without-discarding-duplicates-td2225377.html
        rep(sort(intersect(a, b)), pmin(table(a[a %in% b]), table(b[b %in% a])))
      }
      output <- Reduce(intersect_with_dups, selection)

      # Name by taxon id
      names(output) <- self$taxon_ids()[output]

      return(output)
    },

    # --------------------------------------------------------------------------
    # check if a set of putative taxon ids are valid.
    # returns TRUE/FALSE
    valid_taxon_ids = function(ids) {
      !is.null(ids) && all(ids %in% self$taxon_ids() | is.na(ids))
    }
  )
)
