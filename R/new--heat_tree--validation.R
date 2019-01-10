#' Validate arguments given to heat_tree
#'
#' Validate arguments given to \code{\link{heat_tree}} or
#' \code{\link{heat_tree_data}} and issue warnings or errors if things dont make
#' sense.
#' 
#' @param obj taxmap object
#' @param args A named list of argument values
#' 
#' @return A named list of argument values, potentially modified
#' 
#' @keywords internal
heat_tree_validate_arguments <- function(obj, args)
{
  
  # Check that inputs are named by taxon ID or are the same length as the number of taxa
  args <- check_length_or_id(obj,
                             args[c("node_size", 
                                    "edge_size",
                                    "node_label_size", 
                                    "edge_label_size",  
                                    "tree_label_size",
                                    "node_color", 
                                    "edge_color", 
                                    "tree_color",
                                    "node_label_color", 
                                    "edge_label_color", 
                                    "tree_label_color",
                                    "node_label", 
                                    "edge_label", 
                                    "tree_label")])
  args <- check_for_na(args[c("node_size",
                              "edge_size",
                              "node_label_size",
                              "edge_label_size",
                              "tree_label_size",
                              "node_color", 
                              "edge_color", 
                              "tree_color",
                              "node_label_color",
                              "edge_label_color", 
                              "tree_label_color")],
                       warn = TRUE)
  args <- check_for_na(args[c("node_label",
                              "edge_label",
                              "tree_label")],
                       warn = FALSE)
  args <- check_size(args[c("node_size",
                            "edge_size",
                            "node_label_size", 
                            "edge_label_size", 
                            "tree_label_size")])
  args <- check_size_range(args[c("node_size_range",  
                                  "edge_size_range", 
                                  "node_label_size_range",
                                  "edge_label_size_range", 
                                  "tree_label_size_range",
                                  "node_size_interval", 
                                  "edge_size_interval")])
  args <- check_trans(args[c("node_size_trans", 
                             "edge_size_trans", 
                             "node_color_trans",
                             "edge_color_trans",
                             "tree_color_trans",
                             "node_label_size_trans", 
                             "edge_label_size_trans", 
                             "tree_label_size_trans",
                             "node_label_color_trans",
                             "edge_label_color_trans", 
                             "tree_label_color_trans")])
  args <- check_color_range(args[c("node_color_range",
                                   "edge_color_range", 
                                   "tree_color_range",
                                   "node_label_color_range", 
                                   "edge_label_color_range", 
                                   "tree_label_color_range")])
  
  if (length(overlap_avoidance) == 0 || ! is.numeric(overlap_avoidance)) {
    stop("Argument 'overlap_avoidance' must be a numeric of length 1.")
  }
  layout <- match.arg(layout, layout_functions())
  if (!is.null(initial_layout)) {
    initial_layout <- match.arg(initial_layout, layout_functions())
  }
  
}


#' Check for length or taxon id
#' 
#' Check that inputs are named by taxon ID or are the same length as the number of taxa.
#' Print a message if an argument does not have a taxon ID but is same length as the number of taxa. 
#' 
#' @param obj taxmap object
#' @param args A list of arguments
#' 
#' @return A named list of argument values, potentially modified
#' 
#' @keywords internal
check_length_or_id <- function(obj, args) {
  
  check_one <- function(name, value) {
    if (is.null(value)) {
      return(value)
    }
    if (is.null(names(value))) {
      if (length(value) == 1) {
        return(stats::setNames(rep(value, length(obj$taxon_ids())),
                               obj$taxon_ids()))
      }
      if (length(value) == length(obj$taxon_ids()))
      {
        message('Argument "', name, 
                "' has no taxon IDs. Assuming values are in the same order as taxa.")
        return(stats::setNames(value, obj$taxon_ids()))
      } else {
        stop('Argument "', name, '" has no taxon IDs and is a different length (', 
             length(value), ') than the number of taxa (',  length(obj$taxon_ids()), 
             '). Must either be named by taxon IDs or be of a length equal to the number of taxa or of length 1.', 
             call. = FALSE)
      }
    }
    is_id <- obj$is_taxon_id(names(value))
    if (! all(is_id)) {
      if (length(value) == length(obj$taxon_ids())) {
        message('Argument "', name, 
                "' is named, but not by taxon IDs. Assuming values are in the same order as taxa.",
                ' The following names are not taxon ids:\n', 
                limited_print(names(value)[! is_id], prefix = "  ", type = "silent"))
        return(stats::setNames(value, obj$taxon_ids()))
      } else {
        stop('Argument "', name, '" is not named by taxon ids. The following names are not taxon ids:\n', 
             limited_print(names(value)[! is_id], prefix = "  ", type = "silent"), call. = FALSE)
      }
    }
    return(value)
  }
  
  stats::setNames(purrr::map2(names(args), args, check_one),
                  names(args))
  
}

#' Check for NAs
#' 
#' Check for NAs in the values and IDs of parameters
#' 
#' @param args A list of arguments
#' @param warn If TRUE, warn the user about NAs in argument values
#' 
#' @return A named list of argument values, potentially modified
#' 
#' @keywords internal
check_for_na <- function(args, warn) {
  
  check_one <- function(name, value) {
    if (is.null(value)) {
      return(value)
    }
    if (warn && any(is.na(value))) {
      warning(call. = FALSE, 
              'Argument "', name, "' has NAs. This might cause odd behavior.",
              ' The following taxon ids are NA:\n', 
              limited_print(names(value)[is.na(value)], prefix = "  ", type = "silent"))
    }
    if (any(is.na(names(value)))) {
      message('Argument "', name, "' has NAs in its taxon IDs. These will not be used in the plot.",
              ' The following indexes have NA for IDs:\n', 
              limited_print(which(is.na(names(value))), prefix = "  ", type = "silent"))
      value <- value[! is.na(names(value))]
    }
    return(value)
  }
  
  stats::setNames(purrr::map2(names(args), args, check_one),
                  names(args))
}
