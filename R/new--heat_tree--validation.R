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
  to_check <- c("node_size", 
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
                "tree_label")
  args[to_check] <- check_length_or_id(obj, args[to_check])
  
  # Look for NAs in values and taxon names
  to_check <- c("node_size",
                "edge_size",
                "node_label_size",
                "edge_label_size",
                "tree_label_size",
                "node_color", 
                "edge_color", 
                "tree_color",
                "node_label_color",
                "edge_label_color", 
                "tree_label_color")
  args[to_check] <- check_for_na(args[to_check], warn = TRUE)
  
  to_check <- c("node_label",
                "edge_label",
                "tree_label")
  args[to_check] <- check_for_na(args[to_check], warn = FALSE)
  
  # Check values of arguments specifying node sizes
  to_check <- c("node_size",
                "edge_size",
                "node_label_size", 
                "edge_label_size", 
                "tree_label_size")
  args[to_check] <- check_size(args[to_check])
  
  # Check values of arguments specifying size ranges  
  to_check <- c("node_size_range",  
                "edge_size_range", 
                "node_label_size_range",
                "edge_label_size_range", 
                "tree_label_size_range",
                "node_size_interval", 
                "edge_size_interval")
  args[to_check] <- check_size_range(args[to_check])
  
  # Check values of arguments specifying transformations
  to_check <- c("node_size_trans", 
                "edge_size_trans", 
                "node_color_trans",
                "edge_color_trans",
                "tree_color_trans",
                "node_label_size_trans", 
                "edge_label_size_trans", 
                "tree_label_size_trans",
                "node_label_color_trans",
                "edge_label_color_trans", 
                "tree_label_color_trans")
  args[to_check] <- check_trans(args[to_check])
  
  # Check values of arguments specifying color ranges  
  to_check <- c("node_color_range",
                "edge_color_range", 
                "tree_color_range",
                "node_label_color_range", 
                "edge_label_color_range", 
                "tree_label_color_range")
  args[to_check] <- check_color_range(args[to_check])
  
  # check overlap avoidance argument
  if (! is.null(args[["overlap_avoidance"]])) {
    warning(call. = FALSE, 'The "overlap_avoidance" argument is depreciated and will be ignored. Use the "*_size_range" arguments to manually specify a size range if the default does not work well.')
  }
  
  # Check layout arguments
  if (! args[["layout"]] %in% layout_functions()) {
    stop(call. = FALSE, 'The layout argument must be one of the following:\n',
         limited_print(layout_functions(), prefix = "  ", type = "silent"))
  }
  if (! is.null(args[["initial_layout"]]) && ! args[["initial_layout"]] %in% layout_functions()) {
    stop(call. = FALSE, 'The initial_layout argument must be one of the following:\n',
         limited_print(layout_functions(), prefix = "  ", type = "silent"))
  }
  if (args[["initial_layout"]] == args[["layout"]]) {
    message('The "layout" and "initial_layout" arguments are the same, so the "initial_layout" will be ignored.')
    args["initial_layout"] <- NULL
  }
  
  return(args)
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


#' Check size arguments
#' 
#' Check that arguments specifying size make sense
#' 
#' @param args A list of arguments
#' 
#' @return A named list of argument values, potentially modified
#' 
#' @keywords internal
check_size <- function(args) {
  
  check_one <- function(name, value) {
    if (is.null(value)) {
      return(value)
    }
    if (any(!is.na(value) & is.na(as.numeric(value)))) {
      stop("Argument '", name, "' is not numeric.", call. = FALSE)
    }
    if (any(is.infinite(value))) {
      is_real <- ! is.infinite(value) & ! is.na(value)
      if (any(is_real))
      {
        message('Infinite values found for "', name,
                '". These will be graphed in the same way as the largest (Inf) or smallest (-Inf) real number supplied.')
        value[is.infinite(value) & value < 0] <- min(value[is_real], na.rm = TRUE)
        value[is.infinite(value) & value > 0] <- max(value[is_real], na.rm = TRUE)
      } else {
        stop('Argument "', name, '" has no finite, non-NA values.', call. = FALSE)
      }
    }
    return(value)
  }
  
  
  stats::setNames(purrr::map2(names(args), args, check_one),
                  names(args))
}


#' Verify size range parameters
#' 
#' Verify size range parameters
#' 
#' @param args A list of arguments
#' 
#' @return A named list of argument values, potentially modified
#' 
#' @keywords internal
check_size_range <- function(args) {
  
  check_one <- function(name, value) {
    if (is.null(value)) {
      return(value)
    }
    if (length(value) != 2) {
      stop(call. = FALSE, 'Size range argument "', name, '" must be of length 2.')
    }
    if (all(!is.na(value)) && value[2] < value[1]) {
      stop(call. = FALSE, 'The min value of size range argument "', name, '" is greater than its max.')
    }
    return(value)
  }
  
  stats::setNames(purrr::map2(names(args), args, check_one),
                  names(args))
}


#' Verify transformation parameters
#' 
#' Verify transformation parameters
#' 
#' @param args A list of arguments
#' 
#' @return A named list of argument values, potentially modified
#' 
#' @keywords internal
check_trans <- function(args) {
  
  check_one <- function(name, value) {
    if (is.null(value)) {
      return(value)
    }
    if (! is.function(value) && ! value %in% transform_data()) {
      stop(call. = FALSE, 'Transformation argument "', name,
           '" must be a function or the name of a built-in transformation function:\n',
           limited_print(transform_data(), prefix = "  ", type = "silent"))
    }
    return(value)
  }
  
  stats::setNames(purrr::map2(names(args), args, check_one),
                  names(args))
}


#' Verify color range parameters
#' 
#' Verify color range parameters
#' 
#' @param args A list of arguments
#' 
#' @return A named list of argument values, potentially modified
#' 
#' @keywords internal
check_color_range <- function(args) {
  
  check_one <- function(name, value) {
    if (is.null(value)) {
      return(value)
    }
    if (length(value) == 0) {
      stop(call. = FALSE, 'Color range argument "', name, '" has no values.')
    }
    if (any(! grepl("^#[0-9a-fA-F]{3,8}$", value) & ! value %in% grDevices::colors())) {
      stop(call. = FALSE, 'Color range argument "', name, '" must be hex color codes or a name returned by "colors()"')
    }
    return(value)
  }
  
  stats::setNames(purrr::map2(names(args), args, check_one),
                  names(args))
}
