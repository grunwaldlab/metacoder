#' Validate arguments given to heat_tree
#'
#' Validate and standardize arguments given to \code{\link{heat_tree}} or
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
  # Check that the main input is a taxonomy with at least one taxon
  if (! "Taxonomy" %in% class(obj)) {
    stop(call. = FALSE, '"obj" is not a Taxonomy or Taxmap object: ', class(obj)[1])
  }
  if (length(obj$taxon_ids()) <= 0) {
    stop(call. = FALSE, '"obj" contains no taxa to plot.')
  }
  
  # Check that a taxmap/taxonomy object was not accidentally passed to another parameter
  is_taxonomy <- purrr::map_lgl(args, function(x) "Taxonomy" %in% class(x))
  if (any(is_taxonomy)) {
    stop(call. = FALSE, 'Argument "', names(args)[is_taxonomy][1], '" was given a "', class(args[is_taxonomy][[1]])[1], 
         '" object. This usually happens when you try to pass an object with "%>%" and also supply it in the function call.')
  }
  
  # Check that inputs are named by taxon ID or are the same length as the number of taxa
  to_check <- c("node_size", 
                "edge_size",
                "node_label_size", 
                "edge_label_size",  
                "tree_label_size",
                "node_color", 
                "edge_color",
                "node_label_color", 
                "edge_label_color", 
                "tree_label_color",
                "node_label", 
                "edge_label", 
                "tree_label",
                "node_shape",
                "edge_shape")
  args[to_check] <- check_aes_input(obj, args[to_check])
  
  # Look for NAs in values and taxon names
  to_check <- c("node_size",
                "edge_size",
                "node_label_size",
                "edge_label_size",
                "tree_label_size",
                "node_color", 
                "edge_color", 
                "node_label_color",
                "edge_label_color", 
                "tree_label_color",
                "node_shape",
                "edge_shape")
  args[to_check] <- check_for_na(args[to_check], warn = TRUE)
  
  to_check <- c("node_label",
                "edge_label",
                "tree_label")
  args[to_check] <- check_for_na(args[to_check], warn = FALSE)
  
  # Check values of arguments specifying sizes
  to_check <- c("node_size",
                "edge_size",
                "node_label_size", 
                "edge_label_size", 
                "tree_label_size")
  args[to_check] <- check_size(args[to_check])
  
  # Check values of arguments specifying colors
  to_check <- c("node_color",
                "edge_color",
                "node_label_color", 
                "edge_label_color", 
                "tree_label_color")
  args[to_check] <- check_color(args[to_check])
  
  # Check value of arguments specifying shape ranges
  args["node_shape_range"] <- check_shape_range(args["node_shape_range"], node_shape_functions())
  args["edge_shape_range"] <- check_shape_range(args["edge_shape_range"], edge_shape_functions())
  
  # Check values of arguments specifying shapes
  args["node_shape"] <- check_shape(args["node_shape"], args, names(node_shape_functions()))
  args["edge_shape"] <- check_shape(args["edge_shape"], args, names(edge_shape_functions()))
  
  # Check value of arguments for labels 
  to_check <- c("node_label",
                "edge_label",
                "tree_label", 
                "node_color_axis_label", 
                "node_size_axis_label",
                "edge_color_axis_label",
                "edge_size_axis_label",
                "title")
  args[to_check] <- check_labels(args[to_check])
  
  # Check values of arguments specifying size ranges  
  to_check <- c("node_size_range",  
                "edge_size_range", 
                "node_label_size_range",
                "edge_label_size_range", 
                "tree_label_size_range")
  args[to_check] <- check_size_range(args[to_check])
  
  # Check values of arguments specifying intervals
  to_check <- c("node_color_interval", 
                "edge_color_interval",
                "node_size_interval", 
                "edge_size_interval")
  args[to_check] <- check_intervals(args[to_check], args)
  
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
  args[to_check] <- check_color_range(args[to_check], args)
  
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
  if (! is.null(args[["initial_layout"]]) && args[["initial_layout"]] == args[["layout"]]) {
    message('The "layout" and "initial_layout" arguments are the same, so the "initial_layout" will be ignored.')
    args["initial_layout"] <- NULL
  }
  
  return(args)
}


#' Check for length or taxon id
#' 
#' Check that inputs are named by taxon ID or are the same length as the number of taxa.
#' Print a message if an argument does not have a taxon ID but is same length as the number of taxa. 
#' Converts all values to lists.
#' 
#' @param obj taxmap object
#' @param args A list of arguments
#' 
#' @return A named list of argument values, potentially modified
#' 
#' @keywords internal
check_aes_input <- function(obj, args) {
  
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
  
  output <- stats::setNames(purrr::map2(names(args), args, check_one),
                            names(args))
  
  # convert to lists
  output <- stats::setNames(purrr::map(output, function(x) {
    if (is.null(x)) {
      return(NULL)
    } else {
      if (is.factor(x)) {
        return(stats::setNames(as.list(unname(x)), names(x)))
      } else {
        return(as.list(x))
      }
    }
  }),
  names(args))
  return(output)
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
#' Check that arguments specifying size make sense. Parmeters specifying size must be numeric and
#' not NA and at least one must be finite. Infinite values will be converted to the max/min of finite values.
#'
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
    is_infinite <- purrr::map_lgl(value, is.infinite)
    if (any(is_infinite)) {
      is_real <- ! is_infinite & ! is.na(value)
      if (any(is_real))
      {
        message('Infinite values found for "', name,
                '". These will be graphed in the same way as the largest (Inf) or smallest (-Inf) real number supplied.')
        value[is_infinite & value < 0] <- min(value[is_real], na.rm = TRUE)
        value[is_infinite & value > 0] <- max(value[is_real], na.rm = TRUE)
      } else {
        stop('Argument "', name, '" has no finite, non-NA values.', call. = FALSE)
      }
    }
    return(as.list(value))
  }
  
  
  stats::setNames(purrr::map2(names(args), args, check_one),
                  names(args))
}


#' Check shape arguments
#'
#' Check that arguments specifying shapes make sense.
#'
#' @param shape_args A list of arguments specifying shapes
#' @param args All args
#' @param default_shapes A character vector of built-in function names for shapes
#'
#' @return A named list of argument values, potentially modified
#'
#' @keywords internal
check_shape <- function(shape_args, args, default_shapes) {
  
  check_one <- function(name, value) {
    if (is.null(value) || length(value) == 0) {
      stop('Argument "', name, '" has no values.', call. = FALSE)
    }
    
    # Check that character values are specified in the range functions or are a built in function
    range_param <- args[[paste0(name, "_range")]]
    char_values <- get_characters(value)
    valid_values <- c(default_shapes, names(range_param))
    invalid_values <- char_values[! char_values %in% names(range_param)]
    if (length(invalid_values) > 0) {
      stop('Shape parameter "', names, '" has ', length(invalid_values), 
           ' invalid character value', ifelse(length(invalid_values) > 1, 's', ''), ':\n',
           limited_print(unique(invalid_values), prefix = '  ', type = 'silent'), '\n',
           'Valid choices include:\n',
           limited_print(valid_values, prefix = '  ', type = 'silent'))
    }
    
    return(as.list(value))
  }
  
  
  stats::setNames(purrr::map2(names(shape_args), shape_args, check_one),
                  names(shape_args))
}


#' Check color arguments
#'
#' Check that arguments specifying color make sense. They must be either numbers or a character
#' representing a hex color code and a color name in colors()
#'
#' @param args A list of arguments
#'
#' @return A named list of argument values, potentially modified
#'
#' @keywords internal
check_color <- function(args) {
  
  check_one <- function(name, value) {
    if (is.null(value)) {
      return(value)
    }
    is_char <- purrr::map_lgl(value, is.character)
    is_num <- purrr::map_lgl(value, is.numeric)
    if (any(is_num) && (any(is_char & ! is_color_char(value)))) {
      stop(call. = FALSE,
           'Values given to color argument "', name, 
           '" have numbers and non-color character values. Color values must be either numbers or categories, perhaps mixed with explicit colors, but not both numbers and categories.')
    }
    return(as.list(value))
  }
  
  stats::setNames(purrr::map2(names(args), args, check_one),
                  names(args))
}


#' Verify size range parameters
#' 
#' Verify size range parameters. 
#' 
#' @param args A list of arguments
#' 
#' @return A named list of argument values, potentially modified
#' 
#' @keywords internal
check_size_range <- function(args) {
  
  check_one <- function(name, value) {
    if (is.null(value)) {
      return(c(NA_real_, NA_real_))
    }
    if (length(value) != 2) {
      stop(call. = FALSE, 'Size range argument "', name, '" must be of length 2.')
    }
    if (all(!is.na(value)) && value[2] < value[1]) {
      stop(call. = FALSE, 'The min value of size range argument "', name, '" is greater than its max.')
    }
    if (any(value > 1) || any(value < 0)) {
      stop(call. = FALSE, 'A value of size range argument "', name, '" is not between 0 and 1. A value larger than 1 means that it should be bigger than the plotted area and less than 0 does not make sense.')
    }
    if (all(is.na(value))) {
      return(c(NA_real_, NA_real_))
    }
    return(value)
  }
  
  stats::setNames(purrr::map2(names(args), args, check_one),
                  names(args))
}


#' Verify interval parameters
#' 
#' Verify interval parameters
#' 
#' @param intervals A list of arguments specifying intervals
#' @param args A list of all the args
#' 
#' @return A named list of argument values, potentially modified
#' 
#' @keywords internal
check_intervals <- function(intervals, args) {
  
  check_one <- function(name, value) {
    aes_values <- get_numerics(args[[sub(name, pattern = "_interval$", replacement = "")]])
    if (is.null(value)) {
      if (length(aes_values) > 0) {
        return(range(aes_values, na.rm = TRUE))
      } else {
        return(NULL)
      }
    }
    if (length(value) != 2) {
      stop(call. = FALSE, 'Interval argument "', name, '" must be of length 2.')
    }
    if (all(!is.na(value)) && value[2] < value[1]) {
      stop(call. = FALSE, 'The min value of interval argument "', name, '" is greater than its max.')
    }
    if (length(aes_values) > 0) {
      if (is.na(value[1])) {
        value[1] <- min(aes_values, na.rm = TRUE)
      }
      if (is.na(value[2])) {
        value[2] <- max(aes_values, na.rm = TRUE)
      }
    }
    return(value)
  }
  
  stats::setNames(purrr::map2(names(intervals), intervals, check_one),
                  names(intervals))
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
    if (! value %in% names(heat_tree_transform_funcs())) {
      stop(call. = FALSE, 'Transformation argument "', name,
           '" must be the name of a built-in transformation function:\n',
           limited_print(transform_data(), prefix = "  ", type = "silent"))
    }
    return(heat_tree_transform_funcs()[[value]])
  }
  
  stats::setNames(purrr::map2(names(args), args, check_one),
                  names(args))
}


#' Verify color range parameters
#' 
#' Verify color range parameters
#' 
#' @param ranges A list of range arguments to check
#' @param args A list of all arguments
#' 
#' @return A named list of argument values, potentially modified
#' 
#' @keywords internal
check_color_range <- function(ranges, args) {
  
  # Identify color range arguments
  aes_args <- args[sub(names(ranges), pattern = "_range$", replacement = "")]
  names(aes_args) <- names(ranges)
  
  check_one <- function(name, value) {
    if (is.null(value)) {
      if (is_categorical(aes_args[[name]])) {
        return(qualitative_palette())
      } else {
        return(quantative_palette())
      }
    }
    if (length(value) == 0) {
      stop(call. = FALSE, 'Color range argument "', name, '" has no values.')
    }
    if (any(! is_color_char(value))) {
      stop(call. = FALSE, 'Color range argument "', name, '" must be hex color codes or a name returned by "colors()"')
    }
    return(value)
  }
  
  stats::setNames(purrr::map2(names(ranges), ranges, check_one),
                  names(ranges))
}


#' Verify label parameters
#' 
#' Verify label parameters
#' 
#' @param args A list of arguments
#' 
#' @return A named list of argument values, potentially modified
#' 
#' @keywords internal
check_labels <- function(args) {
  
  check_one <- function(name, value) {
    if (is.null(value)) {
      return(value)
    }
    if (length(value) == 0) {
      stop(call. = FALSE, 'Label argument "', name, '" has no values.')
    }
    value <- stats::setNames(as.character(value),
                             names(value))
    return(value)
  }
  
  stats::setNames(purrr::map2(names(args), args, check_one),
                  names(args))
}


#' Verify shape range parameters
#' 
#' Verify shape range parameters. These should be either character or named functions in a list.
#' 
#' @param args A list of arguments
#' 
#' @return A named list of argument values, potentially modified
#' 
#' @keywords internal
check_shape_range <- function(args, valid_shapes) {
  
  check_one <- function(name, value) {
    
    # Check that there is at least one value
    if (is.null(value)|| length(value) == 0) {
      stop(call. = FALSE, 'Shape range argument "', name, '" has no values.')
    }
    
    # Check that all inputs are either characters or named functions
    is_valid_format <- purrr::map_lgl(seq_len(length(value)), function(i) {
      is.character(value[[i]]) || (is.function(value[[i]]) && ! is.null(names(value)[[i]]))
    })
    if (any(! is_valid_format)) {
      stop('Shape range parameter "', names, '" has ', sum(! is_valid_format),' invalid values.',
           'Values must be either named functions or the name of a built-in function.')
           
    }
    
    # Check that character values are the name of a built-in function
    char_values <- get_characters(value)
    invalid_values <- char_values[! char_values %in% names(valid_shapes)]
    if (length(invalid_values) > 0) {
      stop('Shape range parameter "', names, '" has ', length(invalid_values), 
           ' invalid character value', ifelse(length(invalid_values) > 1, 's', ''), ':\n',
           limited_print(unique(invalid_values), prefix = '  ', type = 'silent'), '\n',
           'Valid choices include:\n',
           limited_print(names(valid_shapes), prefix = '  ', type = 'silent'))
    }
    
    # Convert characeter values to the functions they specify
    to_replace <- purrr::map_lgl(value, is.character)
    names(value)[to_replace] <- value[to_replace]
    value[to_replace] <- valid_shapes[unlist(value[to_replace])]
    
    # Convert to list of functions
    return(as.list(value))
  }
  
  stats::setNames(purrr::map2(names(args), args, check_one),
                  names(args))
  
}