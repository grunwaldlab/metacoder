#' Verify size range parameters
#' 
#' Verify size range parameters
#' 
#' @param args (\code{character}) The names of arguments to verify.
#' 
#' @keywords internal
verify_size_range <- function(args) {
  for (arg in args) {
    value <- get(arg, pos = parent.frame())
    if (length(value) != 2) {
      stop(paste0("Argument ", arg, " must be of length 2."))
    }
    if (all(!is.na(value)) && value[2] < value[1]) {
      stop(paste0("The min value of ", arg, " is greater than its max."))
    }
  }
}


#' Verify transformation function parameters
#' 
#' Verify transformation function parameters
#' 
#' @param args (\code{character}) The names of arguments to verify.
#' 
#' @keywords internal
verify_trans <- function(args) {
  for (arg in args) {
    value <- get(arg, pos = parent.frame())
    if (! is.function(value) && ! value %in% transform_data()) {
      stop(paste0("Argument '", arg,
                  "' must be a function or the name of a built-in transformation function."))
    }
  }
}

#' Verify size parameters
#' 
#' Verify size parameters
#' 
#' @param args (\code{character}) The names of arguments to verify.
#' 
#' @keywords internal
verify_size <- function(args) {
  for (arg in args) {
    value <- get(arg, pos = parent.frame())
    if (any(!is.na(value) & is.na(as.numeric(value)))) {
      stop(paste0("Argument '", arg, "' is not numeric."))
    }
  }
}


#' Verify color range parameters
#' 
#' Verify color range parameters
#' 
#' @param args (\code{character}) The names of arguments to verify.
#' 
#' @keywords internal
verify_color_range <- function(args) {
  for (arg in args) {
    value <- get(arg, pos = parent.frame())
    if (any(! grepl("^#[0-9a-fA-F]{3,8}$", value) & ! value %in% grDevices::colors())) {
      stop(paste0("Argument '", arg, "' must be hex color codes or a name returned by 'colors()'"))
    }
  }
}



#' Verify label count
#' 
#' Verify label count
#' 
#' @param args (\code{character}) The names of arguments to verify.
#' 
#' @keywords internal
verify_label_count <- function(args) {
  for (arg in args) {
    value <- get(arg, pos = parent.frame())
    if (value %% 1 != 0 || value < 0) {
      stop(paste0("Argument '", arg, "' must be a positive integer."))
    }
  }
}


#' Check length of graph attributes
#' 
#' Length should divind evenly into the number of taxon/parent IDs
#' 
#' @keywords internal
check_element_length <- function(args) {
  for (arg in args) {
    observed_length <- length(get(arg, pos = parent.frame()))
    correct_length <- length(get("taxon_id", pos = parent.frame()))
    if (observed_length < 1) {
      stop(paste0("Argument '", arg, "' is empty."))
    }
    if (correct_length %% observed_length != 0) {
      stop(paste0("Length of argument'", arg, "' must be a factor of the length of 'taxon_id'"))
    }
  }
}


#' Transformation functions
#' 
#' Functions used by plotting funtions to transform data.
#' Calling the function with no parameters returns available function names.
#' Calling with just the function name returns the transformation function
#' 
#' @param data (\code{numeric}) Data to transform
#' @param func (\code{character}) Name of transformation to apply.
#' 
#' @keywords internal
transform_data <- function(func = NULL, data = NULL) {
  sign <- function(x) {
    ifelse(x < 0, -1, 1)
  }
  
  funcs <- list("radius" = function(x) {x},
                "area" = function(x) {sign(x) * (abs(x)/pi)^(1/2)},
                "log10 radius" = function(x) {log(x, base = 10)},
                "log2 radius" = function(x) {log(x, base = 2)},
                "ln radius" = function(x) {log(x)},
                "log10 area" = function(x) {log((x/pi)^(1/2), base = 10)},
                "log2 area" = function(x) {log((x/pi)^(1/2), base = 2)},
                "ln area" =  function(x) {log((x/pi)^(1/2))})
  
  if (is.null(data) & is.null(func)) {
    return(names(funcs))
  } else if (is.null(data)) {
    return(funcs[[func]])
  } else if (is.numeric(data)) {
    return(vapply(X = data, FUN = funcs[[func]], FUN.VALUE = numeric(1)))
  } else {
    return(data)
  }
}


#' Pick labels to show
#' 
#' Pick labels to show based off a column name to sort by and a maximum number
#' 
#' @param my_data \code{data.frame}
#' @param label_max \code{numeric} of length 1
#' @param sort_by_column \code{character} of length 1; the name of a column in \code{my_data}
#' @param label_column \code{character} of length 1; the name of a column in \code{my_data}
#' containing labels
#' 
#' @return \code{character} IDs of rows with labels to show
#' 
#' @keywords internal
select_labels <- function(my_data, label_max, sort_by_column, label_column) {
  if (is.null(label_max) || is.na(label_max) || nrow(my_data) <= label_max) {
    labels_shown <- rownames(my_data)
  } else {
    index_order <- rev(do.call(order, my_data[ , sort_by_column]))
    top_indexes <- index_order[1:label_max]
    labels_shown <- rownames(my_data)[top_indexes]
  }
  labels_shown <- labels_shown[!is.na(my_data[labels_shown, label_column])]  # Do not make grobs for NA
  return(rownames(my_data) %in% labels_shown)
}

