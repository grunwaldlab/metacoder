#' Bounding box coords for labels
#' 
#' Given a position, size, rotation, and justification of a lable, calculate the bounding box coordinants
#' 
#' @param x Horizontal position of center of text grob
#' @param y Vertical position of center of text grob
#' @param height Height of text grob
#' @param rotation Rotation in radians
#' @param just Justification. e.g. "left-top"
#' 
#' @keywords internal
label_bounds <- function(label, x, y, height, rotation, just) {
  
  process_one <- function(label, x, y, height, rotation, just) {
    # Calculate some handy values used later
    width <- height * text_grob_length(label) # The length of the text
    from_center_to_corner <- sqrt(height ^ 2 + width ^ 2) * 0.5 # The length between the center of the text box and a corner
    angle_to_corner <- atan2(height, width) # The angle between the center of the text box and a corner
    
    # Find the coordinates for the four corners, assuming a central justification
    coords <- data.frame(stringsAsFactors = FALSE,
                         x = c(- cos(rotation - angle_to_corner),  # top left
                               cos(rotation + angle_to_corner),  # top right
                               cos(rotation - angle_to_corner), # bottom right
                               - cos(rotation + angle_to_corner)),  # bottom left
                         y = c(- sin(rotation - angle_to_corner),  # top left
                               sin(rotation + angle_to_corner),  # top right
                               sin(rotation - angle_to_corner), # bottom right
                               - sin(rotation + angle_to_corner))  # bottom left
    ) * from_center_to_corner
    
    # Offset based on justification if not centered
    if        (just == "left-top") {
      coords$x <- coords$x - coords$x[1]
      coords$y <- coords$y - coords$y[1]
    } else if (just == "center-top") {
      coords$x <- coords$x - cos(rotation + pi * 0.5) * height * 0.5
      coords$y <- coords$y - sin(rotation + pi * 0.5) * height * 0.5
    } else if (just == "right-top") {
      coords$x <- coords$x - coords$x[2]
      coords$y <- coords$y - coords$y[2]
    } else if (just == "left-center" || just == "left") {
      coords$x <- coords$x + cos(rotation) * width * 0.5
      coords$y <- coords$y + sin(rotation) * width * 0.5
    } else if (just == "right-center" || just == "right") {
      coords$x <- coords$x - cos(rotation) * width * 0.5
      coords$y <- coords$y - sin(rotation) * width * 0.5
    } else if (just == "left-bottom") {
      coords$x <- coords$x - coords$x[4]
      coords$y <- coords$y - coords$y[4]
    } else if (just == "center-bottom") {
      coords$x <- coords$x + cos(rotation + pi * 0.5) * height * 0.5
      coords$y <- coords$y + sin(rotation + pi * 0.5) * height * 0.5
    } else if (just == "right-bottom") {
      coords$x <- coords$x - coords$x[3]
      coords$y <- coords$y - coords$y[3]
    }
    
    # Adjust for input coordinates
    coords$x <- coords$x + x
    coords$y <- coords$y + y
    
    # Add input label and return
    cbind(data.frame(label = label, stringsAsFactors = FALSE),  coords)
  }
  
  output <- do.call(rbind, mapply(FUN = process_one, SIMPLIFY = FALSE,
                                  label, x, y, height, rotation, just))
  rownames(output) <- NULL
  return(output)
}


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
#' @param inverse (\code{logical} of length 1) If \code{TRUE}, return the inverse of the selected function.
#' 
#' @keywords internal
transform_data <- function(func = NULL, data = NULL, inverse = FALSE) {
  sign <- function(x) {
    ifelse(x < 0, -1, 1)
  }
  
  funcs <- list("linear" = function(x) {x},
                "area" = function(x) {sign(x) * (abs(x)/pi)^(1/2)},
                "log10" = function(x) {log(x, base = 10)},
                "log2" = function(x) {log(x, base = 2)},
                "ln" = function(x) {log(x)},
                "log10 area" = function(x) {log((x/pi)^(1/2), base = 10)},
                "log2 area" = function(x) {log((x/pi)^(1/2), base = 2)},
                "ln area" =  function(x) {log((x/pi)^(1/2))})
  
  inverse_funcs <- list("linear" = function(x) {x},
                        "area" = function(x) {sign(x) * pi * (x ^ 2)},
                        "log10" = function(x) {10 ^ x},
                        "log2" = function(x) {2 ^ x},
                        "ln" = function(x) {exp(1) ^ x},
                        "log10 area" = function(x) {pi * (abs(10 ^ x) ^ 2)},
                        "log2 area" = function(x) {pi * (abs(2 ^ x) ^ 2)},
                        "ln area" =  function(x) {pi * (abs(exp(1) ^ x) ^ 2)})
  
  
  if (is.null(data) & is.null(func)) {
    return(names(funcs))
  }
  
  if (inverse) {
    returned_funcs <- inverse_funcs
  } else {
    returned_funcs <- funcs
  }
  
  if (is.null(data)) {
    return(returned_funcs[[func]])
  } else if (is.numeric(data)) {
    return(vapply(X = data, FUN = returned_funcs[[func]], FUN.VALUE = numeric(1)))
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

