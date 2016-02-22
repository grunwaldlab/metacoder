#| ## Tree plotting
#|
#|
#|
#|
#|
#|
#|
#|
#|
#|
#|
#|
#|
#|
#|
#|
#===================================================================================================
#' Plot items on taxonomy
#' 
#' Plots the distribution of values associated with an taxonomic classification.
#' Colors and sizes of vertexes and edges can be mapped to variables.
#' Uses \code{igraph} to make layout and \code{ggplot2} to make plots.
#' 
#' @param taxon_id (\code{character}} The unique ids of the taxon for each row.
#' @param parent_id (\code{character}} The unique id of supertaxon \code{taxon_id} is a part of.
#' @param vertex_size (\code{numeric)} The value to base vertex size on. Default: inverse of depth 
#' of taxa in hierarchy. 
#' @param vertex_size_range (\code{numeric} of length 2) The minimum and maximum size of vertexes.
#' If either value is \code{NA}, the missing value(s) will be optimized to maximize range and
#' minimizing overlaps. Defualt: \code{c(NA, NA)}.
#' @param vertex_size_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{vertex_size}. Alternativly one of the following \code{character} values. Default:
#' \code{"area"}.
#' \describe{
#'   \item{"radius"}{Proprotional to radius/diameter of vertex}
#'   \item{"area"}{circular area; better perceptual accuracy than \code{"radius"}}
#'   \item{"log10 radius"}{Log base 10 of radius}
#'   \item{"log2 radius"}{Log base 2 of radius}
#'   \item{"ln radius"}{Log base e of radius}
#'   \item{"log10 area"}{Log base 10 of area}
#'   \item{"log2 area"}{Log base 2 of area}
#'   \item{"ln area"}{Log base e of area}
#' }
#' @param vertex_color (\code{numeric} OR \code{character}) The value to base vertex color on. 
#' If a numeric vector is given, it is used to  construct a color scale. Hex values or color 
#' names can be used (e.g. \code{#000000} or \code{"black"}). Default: grey.
#' @param vertex_color_range (valid argument of \code{col2rgb}) A series of colors corresponding 
#' to \code{vertex_color}. Default: Color-blind friendly palette. 
#' @param vertex_color_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{vertex_size}. Alternativly, one of the \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: \code{"area"}.
#' @param edge_size (\code{numeric)} The value to base edge size on. Default: relative to vertex
#' size. 
#' @param edge_size_range (\code{numeric} of length 2) The minimum and maximum size of edges.
#' If either value is \code{NA}, the missing value(s) will be optimized relative to vertex size. 
#' Default: relative to vertex size. 
#' @param edge_size_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{edge_size}. Alternativly one of the following \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: same as \code{vertex_size_trans}. 
#' @param edge_color (\code{numeric} OR \code{character}) The value to base edge color on. 
#' If a numeric vector is given, it is used to  construct a color scale. Hex values or color 
#' names can be used (e.g. \code{#000000} or \code{"black"}). Default: same as vertex color.
#' @param edge_color_range (valid argument of \code{col2rgb}) A series of colors corresponding 
#' to low-high statistics supplied to \code{edge_color}. Default: same as vertex color.
#' @param edge_color_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{edge_size}. Alternativly one of the following \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: same as vertex color transformation.
#' @param vertex_label The values of labels over vertcies. Use \code{NA} to exclude labels.
#' Default: no labels.
#' @param vertex_label_size (\code{numeric)} The value to base veterx label size on. Default: 
#' relative to veterx size.
#' @param vertex_label_size_range (\code{numeric} of length 2) The minimum and maximum size of labels.
#' If either value is \code{NA}, the missing value(s) will be optimized relative to vertex size. 
#' Default: relative to vertex size. 
#' @param vertex_label_size_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{vertex_label_size}. Alternativly one of the following \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: same as \code{vertex_size_trans}.
#' @param vertex_label_color (\code{numeric} OR \code{character}) The value to base label color on. 
#' If a numeric vector is given, it is used to  construct a color scale. Hex values or color 
#' names can be used (e.g. \code{#000000} or \code{"black"}). Default: black.
#' @param vertex_label_color_range (valid argument of \code{col2rgb}) A series of colors corresponding 
#' to low-high statistics supplied to \code{vertex_label_color}. Default: Color-blind friendly palette. 
#' @param vertex_label_color_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{vertex_label_size}. Alternativly one of the \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: \code{"area"}.
#' @param vertex_label_max (\code{numeric}) The maximum number of vertex labels. Default: 20.
#' @param edge_label The values of labels over edges. Default: no labels.
#' @param edge_label_size (\code{numeric)} The value to base edge label size on. Default: relative to
#' edge size.
#' @param edge_label_size_range (\code{numeric} of length 2) The minimum and maximum size of labels.
#' If either value is \code{NA}, the missing value(s) will be relative to edge size. Default: relative to
#' edge size.
#' @param edge_label_size_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{edge_label_size}. Alternativly one of the \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: same as \code{edge_size_trans}.
#' @param edge_label_color (\code{numeric} OR \code{character}) The value to base label color on. 
#' If a numeric vector is given, it is used to  construct a color scale. Hex values or color 
#' names can be used (e.g. \code{#000000} or \code{"black"}). Default: black.
#' @param edge_label_color_range (valid argument of \code{col2rgb}) A series of colors corresponding 
#' to low-high statistics supplied to \code{edge_label_color}. Default: Color-blind friendly palette.
#' @param edge_label_color_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{edge_label_size}. Alternativly one of the following \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: \code{"area"}.
#' @param edge_label_max (\code{numeric}) The maximum number of edge labels. Default: 20.
#' @param overlap_bias (\code{numeric}) The relative importance of avoiding overlaps vs maximizing 
#' size range. Default: \code{1}.
#' @param margin_size (\code{numeric} of length 2) The horizontal and vertical margins.
#' Default: \code{0, 0}.
#' @param aspect_ratio (\code{numeric}) The height / width of the plot. Default: Whatever the layout
#' function produces.
#' @param layout (\code{character}) The layout function to use. 
#' 
#' @export
new_plot_taxonomy <- function(taxon_id, parent_id, 
                              vertex_size = 1,
                              vertex_size_range = c(NA, NA),
                              vertex_size_trans = "area",
                              vertex_color = "#999999",
                              vertex_color_range = quantative_palette(),
                              vertex_color_trans = vertex_size_trans,
                              edge_size = vertex_size,
                              egde_size_range = c(NA, NA),
                              edge_size_trans = vertex_size_trans,
                              edge_color = vertex_color,
                              egde_color_range = vertex_color_range,
                              edge_color_trans = vertex_color_trans,
                              vertex_label = NA,
                              vertex_label_size = vertex_size,
                              vertex_label_size_range = c(NA, NA),
                              vertex_label_size_trans = vertex_size_trans,
                              vertex_label_color = "#000000",
                              vertex_label_color_range = quantative_palette(),
                              vertex_label_color_trans = "area",
                              edge_label = NA,
                              edge_label_size = edge_size,
                              edge_label_size_range = c(NA, NA),
                              edge_label_size_trans = edge_size_trans,
                              edge_label_color = "#000000",
                              edge_label_color_range = quantative_palette(),
                              edge_label_color_trans = "area",
                              vertex_label_max = 20,
                              edge_label_max = 20, 
                              overlap_bias = 1,
                              margin_size = c(0, 0),
                              aspect_ratio = NULL,
                              layout = igraph::layout.reingold.tilford) {
  #| ### Verify arguments =========================================================================
  if (length(taxon_id) != length(parent_id)) {
    stop("'taxon_id' and 'parent_id' must be of equal length.")
  }
  if (length(taxon_id) == 0) {
    stop("'taxon_id' and 'parent_id' are empty.")
  }
  check_element_length(c("vertex_size", "edge_size", "vertex_label_size", "edge_label_size",
                         "vertex_color", "edge_color", "vertex_label_color", "edge_label_color",
                         "vertex_label", "edge_label"))
  verify_size(c("vertex_size", "edge_size", "vertex_label_size", "edge_label_size"))
  verify_size_range(c("vertex_size_range", "egde_size_range",
                      "vertex_label_size_range", "edge_label_size_range"))
  verify_trans(c("vertex_size_trans", "vertex_color_trans", 
                 "edge_size_trans", "edge_color_trans",
                 "vertex_label_size_trans", "vertex_label_color_trans",
                 "edge_label_size_trans", "edge_label_color_trans"))
  verify_color_range(c("vertex_color_range", "egde_color_range",
                       "vertex_label_color_range", "edge_label_color_range"))
  verify_label_count(c("vertex_label_max", "edge_label_max"))
  if (length(overlap_bias) == 0 || ! is.numeric(overlap_bias)) {
    stop("Argument 'overlap_bias' must be a numeric of length 1.")
  }
  if (length(margin_size) != 2 || ! is.numeric(margin_size)) {
    stop("Argument 'margin_size' must be a numeric of length 2.")
  }
  if (! is.null(aspect_ratio) && ! is.numeric(aspect_ratio)) {
    stop("Argument 'aspect_ratio' must be a numeric of length 1.")
  }
  # TODO: verify layout
  #| ### Standardize source data ==================================================================
  data <- data.frame(stringsAsFactors = FALSE,
                     tid = as.character(taxon_id),
                     pid = as.character(parent_id),
                     vs = as.numeric(vertex_size),
                     vc = as.character(vertex_color),
                     es = as.numeric(edge_size),
                     ec = as.character(edge_color),
                     vls = as.numeric(vertex_label_size),
                     vlc = as.character(vertex_label_color),
                     els = as.numeric(edge_label_size),
                     elc = as.character(edge_label_color),
                     vl = as.character(vertex_label),
                     el = as.character(edge_label))
  
  
  #| ### Core plot data ===========================================================================
  
  
  
  #| ### Secondary plot data ======================================================================
  
  
  
  #| ### Draw plot ================================================================================
  
}



#' Verify size range parameters
#' 
#' Verify size range parameters
#' 
#' @param args (\code{character}) The names of arguments to verify.
verify_size_range <- function(args) {
  for (arg in args) {
    value <- get(arg, pos = parent.frame())
    if (length(value) != 2) {
      stop(paste0("Argument ", arg, " must be of length 2."))
    }
    if (! all(is.na(value)) && value[2] < value[1]) {
      stop(paste0("The min value of ", arg, " is greater than its max."))
    }
  }
}


#' Verify transformation function parameters
#' 
#' Verify transformation function parameters
#' 
#' @param args (\code{character}) The names of arguments to verify.
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
verify_size <- function(args) {
  for (arg in args) {
    value <- get(arg, pos = parent.frame())
    if (! is.numeric(value)) {
      stop(paste0("Argument '", arg, "' is not numeric."))
    }
  }
}


#' Verify color range parameters
#' 
#' Verify color range parameters
#' 
#' @param args (\code{character}) The names of arguments to verify.
verify_color_range <- function(args) {
  for (arg in args) {
    value <- get(arg, pos = parent.frame())
    if (! all(grepl("^#(?:[0-9a-fA-F]{3}){1,2}$", value) | value %in% colors())) {
      stop(paste0("Argument '", arg, "' must be hex color codes or a name returned by 'colors()'"))
    }
  }
}



#' Verify label count
#' 
#' Verify label count
#' 
#' @param args (\code{character}) The names of arguments to verify.
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

# 
# #' Standardize color input 
# #' 
# #' This forces user input to be the correct type.
# standardize_color <- function(data) {
#   if (all(grepl("^#(?:[0-9a-fA-F]{3}){1,2}$", data) | data %in% colors())) {
#     return(as.character(data))
#   } else 
#   
#   if (is.factor(data)) {
#     if (length)
#   }
# }
# 

#' Transformation functions
#' 
#' Functions used by plotting funtions to transform data.
#' Calling the function with no parameters returns available function names.
#' 
#' @param data (\code{numeric}) Data to transform
#' @param func (\code{character}) Name of transformation to apply.
transform_data <- function(data = NULL, func = NULL) {
  funcs <- list("radius" = function(x) {x},
                "area" = function(x) {(x/pi)^(1/2)},
                "log10 radius" = function(x) {log(x, base = 10)},
                "log2 radius" = function(x) {log(x, base = 2)},
                "ln radius" = function(x) {log(x)},
                "log10 area" = function(x) {log((x/pi)^(1/2), base = 10)},
                "log2 area" = function(x) {log((x/pi)^(1/2), base = 2)},
                "ln area" =  function(x) {log((x/pi)^(1/2))})
  if (is.null(data) & is.null(func)) {
    return(names(funcs))
  } else {
    return(vapply(X = data, FUN = funcs[[func]], FUN.VALUE = numeric(1)))
  }
}



#' The defualt quantative color palette
#' 
#' The default color palette for quantative data
quantative_palette <- function() {
  return(c("grey", "#018571", "#80cdc1", "#dfc27d", "#a6611a"))
}


#' The defualt qualitative color palette
#' 
#' The default color palette for qualitative data
qualitative_palette <- function() {
  return(c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(9, "Pastel1")))
}