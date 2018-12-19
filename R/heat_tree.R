#' @rdname heat_tree
#' @export
heat_tree <- function(...) {
  UseMethod("heat_tree")
}


#' @param .input An object of type \code{\link{taxmap}}
#' 
#' @method heat_tree Taxmap
#' @export
#' @rdname heat_tree
heat_tree.Taxmap <- function(.input, ...) {
  # Non-standard argument evaluation
  data <- .input$data_used(...)
  data <- lapply(data,
                 function(x) { # orders everthing the same
                   if (is.null(names(x))) {
                     return(x)
                   } else {
                     return(x[.input$edge_list$to])
                   }
                 })
  arguments <- c(list(taxon_id = .input$edge_list$to, supertaxon_id = .input$edge_list$from),
                 lazyeval::lazy_eval(lazyeval::lazy_dots(...), data = data))
  
  # Check for common mistakes
  # func_not_vars <- c("taxon_name", "taxon_rank")
  invalid_input <- vapply(arguments, is.function, logical(1))
  if (any(invalid_input)) {
    stop(paste0("Function given to parameter '", names(invalid_input)[invalid_input][1], "'",
                " Did you use taxon_name/taxon_rank instead of taxon_names/taxon_ranks perhaps?"))
  }
  
  # Use variable name for scale axis labels
  if (! "node_color_axis_label" %in% names(arguments)) {
    arguments$node_color_axis_label <- deparse(as.list(match.call())$node_color)
  }
  if (! "node_size_axis_label" %in% names(arguments)) {
    arguments$node_size_axis_label <- deparse(as.list(match.call())$node_size)
  }
  if (! "edge_color_axis_label" %in% names(arguments)) {
    arguments$edge_color_axis_label <- deparse(as.list(match.call())$edge_color)
  }
  if (! "edge_size_axis_label" %in% names(arguments)) {
    arguments$edge_size_axis_label <- deparse(as.list(match.call())$edge_size)
  }
  
  # Call heat_tree
  do.call(heat_tree, arguments)
}




#' Plot a taxonomic tree
#' 
#' Plots the distribution of values associated with a taxonomic classification/heirarchy.
#' Taxonomic classifications can have multiple roots, resulting in multiple trees on the same plot.
#' A tree consists of elements, element properties, conditions, and mapping properties which are
#' represented as parameters in the heat_tree object.
#' The elements (e.g. nodes, edges, lables, and individual trees) are the infrastructure of the heat tree.
#' The element properties (e.g. size and color) are characteristics that are manipulated by various 
#' data conditions and mapping properties.  The element properties can be explicitly defined or automatically generated.
#' The conditions are data (e.g. taxon statistics, such as abundance) represented in the taxmap/metacoder object.
#' The mapping properties are parameters (e.g. transformations, range, interval, and layout) used to change the 
#' elements/element properties and how they are used to represent (or not represent) the various conditions.
#' 
#' @param taxon_id The unique ids of taxa.
#' @param supertaxon_id The unique id of supertaxon \code{taxon_id} is a part of.
#' 
#' @param node_label See details on labels.
#' Default: no labels.
#' @param edge_label See details on labels.
#' Default: no labels.
#' @param tree_label See details on labels.
#' The label to display above each graph.
#' The value of the root of each graph will be used.
#' Default: None.
#' 
#' @param node_size See details on size.
#' Default: constant size.
#' @param edge_size See details on size.
#' Default: relative to node size. 
# #' @param tree_size See details on size.
# #' The value of the root of each graph will be used.
# #' This scales the space used to display graphs, but does not effect node/edge size.
# #' Default: Not used. 
#' 
#' @param node_label_size See details on size.
#' Default: relative to vertex size.
#' @param edge_label_size See details on size.
#' Default: relative to edge size.
#' @param tree_label_size See details on size.
#' Default: relative to graph size.
#' 
#' @param node_color See details on colors.
#' Default: grey.
#' @param edge_color See details on colors.
#' Default: same as node color.
#' @param tree_color See details on colors.
#' The value of the root of each graph will be used.
#' Overwrites the node and edge color if specified.
#' Default: Not used.
#' 
#' @param node_label_color See details on colors.
#' Default: black.
#' @param edge_label_color See details on colors.
#' Default: black.
#' @param tree_label_color See details on colors.
#' Default: black.
#' 
#' @param node_size_trans See details on transformations.
#' Default: \code{"area"}.
#' @param edge_size_trans See details on transformations. 
#' Default: same as \code{node_size_trans}. 
# #' @param tree_size_trans See details on transformations.
# #' Default: \code{"area"}.
#' 
#' @param node_label_size_trans See details on transformations. 
#' Default: same as \code{node_size_trans}.
#' @param edge_label_size_trans See details on transformations. 
#' Default: same as \code{edge_size_trans}.
#' @param tree_label_size_trans See details on transformations.
#' Default: \code{"area"}.
#' 
#' @param node_color_trans See details on transformations. 
#' Default: \code{"area"}.
#' @param edge_color_trans See details on transformations.
#' Default: same as node color transformation.
#' @param tree_color_trans See details on transformations.
#' Default: \code{"area"}.
#' 
#' @param node_label_color_trans See details on transformations.
#' Default: \code{"area"}.
#' @param edge_label_color_trans See details on transformations.
#' Default: \code{"area"}.
#' @param tree_label_color_trans See details on transformations. 
#' Default: \code{"area"}.
#' 
#' @param node_size_range See details on ranges.
#' Default: Optimize to balance overlaps and range size.
#' @param edge_size_range See details on ranges.
#' Default: relative to node size range. 
# #' @param tree_size_range See details on ranges.
# #' Default: Not set.
#' 
#' @param node_label_size_range See details on ranges.
#' Default: relative to node size. 
#' @param edge_label_size_range See details on ranges.
#' Default: relative to edge size.
#' @param tree_label_size_range See details on ranges.
#' Default: relative to tree size.
#' 
#' @param node_color_range See details on ranges.
#' Default: Color-blind friendly palette. 
#' @param edge_color_range See details on ranges.
#' Default: same as node color.
#' @param tree_color_range See details on ranges.
#' Default: Color-blind friendly palette. 
#' 
#' @param node_label_color_range See details on ranges.
#' Default: Color-blind friendly palette. 
#' @param edge_label_color_range See details on ranges.
#' Default: Color-blind friendly palette.
#' @param  tree_label_color_range See details on ranges.
#' Default: Color-blind friendly palette. 
#' 
#' @param node_size_interval See details on intervals.
#' Default: The range of values in \code{node_size}. 
#' @param node_color_interval See details on intervals.
#' Default: The range of values in \code{node_color}. 
#' @param edge_size_interval See details on intervals.
#' Default: The range of values in \code{edge_size}. 
#' @param edge_color_interval See details on intervals.
#' Default: The range of values in \code{edge_color}. 
#' 
#' 
#' @param node_label_max The maximum number of node labels.
#' Default: 20.
#' @param edge_label_max The maximum number of edge labels.
#' Default: 20.
#' @param tree_label_max The maximum number of tree labels.
#' Default: 20.
#' 
#' @param overlap_avoidance (\code{numeric})
#' The relative importance of avoiding overlaps vs maximizing size range.
#' Higher numbers will cause node size optimization to avoid overlaps more.
#' Default: \code{1}.
#' 
#' @param margin_size (\code{numeric} of length 2)
#' The horizontal and vertical margins. c(left, right, bottom, top).
#' Default: \code{0, 0, 0, 0}.
#' 
#' @param layout The layout algorithm used to position nodes.
#' See details on layouts.
#' Default: \code{"reingold-tilford"}.
#' @param initial_layout he layout algorithm used to set the initial position
#' of nodes, passed as input to the \code{layout} algorithm.
#' See details on layouts.
#' Default: Not used.
#' @param make_node_legend if TRUE, make legend for node size/color mappings. 
#' @param make_edge_legend if TRUE, make legend for edge size/color mappings. 
#' @param title Name to print above the graph.
#' @param title_size The size of the title relative to the rest of the graph. 
#' 
#' @param node_color_axis_label The label on the scale axis corresponding to \code{node_color}.
#' Default: The expression given to \code{node_color}.
#' @param node_size_axis_label The label on the scale axis corresponding to \code{node_size}.
#' Default: The expression given to \code{node_size}.
#' @param edge_color_axis_label The label on the scale axis corresponding to \code{edge_color}.
#' Default: The expression given to \code{edge_color}.
#' @param edge_size_axis_label The label on the scale axis corresponding to \code{edge_size}.
#' Default: The expression given to \code{edge_size}.
#' 
#' @param background_color The background color of the plot.
#' Default: Transparent
#' @param output_file The path to one or more files to save the plot in using \code{\link[ggplot2]{ggsave}}. 
#' The type of the file will be determined by the extension given.
#' Default: Do not save plot.
#' 
#' @param aspect_ratio The aspect_ratio of the plot.
#' @param repel_labels If \code{TRUE} (Default), use the ggrepel package to spread out labels.
#' @param repel_force The force of which overlapping labels will be repelled from eachother. 
#' @param repel_iter The number of iterations used when repelling labels
#' @param verbose If \code{TRUE} print progress reports as the function runs.
#' 
#' @param ... (other named arguments)
#' Passed to the \code{\link{igraph}} layout function used.
#' 
#' 
#' @section labels:
#' 
#' The labels of nodes, edges, and trees can be added.
#' Node labels are centered over their node.
#' Edge labels are displayed over edges, in the same orientation.
#' Tree labels are displayed over their tree.
#' 
#' Accepts a vector, the same length \code{taxon_id} or a factor of its length.
#'
#' @section sizes:
#' 
#' The size of nodes, edges, labels, and trees can be mapped to various conditions.
#' This is useful for displaying statistics for taxa, such as abundance.
#' Only the relative size of the condition is used, not the values themselves.
#' The <element>_size_trans (transformation) parameter can be used to make the size mapping non-linear.
#' The <element>_size_range parameter can be used to proportionately change the size of an
#' element based on the condition mapped to that element.
#' The <element>_size_interval parameter can be used to change the limit at which a condition
#' will be graphically represented as the same size as the minimum/maximum <element>_size_range.
#' 
#' Accepts a \code{numeric} vector, the same length \code{taxon_id} or a
#' factor of its length.
#' 
#' @section colors:
#' 
#' The colors of nodes, edges, labels, and trees can be mapped to various conditions.
#' This is useful for visually highlighting/clustering groups of taxa.
#' Only the relative size of the condition is used, not the values themselves.
#' The <element>_color_trans (transformation) parameter can be used to make the color mapping non-linear.
#' The <element>_color_range parameter can be used to proportionately change the color of an
#' element based on the condition mapped to that element.
#' The <element>_color_interval parameter can be used to change the limit at which a condition
#' will be graphically represented as the same color as the minimum/maximum <element>_color_range.
#' 
#' Accepts a vector, the same length \code{taxon_id} or a factor of its length.
#' If a numeric vector is given, it is mapped to a color scale.
#' Hex values or color names can be used (e.g. \code{#000000} or \code{"black"}).
#' 
#' Mapping Properties
#'
#' @section transformations:
#' 
#' Before any conditions specified are mapped to an element property (color/size), they 
#' can be transformed to make the mapping non-linear. 
#' Any of the transformations listed below can be used by specifying their name.
#' A customized function can also be supplied to do the transformation.
#' 
#' \describe{
#'   \item{"linear"}{Proportional to radius/diameter of node}
#'   \item{"area"}{circular area; better perceptual accuracy than \code{"linear"}}
#'   \item{"log10"}{Log base 10 of radius}
#'   \item{"log2"}{Log base 2 of radius}
#'   \item{"ln"}{Log base e of radius}
#'   \item{"log10 area"}{Log base 10 of circular area}
#'   \item{"log2 area"}{Log base 2 of circular area}
#'   \item{"ln area"}{Log base e of circular area}
#' }
#' 
#' @section ranges:
#' 
#' The displayed range of colors and sizes can be explicitly defined or automatically generated.
#' When explicitly used, the size range will proportionately increase/decrease the size of a particular element.
#' Size ranges are specified by supplying a \code{numeric} vector with two values: the minimum and maximum.
#' The units used should be between 0 and 1, representing the proportion of a dimension of the graph.
#' Since the dimensions of the graph are determined by layout, and not always square, the value
#' that \code{1} corresponds to is the square root of the graph area (i.e. the side of a square with 
#' the same area as the plotted space).
#' Color ranges can be any number of color values as either HEX codes (e.g. \code{#000000}) or
#' color names (e.g. \code{"black"}).
#' 
#' @section layout:
#' 
#' Layouts determine the position of node elements on the graph.
#' They are implemented using the \code{\link{igraph}} package.
#' Any additional arguments passed to \code{heat_tree} are passed to the  \code{\link{igraph}}
#' function used.
#' The following \code{character} values are understood:
#' 
#' \describe{
#'   \item{"automatic"}{Use \code{\link[igraph]{nicely}}. Let \code{\link{igraph}} choose the layout.}
#'   \item{"reingold-tilford"}{Use \code{\link[igraph]{as_tree}}. A circular tree-like layout.}
#'   \item{"davidson-harel"}{Use \code{\link[igraph]{with_dh}}. A type of simulated annealing.}
#'   \item{"gem"}{Use \code{\link[igraph]{with_gem}}. A force-directed layout.}
#'   \item{"graphopt"}{Use \code{\link[igraph]{with_graphopt}}. A force-directed layout.}
#'   \item{"mds"}{Use \code{\link[igraph]{with_mds}}. Multidimensional scaling.}
#'   \item{"fruchterman-reingold"}{Use \code{\link[igraph]{with_fr}}. A force-directed layout.}
#'   \item{"kamada-kawai"}{Use \code{\link[igraph]{with_kk}}. A layout based on a physical model of springs.}
#'   \item{"large-graph"}{Use \code{\link[igraph]{with_lgl}}. Meant for larger graphs.}
#'   \item{"drl"}{Use \code{\link[igraph]{with_drl}}. A force-directed layout.}
#' }
#' 
#' 
#' @section intervals:
#' 
#' This is the minimum and maximum of values displayed on the legend scales.
#' Intervals are specified by supplying a \code{numeric} vector with two values: the minimum and maximum.
#' When explicitly used, the <element>_<property>_interval will redefine the way the actual conditional values are being represented
#' by setting a limit for the <element>_<property>.
#' Any condition below the minimum <element>_<property>_interval will be graphically represented the same as a condition AT the
#' minimum value in the full range of conditional values.  Any value above the maximum <element>_<property>_interval will be graphically 
#' represented the same as a value AT the maximum value in the full range of conditional values.
#' By default, the minimum and maximum equals the <element>_<property>_range used to infer the value of the <element>_<property>.
#' Setting a custom interval is useful for making <element>_<properties> in multiple graphs correspond to the same conditions,
#' or setting logical boundaries (such as \code{c(0,1)} for proportions.
#' Note that this is different from the <element>_<property>_range mapping property, which determines the size/color of graphed elements.
#' 
#' @section Acknowledgements:
#' 
#' This package includes code from the R package ggrepel to handle label overlap
#' avoidance with permission from the author of ggrepel Kamil Slowikowski. We
#' included the code instead of depending on ggrepel because we are using
#' internal functions to ggrepel that might change in the future. We thank Kamil
#' Slowikowski for letting us use his code and would like to acknowledge his
#' implementation of the label overlap avoidance used in metacoder.
#' 
#' @examples
#' \dontrun{
#' # Parse dataset for plotting
#' x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
#'                    class_key = c(tax_rank = "info", tax_name = "taxon_name"),
#'                    class_regex = "^(.+)__(.+)$")
#'                    
#' # Default appearance:
#' #  No parmeters are needed, but the default tree is not too useful
#' heat_tree(x)
#' 
#' # A good place to start:
#' #  There will always be "taxon_names" and "n_obs" variables, so this is a 
#' #  good place to start. This will shown the number of OTUs in this case. 
#' heat_tree(x, node_label = taxon_names, node_size = n_obs, node_color = n_obs)
#' 
#' # Plotting read depth:
#' #  To plot read depth, you first need to add up the number of reads per taxon.
#' #  The function `calc_taxon_abund` is good for this. 
#' x$data$taxon_counts <- calc_taxon_abund(x, data = "tax_data")
#' x$data$taxon_counts$total <- rowSums(x$data$taxon_counts[, -1]) # -1 = taxon_id column
#' heat_tree(x, node_label = taxon_names, node_size = total, node_color = total)
#' 
#' # Plotting multiple variables:
#' #  You can plot up to 4 quantative variables use node/edge size/color, but it
#' #  is usually best to use 2 or 3. The plot below uses node size for number of
#' #  OTUs and color for number of reads and edge size for number of samples
#' x$data$n_samples <- calc_n_samples(x, data = "taxon_counts")
#' heat_tree(x, node_label = taxon_names, node_size = n_obs, node_color = total,
#'           edge_color = n_samples)
#' 
#' # Different layouts:
#' #  You can use any layout implemented by igraph. You can also specify an
#' #  initial layout to seed the main layout with.
#' heat_tree(x, node_label = taxon_names, node_size = n_obs, node_color = n_obs,
#'           layout = "davidson-harel")
#' heat_tree(x, node_label = taxon_names, node_size = n_obs, node_color = n_obs,
#'           layout = "davidson-harel", initial_layout = "reingold-tilford")
#' 
#' # Axis labels:
#' #  You can add custom labeles to the legends
#' heat_tree(x, node_label = taxon_names, node_size = n_obs, node_color = total,
#'           edge_color = n_samples, node_size_axis_label = "Number of OTUs", 
#'           node_color_axis_label = "Number of reads",
#'           edge_color_axis_label = "Number of samples")
#'           
#' # Overlap avoidance:
#' #  You can change how much node overlap avoidance is used.
#' heat_tree(x, node_label = taxon_names, node_size = n_obs, node_color = n_obs,
#'           overlap_avoidance = .5)
#'           
#' # Label overlap avoidance
#' #  You can modfiy how label scattering is handled using the `replel_force` and
#' `repel_iter` options. You can turn off label scattering using the `repel_labels` option.
#' heat_tree(x, node_label = taxon_names, node_size = n_obs, node_color = n_obs,
#'           repel_force = 2, repel_iter = 20000)
#' heat_tree(x, node_label = taxon_names, node_size = n_obs, node_color = n_obs,
#'           repel_labels = FALSE)
#' 
#' # Setting the size of graph elements: 
#' #  You can force nodes, edges, and lables to be a specific size/color range instead
#' #  of letting the function optimize it. These options end in `_range`.
#' heat_tree(x, node_label = taxon_names, node_size = n_obs, node_color = n_obs,
#'           node_size_range = c(0.01, .1))
#' heat_tree(x, node_label = taxon_names, node_size = n_obs, node_color = n_obs,
#'           edge_color_range = c("black", "#FFFFFF"))
#' heat_tree(x, node_label = taxon_names, node_size = n_obs, node_color = n_obs,
#'           node_label_size_range = c(0.02, 0.02))
#' 
#' # Setting the transformation used:
#' #  You can change how raw statistics are converted to color/size using options
#' #  ending in _trans.
#' heat_tree(x, node_label = taxon_names, node_size = n_obs, node_color = n_obs,
#'           node_size_trans = "log10 area")
#' 
#' # Setting the interval displayed:
#' #  By default, the whole range of the statistic provided will be displayed.
#' #  You can set what range of values are displayed using options ending in `_interval`.
#' heat_tree(x, node_label = taxon_names, node_size = n_obs, node_color = n_obs,
#'           node_size_interval = c(10, 100))
#' 
#' }
#' @method heat_tree default
#' @rdname heat_tree
heat_tree.default <- function(taxon_id, supertaxon_id, 
                              node_label = NA,
                              edge_label = NA,
                              tree_label = NA,
                              
                              node_size = 1,
                              edge_size = node_size,
                              # tree_size = 1,
                              
                              node_label_size = node_size,
                              edge_label_size = edge_size,
                              tree_label_size = as.numeric(NA), 
                              
                              node_color = "#999999",
                              edge_color = node_color,
                              tree_color = NA,
                              
                              node_label_color = "#000000",
                              edge_label_color = "#000000",
                              tree_label_color = "#000000",
                              
                              node_size_trans = "area",
                              edge_size_trans = node_size_trans,
                              # tree_size_trans = "area",
                              
                              node_label_size_trans = node_size_trans,
                              edge_label_size_trans = edge_size_trans,
                              tree_label_size_trans = "area",
                              
                              node_color_trans = "area",
                              edge_color_trans = node_color_trans,
                              tree_color_trans = "area",
                              
                              node_label_color_trans = "area",
                              edge_label_color_trans = "area",
                              tree_label_color_trans = "area",
                              
                              node_size_range = c(NA, NA),
                              edge_size_range = c(NA, NA),
                              # tree_size_range = c(NA, NA),
                              
                              node_label_size_range = c(NA, NA),
                              edge_label_size_range = c(NA, NA),
                              tree_label_size_range = c(NA, NA),
                              
                              node_color_range = quantative_palette(),
                              edge_color_range = node_color_range,
                              tree_color_range = quantative_palette(),
                              
                              node_label_color_range = quantative_palette(),
                              edge_label_color_range = quantative_palette(),
                              tree_label_color_range = quantative_palette(),
                              
                              node_size_interval = range(node_size, na.rm = TRUE, finite = TRUE),
                              node_color_interval = NULL,
                              edge_size_interval = range(edge_size, na.rm = TRUE, finite = TRUE),
                              edge_color_interval = NULL,
                              
                              node_label_max = 500,
                              edge_label_max = 500,
                              tree_label_max = 500,
                              
                              overlap_avoidance = 1,
                              margin_size = c(0, 0, 0, 0),
                              layout = "reingold-tilford",
                              initial_layout = "fruchterman-reingold",
                              make_node_legend = TRUE,
                              make_edge_legend = TRUE,
                              title = NULL,
                              title_size = 0.08,
                              
                              node_color_axis_label = NULL, 
                              node_size_axis_label = NULL,
                              edge_color_axis_label = NULL, 
                              edge_size_axis_label = NULL,
                              
                              background_color = "#FFFFFF00",
                              output_file = NULL,
                              
                              aspect_ratio = 1,
                              repel_labels = TRUE,
                              repel_force = 1,
                              repel_iter = 1000,
                              
                              verbose = FALSE,
                              
                              ...) {
  #| ### Verify arguments =========================================================================
  if (length(taxon_id) != length(supertaxon_id)) {
    stop("'taxon_id' and 'supertaxon_id' must be of equal length.")
  }
  if (length(taxon_id) == 0) {
    warning("'taxon_id' and 'supertaxon_id' are empty. Returning NULL.")
    return(NULL)
  }
  if (length(unique(taxon_id)) != length(taxon_id)) {
    stop("All values of 'taxon_id' are not unique.")
  }
  check_element_length(c("node_size", "edge_size",# "tree_size",
                         "node_label_size", "edge_label_size",  "tree_label_size",
                         "node_color", "edge_color", "tree_color",
                         "node_label_color", "edge_label_color", "tree_label_color",
                         "node_label", "edge_label", "tree_label"))
  look_for_na(taxon_id, 
              c("node_size", "edge_size",
                "node_label_size", "edge_label_size",  "tree_label_size",
                "node_color", "edge_color", "tree_color",
                "node_label_color", "edge_label_color", "tree_label_color",
                "node_label", "edge_label", "tree_label"))
  verify_size(c("node_size", "edge_size", #"tree_size",
                "node_label_size", "edge_label_size", "tree_label_size"))
  verify_size_range(c("node_size_range",  "edge_size_range", # "tree_size_range",
                      "node_label_size_range", "edge_label_size_range", "tree_label_size_range",
                      "node_size_interval", "edge_size_interval"))
  verify_trans(c("node_size_trans", "edge_size_trans", #"tree_size_trans",
                 "node_color_trans", "edge_color_trans", "tree_color_trans",
                 "node_label_size_trans", "edge_label_size_trans", "tree_label_size_trans", 
                 "node_label_color_trans", "edge_label_color_trans", "tree_label_color_trans"))
  verify_color_range(c("node_color_range", "edge_color_range", "tree_color_range",
                       "node_label_color_range", "edge_label_color_range", "tree_label_color_range"))
  verify_label_count(c("node_label_max", "edge_label_max", "tree_label_max"))
  if (length(overlap_avoidance) == 0 || ! is.numeric(overlap_avoidance)) {
    stop("Argument 'overlap_avoidance' must be a numeric of length 1.")
  }
  if (length(margin_size) != 4 || ! is.numeric(margin_size)) {
    stop("Argument 'margin_size' must be a numeric of length 4: c(left, right, bottom, top)")
  }
  layout <- match.arg(layout, layout_functions())
  if (!is.null(initial_layout)) {
    initial_layout <- match.arg(initial_layout, layout_functions())
  }
  
  #| ### Parse arguments
  
  if (is.null(node_color_interval)) {
    if (! is.null(edge_color_interval) && all(node_color == edge_color)) {
      node_color_interval <- edge_color_interval
    } else if (length(get_numerics(node_color)) > 0) {
      node_color_interval <- range(get_numerics(node_color),
                                   na.rm = TRUE, finite = TRUE)
    }
  }
  
  if (is.null(edge_color_interval)) {
    if (! is.null(node_color_interval) && all(node_color == edge_color)) {
      edge_color_interval <- node_color_interval
    } else if (length(get_numerics(edge_color)) > 0) {
      edge_color_interval <- range(get_numerics(edge_color),
                                   na.rm = TRUE, finite = TRUE)
    }
  }
  
  #| ### Standardize source data ==================================================================
  data <- data.frame(stringsAsFactors = FALSE,
                     tid_user = as.character(taxon_id),
                     pid_user = as.character(supertaxon_id),
                     
                     vl_user = as.character(node_label),
                     el_user = as.character(edge_label),
                     tl_user = as.character(tree_label),
                     
                     vs_user = as.numeric(node_size),
                     es_user = as.numeric(edge_size),
                     # ts_user = as.numeric(tree_size),
                     
                     vls_user = as.numeric(node_label_size),
                     els_user = as.numeric(edge_label_size),
                     tls_user = as.numeric(tree_label_size),
                     
                     vc_user = node_color,
                     ec_user = edge_color,
                     tc_user = tree_color,
                     
                     vlc_user = node_label_color,
                     elc_user = edge_label_color,
                     tlc_user = tree_label_color)
  row.names(data) <- data$tid_user
  
  #| #### Apply statistic transformations =========================================================
  trans_key <- c(vs_user = node_size_trans, es_user = edge_size_trans, #ts_user = tree_size_trans,
                 vls_user = node_label_size_trans, els_user = edge_label_size_trans,  tls_user = tree_label_size_trans,
                 vc_user = node_color_trans, ec_user = edge_color_trans, tc_user = edge_color_trans,
                 vlc_user = node_label_color_trans, elc_user = edge_label_color_trans, tlc_user = tree_label_color_trans)
  transformed_names <- gsub(pattern = "_user$", x = names(trans_key), replacement = "_trans")
  apply_trans <- function(col_name) {
    if (is.numeric(data[ , col_name])) { 
      transform_data(trans_key[col_name], data[ , col_name]) # if numbers are supplied
    } else {
      data[ , col_name] # if colors are defined explicitly, then no transformation is done
    }
  }
  data[, transformed_names] <- lapply(names(trans_key), apply_trans)
  # transform intervals
  node_size_interval_trans <- transform_data(node_size_trans, node_size_interval)
  edge_size_interval_trans <- transform_data(edge_size_trans, edge_size_interval)
  node_color_interval_trans <- transform_data(node_color_trans, node_color_interval)
  edge_color_interval_trans <- transform_data(edge_color_trans, edge_color_interval)
  
  
  #| ### Make layout ==============================================================================
  #| The layout is used to generate a list of coordinates to places graph verticies
  #| First the edge list consituted by the `taxon_id` and `supertaxon_id` columns is used to construct 
  #| an `igraph` graph object and then the layout is generated for that object. 
  #|
  #| #### Make a graph for each root in the graph -------------------------------------------------
  my_print("Calculating layout for ", nrow(data), " taxa...", verbose = verbose)
  get_sub_graphs <- function(taxa) {
    if (length(taxa) == 1) {
      # Make a graph with only a single node
      adj_matrix <- matrix(c(0), ncol = 1, dimnames =  list(taxa, taxa))
      sub_graph <- igraph::graph.adjacency(adj_matrix)
    } else {
      # Make edge list from taxon_id and supertaxon_id
      edgelist <- as.matrix(data[taxa, c("pid_user", "tid_user")])
      # Remove edges to taxa that dont exist in this subset of the dataset
      edgelist <- edgelist[! is.na(edgelist[, "pid_user"]), , drop = FALSE]
      #       # Randomly resort if layout is "reingold-tilford". NOTE: This is kinda hackish and should be replaced
      #       if (layout == "reingold-tilford") { 
      #         grouped_index <- split(rownames(edgelist), f = edgelist[, "pid_user"])
      #         grouped_index <- unlist(grouped_index[sample(seq_along(grouped_index))])
      #         edgelist <- edgelist[grouped_index, , drop = FALSE] 
      #       }
      sub_graph <- igraph::graph_from_edgelist(edgelist)
    }
    igraph::V(sub_graph)$weight_factor <- data[taxa, c("vs_trans")]
    edge_end_node <- gsub("^[0-9]+\\|", "", attr(igraph::E(sub_graph), "vnames"))
    igraph::E(sub_graph)$weight_factor <- data[edge_end_node, c("vs_trans")]
    return(sub_graph)
  }
  data$is_root <- !(data$pid_user %in% data$tid_user)
  data[data$is_root, "pid_user"] <- NA # Needed by split_by_level
  sub_graph_taxa <- split_by_level(data$tid_user, data$pid_user, level =  1)
  sub_graphs <- lapply(sub_graph_taxa, get_sub_graphs)
  #|
  #| #### Generate a layout for each graph --------------------------------------------------------
  #|
  get_sub_layouts <- function(graph, backup_layout = 'fruchterman-reingold') {
    # Calculate an initial layout if specified
    if (! is.null(initial_layout) && layout != initial_layout) {
      intitial_coords <- layout_functions(initial_layout, graph)
      if (! any(is.na(intitial_coords) | is.nan(unlist(intitial_coords)))) {
        intitial_coords <- rescale(intitial_coords, to = ((nrow(intitial_coords) ^ 0.65) + 5) * c(-1, 1))
        # intitial_coords <- rescale(intitial_coords, to = c(-100, 100))
      }
    } else {
      intitial_coords <- NULL
    }
    # Calculate the primary layout 
    coords <- layout_functions(layout, graph, intitial_coords = intitial_coords, ...)
    # Calculate backup layout if primary one does not work
    if (any(is.na(coords) | is.nan(unlist(coords)))) {
      coords <- layout_functions(backup_layout, graph)
      warning(paste0("Could not apply layout '", layout,
                     "' to subgraph. Using 'fruchterman-reingold' instead."))
    }
    return(coords)
  }
  
  sub_coords <- lapply(sub_graphs, get_sub_layouts)
  subgraph_key <- stats::setNames(rep(names(sub_graph_taxa), vapply(sub_graph_taxa, length, numeric(1))),
                                  unlist(sub_graph_taxa))
  data$subgraph_root <- subgraph_key[data$tid_user]
  #|
  #| #### Merge layout coordinates into an overall graph ------------------------------------------
  #|
  coords <- igraph::merge_coords(sub_graphs, sub_coords) # merge node coordinates for each tree
  graph <- igraph::disjoint_union(sub_graphs) # merge graphs of each tree
  row.names(coords) <- names(igraph::V(graph))
  data$vx_plot <- coords[data$tid_user, 1]
  data$vy_plot <- coords[data$tid_user, 2]
  
  # Rescale to constant size
  my_range <- range(c(data$vx_plot, data$vy_plot))
  to_scale <- c(0, 1)
  data$vx_plot <- rescale(data$vx_plot, to = to_scale, from = my_range)
  data$vy_plot <- rescale(data$vy_plot, to = to_scale, from = my_range)
  
  #| ### Set aspect ration
  data$vx_plot <- data$vx_plot * aspect_ratio
  
  #| ### Core plot data ===========================================================================
  #|
  #| #### Optimize node size range --------------------------------------------------------------
  #|
  if (any(is.na(node_size_range))) {
    my_print("Optmizing node size range...", verbose = verbose)
  }
  # Get range of potential node size ranges - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (nrow(data) > 1) {
    all_pairwise <- molten_dist(x = data$vx_plot, y = data$vy_plot) # get distance between all nodes
    x_diff <- max(data$vx_plot) - min(data$vx_plot)
    y_diff <- max(data$vy_plot) - min(data$vy_plot)
    square_side_length <- sqrt(x_diff * y_diff)
    if (is.na(node_size_range[1])) { # if minimum node size not set
      min_range <- c(0, min(all_pairwise$distance))
    } else {
      min_range <- rep(node_size_range[1], 2) * square_side_length
    }
    if (is.na(node_size_range[2])) { # if maximum node size not set
      max_range <- c(min_range[1], square_side_length / 5)
    } else {
      max_range <- c(node_size_range[2], 2) * square_side_length
    }
    if (! is.na(node_size_range[1]) && ! is.na(node_size_range[2])) {
      vsr_plot <- node_size_range * square_side_length
    } else {
      # Subset pairwise pairs to increase speed - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      max_important_pairs <- 1000 # Takes into account both size and distance
      max_biggest_pairs <- 1000
      max_closest_pairs <- 1000
      if (nrow(all_pairwise) > sum(max_important_pairs, max_biggest_pairs, max_closest_pairs)) {
        all_pairwise$size_sum <- data$vs_trans[all_pairwise$index_1] + data$vs_trans[all_pairwise$index_2]
        all_pairwise$importance <- all_pairwise$size_sum / all_pairwise$distance
        # all_pairwise <- all_pairwise[order(all_pairwise$importance, decreasing = TRUE), ]
        pair_subset <- c(order(all_pairwise$importance, decreasing = TRUE)[1:max_important_pairs],
                         order(all_pairwise$size_sum, decreasing = TRUE)[1:max_biggest_pairs],
                         order(all_pairwise$distance)[1:max_closest_pairs])
        all_pairwise <- all_pairwise[pair_subset, ]
      }
      find_overlap <- function(a_min, a_max, distance) {
        scaled_vs <- rescale(data$vs_t, to = c(a_min, a_max), from = node_size_interval_trans)
        names(scaled_vs) <- data$tid_user
        gap <- distance$distance - scaled_vs[distance$index_1] - scaled_vs[distance$index_2]
        overlap <- ifelse(gap < 0, abs(gap), 0)
        overlap <- (overlap ^ 2) / (scaled_vs[distance$index_1] ^ 2 + scaled_vs[distance$index_2] ^ 2)
        mean(overlap)
      }
      
      # Choose base range based on optimality criteria  - - - - - - - - - - - - - - - - - - - - - - - -
      optimality_stat <- function(minimum, maximum) {
        if (minimum == 0) {
          overlap <- 0
        } else {
          overlap <- find_overlap(minimum, maximum, all_pairwise)
        }
        ideal_min <- 0.02
        ideal_max <- 0.3
        ideal_range <- .1
        minimum <- minimum / square_side_length
        maximum <- maximum / square_side_length
        min_size_score <- min(c(1, 1 - (ideal_min - minimum) / ideal_min))
        max_size_score <- min(c(1, 1 - (maximum - ideal_max) / maximum))
        range_prop <- minimum / maximum 
        range_size_score <- min(c(1, 1 - abs(range_prop - ideal_range)))
        overlap_score <- min(c(1, 1 - overlap ^ (0.08 / overlap_avoidance) )) # Totally observation based; might need to be rethought
        result <- prod(c(min_size_score, max_size_score, range_size_score, overlap_score))
        # print(c(min_size_score, max_size_score, range_size_score, overlap_score, result))
        return(result)
      }
      
      # Use genetic algorithm to pick range
      ga_result <- GA::ga(type = "real-valued", 
                          fitness =  function(x) optimality_stat(x[1], x[2]),
                          lower = c(min_range[1], max_range[1]), upper = c(min_range[2], max_range[2]),
                          maxiter = 40, run = 30, popSize = 70, monitor = FALSE, parallel = FALSE)
      vsr_plot <- as.vector(ga_result@solution[1, ])
    }
  } else {
    square_side_length = 1
    vsr_plot <- rep(square_side_length / 4, 2)
  }
  data$vs_plot <- rescale(data$vs_t, to = vsr_plot, from = node_size_interval_trans)
  #|
  #| #### Infer edge size range -------------------------------------------------------------------
  #|
  infer_size_range <- function(specified_range, reference_range, default_scale) {
    result <- specified_range * square_side_length
    if (is.na(result[1]) && is.na(result[2])) { # If the user has not set range
      result <- reference_range * default_scale
    } else if (is.na(result[1])) { # If the user has set a maximum but not a minimum
      result[1] <- result[2] * reference_range[1] / reference_range[2]
    } else if (is.na(result[2])) { # If the user has set a minimum but not a maximum
      result[2] <- result[1] * reference_range[2] / reference_range[1]
    }
    return(result)
  }
  
  esr_plot <- infer_size_range(edge_size_range, vsr_plot, default_scale = 0.5)
  data$es_plot <- rescale(data$es_t, to = esr_plot, from = edge_size_interval_trans)
  #|
  #| #### Infer tree size range -------------------------------------------------------------------
  #|
  get_tree_area <- function(a_root) {
    size <- data[data$subgraph_root == a_root, "vs_plot"]
    x <- data[data$subgraph_root == a_root, "vx_plot"]
    x <- c(x + size, x - size)
    y <- data[data$subgraph_root == a_root, "vy_plot"]
    y <- c(y + size, y - size)
    (max(x) - min(x)) * (max(y) - min(y)) 
  }
  tree_area <- vapply(unique(data$subgraph_root), get_tree_area, FUN.VALUE = numeric(1))
  data$tree_area <- tree_area[data$subgraph_root]
  tsr_plot <- range(sqrt(tree_area))
  #|
  #| #### Infer label size ranges -----------------------------------------------------------------
  #|
  if (all(is.na(data$tls_user))) {
    data$tls_user <- sqrt(data$tree_area)
    data$tls_trans <- apply_trans("tls_user") 
  }
  vlsr_plot <- infer_size_range(node_label_size_range, vsr_plot, default_scale = 0.8)
  elsr_plot <- infer_size_range(edge_label_size_range, esr_plot, default_scale = 0.8)
  tlsr_plot <- infer_size_range(tree_label_size_range, tsr_plot, default_scale = 0.1)
  data$vls_plot <- rescale(data$vls_trans, to = vlsr_plot)
  data$els_plot <- rescale(data$els_trans, to = elsr_plot)
  data$tls_plot <- rescale(data$tls_trans, to = tlsr_plot)
  #|
  #| #### Assign color scales ---------------------------------------------------------------------
  #|
  
  color_colume_key <- list("ec_trans" = edge_color_range, "vc_trans" = node_color_range, 
                           "tc_trans" = tree_color_range, "vlc_trans" = node_label_color_range,
                           "elc_trans" = edge_label_color_range, "tlc_trans" = tree_label_color_range)
  color_interval_key <- list("ec_trans" = edge_color_interval_trans, "vc_trans" = node_color_interval_trans)
  plot_value_names <- gsub(pattern = "_trans$", x = names(color_colume_key), replacement = "_plot")
  data[, plot_value_names] <- lapply(names(color_colume_key),
                                     function(x) apply_color_scale(data[ , x],
                                                                   color_colume_key[[x]],
                                                                   interval = color_interval_key[[x]]))
  # If tree_color is used, overwrite other colors - - - - - - - - - - - - - - - - - - - - - - - - -
  data$tc_plot <- data[data$subgraph_root, "tc_plot"]
  to_replace <- ! is.na(data$tc_plot)
  data[to_replace, "vc_plot"] <- data[to_replace, "tc_plot"]
  data[to_replace, "ec_plot"] <- data[to_replace, "tc_plot"]
  
  #| ### Secondary plot data ======================================================================
  #|
  #| #### Calculate coordinants of graph elements -------------------------------------------------
  #| The nodes and edges must be specified by a dataframe of coordinates, with a colume 
  #| grouping the coordinates of each shape.
  #| These shapes must be added to the graph in a specific order.
  #| A list of nodes is sorted by first node depth in the heirarchy and then by node size.
  taxon_elements <- function(tid) {
    circle_resolution <- 35
    edge_data <- line_coords(x1 = data[tid, 'vx_plot'],
                             y1 = data[tid, 'vy_plot'],
                             x2 = data[data[tid, 'pid_user'], "vx_plot"],
                             y2 = data[data[tid, 'pid_user'], "vy_plot"],
                             width = data[tid, 'es_plot'] * 2)
    edge_data$group <- paste0(tid, "_edge")
    edge_data$color <- rep(data[tid, 'ec_plot'], each = 4)
    node_data <- polygon_coords(n = circle_resolution,
                                x = data[tid, 'vx_plot'],
                                y = data[tid, 'vy_plot'],
                                radius = data[tid, 'vs_plot'])
    node_data$group <- paste0(tid, "_node")
    node_data$color <- rep(data[tid, 'vc_plot'], each = circle_resolution + 1)
    output <- rbind(edge_data, node_data)
    # output$tid_user <- tid
    return(output[stats::complete.cases(output),])
  }
  data$level = edge_list_depth(data$tid_user, data$pid_user)
  element_order <- data$tid_user[order(data$level, 1 / data$vs_plot, decreasing = TRUE)]
  element_data <- do.call(rbind, lapply(element_order, taxon_elements))
  element_data$group <- factor(element_data$group, levels = unique(element_data$group))
  #|
  #| #### Make text data ------------------------------------------------------------------
  #|
  # Get node label data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data$vl_is_shown <- select_labels(data, node_label_max,
                                    sort_by_column = c("vls_plot", "vs_plot"),
                                    label_column = "vl_user")
  if (any(data$vl_is_shown)) {
    vl_data <- data[data$vl_is_shown, , drop = FALSE]
    text_data <- data.frame(stringsAsFactors = FALSE,
                            label = vl_data$vl_user,
                            x = vl_data$vx_plot,
                            y = vl_data$vy_plot,
                            size = vl_data$vls_plot,
                            color = vl_data$vlc_plot,
                            rotation = 0,
                            justification = "center",
                            group = "nodes")
  } else {
    text_data <- NULL
  }
  # Get edge label data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data$el_is_shown <- select_labels(data, edge_label_max,
                                    sort_by_column = c("els_plot", "es_plot"),
                                    label_column = "el_user")
  data[is.na(data$pid_user), "el_is_shown"] <- FALSE # taxa with no parents get no line label
  if (any(data$el_is_shown)) {
    el_data <- data[data$el_is_shown, ]
    # edge label rotation 
    el_data$el_slope <- (el_data$vy_plot - data[el_data$pid_user, "vy_plot"]) / (el_data$vx_plot - data[el_data$pid_user, "vx_plot"])
    el_data$el_slope[is.na(el_data$el_slope)] <- 0
    el_data$el_rotation <- atan(el_data$el_slope)
    # edge label coordinate 
    line_label_offset = 1
    justify <- data[el_data$pid_user, "vx_plot"] > el_data$vx_plot
    justify[is.na(justify)] <- TRUE
    justification <- ifelse(justify, "left-center", "right-center")
    line_label_x_offset <- line_label_offset * el_data$vs_plot * cos(el_data$el_rotation)
    line_label_y_offset <- line_label_offset * el_data$vs_plot * sin(el_data$el_rotation)
    el_data$elx_plot <- el_data$vx_plot + ifelse(justify, 1, -1) * line_label_x_offset
    el_data$ely_plot <- el_data$vy_plot + ifelse(justify, 1, -1) * line_label_y_offset
    # create text data   
    text_data <- rbind(text_data,
                       data.frame(stringsAsFactors = FALSE, 
                                  label = el_data$el_user,
                                  x = el_data$elx_plot,
                                  y = el_data$ely_plot,
                                  size = el_data$els_plot,
                                  color = el_data$elc_plot,
                                  rotation = el_data$el_rotation,
                                  justification = justification,
                                  group = "edges"))
  }
  # Get tree label data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data$tl_is_shown <- FALSE
  data[data$is_root, "tl_is_shown"] <- select_labels(data[data$is_root, ], tree_label_max,
                                                     sort_by_column = c("tls_plot", "vs_plot"), label_column = "tl_user")
  if (any(data$tl_is_shown)) {
    title_data <- data[data$tl_is_shown, , drop = FALSE]
    tx_plot <- vapply(split(data$vx_plot, data$subgraph_root), FUN.VALUE = numeric(1),
                      function(x) mean(range(x)))
    title_data$tx_plot <- tx_plot[title_data$subgraph_root]
    ty_plot <- vapply(split(data$vy_plot, data$subgraph_root), FUN.VALUE = numeric(1),
                      function(y) mean(range(y)))
    title_data$ty_plot <- ty_plot[title_data$subgraph_root]
    title_data$tlx_plot <- title_data$tx_plot 
    tly_plot <- mapply(function(y, size) max(y + size),
                       y = split(data$vy_plot, data$subgraph_root),
                       size = split(data$vs_plot, data$subgraph_root))
    title_data$tly_plot <- tly_plot[title_data$subgraph_root] + title_data$tls_plot * 1.1
    text_data <- rbind(text_data,
                       data.frame(stringsAsFactors = FALSE, 
                                  label = title_data$tl_user,
                                  x = title_data$tlx_plot,
                                  y = title_data$tly_plot,
                                  size = title_data$tls_plot,
                                  color = title_data$tlc_plot,
                                  rotation = 0,
                                  justification = "center",
                                  group = "trees"))
  }
  
  
  
  # Get range data ---------------------------------------------------------------------------------
  get_limits <- function() {
    label_corners <- label_bounds(label = text_data$label, x = text_data$x, y = text_data$y,
                                  height = text_data$size, rotation = text_data$rotation,
                                  just = text_data$justification)
    x_points <- c(element_data$x, label_corners$x)
    y_points <- c(element_data$y, label_corners$y)
    margin_size_plot <- margin_size * square_side_length
    x_range <- c(min(x_points) - margin_size_plot[1], max(x_points) + margin_size_plot[2]) 
    y_range <- c(min(y_points) - margin_size_plot[3], max(y_points) + margin_size_plot[4]) 
    return(list(x = x_range, y = y_range))
  }
  
  # Add tree title data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ranges <- get_limits() # UGLY HACK! FIX!
  if (! is.null(title)) {
    title_size <- diff(ranges$x) * title_size
    text_data <- rbind(text_data,
                       data.frame(stringsAsFactors = FALSE, 
                                  label = title,
                                  x = mean(ranges$x),
                                  y = max(ranges$y) + title_size * 0.5,
                                  size = title_size,
                                  color = "#000000",
                                  rotation = 0,
                                  justification = "center-bottom",
                                  group = "title"))
  }
  
  
  # Repel labels
  ranges <- get_limits() # UGLY HACK! FIX!
  reformat_bounds <- function(bounds) {
    if (is.null(bounds)) {
      return(NULL)
    }
    bounds$label <- factor(bounds$label, levels=unique(bounds$label)) # keep order when split
    x_coords <- split(bounds$x, rep(seq_len(nrow(text_data)), each = 4))
    y_coords <- split(bounds$y, rep(seq_len(nrow(text_data)), each = 4))
    data.frame(label = text_data$label,
               color = text_data$color,
               rotation = rad_to_deg(text_data$rotation),
               group = text_data$group,
               xmin = vapply(x_coords, min, numeric(1)),
               xmax = vapply(x_coords, max, numeric(1)),
               ymin = vapply(y_coords, min, numeric(1)),
               ymax = vapply(y_coords, max, numeric(1)),
               stringsAsFactors = FALSE)
  }
  
  if (!is.null(text_data)) {
    bounds <- label_bounds(label = text_data$label, x = text_data$x, y = text_data$y,
                           height = text_data$size, rotation = text_data$rotation,
                           just = text_data$justification)
    bounds <- reformat_bounds(bounds)
    
    if (repel_labels) {
      movable <- text_data$group != "legend"
      text_data[movable, c("x", "y")] <- repel_boxes(data_points = as.matrix(text_data[movable, c("x", "y")]),
                                                     boxes = as.matrix(bounds[movable, c("xmin", "ymin", "xmax", "ymax")]),
                                                     point_padding_x = 0, point_padding_y = 0,
                                                     xlim = ranges$x,
                                                     ylim = ranges$y,
                                                     # hjust = 0.5,
                                                     # vjust = 0.5,
                                                     force = 1e-06 * repel_force,
                                                     maxiter = repel_iter,
                                                     direction = "both")
      bounds <- label_bounds(label = text_data$label, x = text_data$x, y = text_data$y,
                             height = text_data$size, rotation = text_data$rotation,
                             just = text_data$justification)
      bounds <- reformat_bounds(bounds)
    }
  } else {
    bounds <- NULL
  }
  
  
  
  #|
  #| #### Make node legend -----------------------------------------------------------------------
  #|
  my_print("Making legends...", verbose = verbose)
  if (make_node_legend | make_edge_legend) {
    legend_length <- square_side_length * 0.3 
    
    # right_plot_boundry <- max(c(element_data[element_data$y <= legend_length + min(element_data$y), "x"],
    #       bounds[bounds$ymin <= legend_length +  min(element_data$y), "xmax"]))
    
    if (is.null(bounds)) {
      right_plot_boundry <- max(element_data$x)
    } else {
      right_plot_boundry <- max(c(element_data$x, bounds$xmax))
    }
    
    if (make_node_legend) {
      node_legend <- make_plot_legend(x = right_plot_boundry,
                                      y = min(element_data$y) * 0.9, 
                                      length = legend_length, 
                                      width_range = vsr_plot * 2, 
                                      width_trans_range = range(data$vs_trans) * 2,
                                      width_stat_range =  node_size_interval,
                                      group_prefix = "node_legend",
                                      width_stat_trans = transform_data(func = node_size_trans, inverse = TRUE),
                                      color_range = node_color_range,
                                      color_trans_range = node_color_interval_trans,
                                      color_stat_range = node_color_interval, 
                                      color_stat_trans =  transform_data(func = node_color_trans, inverse = TRUE),
                                      title = "Nodes",
                                      color_axis_label = node_color_axis_label,
                                      size_axis_label = node_size_axis_label,
                                      hide_size = missing(node_size),
                                      hide_color = missing(node_color))
      element_data <- rbind(element_data, node_legend$shapes)
      text_data <- rbind(text_data, node_legend$labels)
    }
    #|
    #| #### Make edge legend -----------------------------------------------------------------------
    #|
    
    # right_plot_boundry <- max(c(element_data[element_data$y >= max(element_data$y) - legend_length, "x"],
    #                             bounds[bounds$ymax >= max(element_data$y) - legend_length, "xmin"]))
    
    if (make_edge_legend) {
      edge_legend <- make_plot_legend(x = right_plot_boundry,
                                      y = max(element_data$y) - legend_length * 1.3, 
                                      length = legend_length, 
                                      width_range = esr_plot * 2, 
                                      width_trans_range = range(data$vs_trans) * 2,
                                      width_stat_range =  edge_size_interval,
                                      group_prefix = "edge_legend",
                                      width_stat_trans = transform_data(func = edge_size_trans, inverse = TRUE),
                                      color_range = edge_color_range,
                                      color_trans_range = edge_color_interval_trans,
                                      color_stat_range = edge_color_interval, 
                                      color_stat_trans =  transform_data(func = edge_color_trans, inverse = TRUE),
                                      title = "Edges",
                                      color_axis_label = edge_color_axis_label,
                                      size_axis_label = edge_size_axis_label,
                                      hide_size = missing(edge_size),
                                      hide_color = missing(edge_color))
      element_data <- rbind(element_data, edge_legend$shapes)
      text_data <- rbind(text_data, edge_legend$labels)
    }
    bounds <- label_bounds(label = text_data$label, x = text_data$x, y = text_data$y,
                           height = text_data$size, rotation = text_data$rotation,
                           just = text_data$justification)
    bounds <- reformat_bounds(bounds)
  } else {
    legend_data <- NULL
  }
  
  
  
  
  #| ### Draw plot ================================================================================
  # text_boxes <-  label_bounds(label = text_data$label, x = text_data$x, y = text_data$y,  # debug
  #                                       height = text_data$size, rotation = text_data$rotation,
  #                                       just = text_data$justification)
  # text_boxes$group <- rep(seq_along(text_data$label), each = 4)  # debug
  my_print("Plotting graph...", verbose = verbose)
  result = tryCatch({
    
    ranges <- get_limits()
    the_plot <- ggplot2::ggplot(data = data) +
      ggplot2::geom_polygon(data = element_data, ggplot2::aes_string(x = "x", y = "y", group = "group"),
                            fill = element_data$color) +
      ggplot2::guides(fill = "none") +
      ggplot2::coord_fixed(xlim = ranges$x, ylim = ranges$y) +
      ggplot2::scale_y_continuous(expand = c(0,0), limits = ranges$y) + 
      ggplot2::scale_x_continuous(expand = c(0,0), limits = ranges$x) +
      # ggplot2::geom_polygon(data = text_boxes, mapping = ggplot2::aes(x = x, y = y, group = group), color = "black", fill = NA) + # debug
      ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                     panel.background = ggplot2::element_rect(fill = background_color, colour = background_color),
                     plot.background = ggplot2::element_rect(fill = background_color, colour = background_color),
                     axis.title = ggplot2::element_blank(),
                     axis.text  = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(), 
                     axis.line  = ggplot2::element_blank(),
                     plot.margin = grid::unit(c(0,0,0,0) , "in"))
    
    # Plot text..
    
    if (! is.null(bounds)) {
      bounds <- bounds[bounds$label != "" & ! is.na(bounds$label), ]
      if (nrow(bounds) > 0) {
        the_plot <- the_plot + ggfittext::geom_fit_text(data = bounds, 
                                                        grow = TRUE,
                                                        min.size = 0,
                                                        # reflow = TRUE,
                                                        color = bounds$color,
                                                        padding.x = grid::unit(0, "mm"),
                                                        padding.y = grid::unit(0, "mm"),
                                                        ggplot2::aes_string(label = "label",
                                                                            xmin = "xmin",
                                                                            xmax = "xmax",
                                                                            ymin = "ymin",
                                                                            ymax = "ymax",
                                                                            angle = "rotation")) 
      }
    }
    
    #| ### Save output file
    if (!is.null(output_file)) {
      img_width <- diff(ranges$x)
      img_height <- diff(ranges$y)
      for (path in output_file) {
        ggplot2::ggsave(path, the_plot, bg = "transparent", width = 10, height = 10 * (img_height / img_width))
      }
    }
    
    
  }, error = function(msg) {
    if (grepl(msg$message, "Error: evaluation nested too deeply: infinite recursion / options(expressions=)?", fixed = TRUE)) {
      stop(paste(msg, sep = "\n", 
                 "NOTE: This error typically occurs because of too many text labels being printed.", 
                 "You can avoid it by increasing the value of `expressions` in the global options:",
                 "    * How to see the current value: options('expressions')",
                 "    * How to increase the value:    options(expressions = 100000)"))
    } else {
      stop(msg)
    }
  })
  
  
  
  return(the_plot)
}

