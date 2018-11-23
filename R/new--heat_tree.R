#' Plot a taxonomic tree
#'
#' Plots the distribution of values associated with a taxonomic classification/heirarchy. Taxonomic
#' classifications can have multiple roots, resulting in multiple trees on the same plot. The color
#' and size of tree elements, such as nodes (circles), edges (lines between circles), and lables can
#' used to display any data associated with taxa. These are similar to aesthetics in ggplot2. There
#' are many options, perhaps an overwheliming number, but typically only a few need to be used to
#' make an effective plot.
#'
#' @param obj An object of type \code{\link[taxa]{Taxmap}}. This typically contains all of the data
#'   that will be plotted. Any data that is to be plotted that is not contained in this object will
#'   need to named by the taxon ids used in this object or be the same length as the number of taxa
#'   in this object.
#'
#' @param node_label Labels associated with nodes. See the "labels" section below for more
#'   information. Default: no node labels.
#' @param edge_label Labels associated with edges. See the "labels" section below for more
#'   information. Default: no edge labels.
#' @param tree_label Labels associated with subtrees. When there are multiple roots to the taxonomy,
#'   multiple trees are printed, one for each root, and they are automatically arranged on the same
#'   graph using the same color/size scales. This option allows each subtree to be labeled. Values
#'   assocaited with the roots and given to this option will be used for the subtree labels. See the
#'   "labels" section below for more information. Default: no subtree labels.
#'
#' @param node_size A numeric value used to determine relative node size. See the "size" section
#'   below for more information. Default: constant size.
#' @param edge_size A numeric value used to determine relative edge size. Although edges connect two
#'   nodes, in this function, their relative width is associated with a single taxon, represented by
#'   the node the edge leads to. For this reason, root taxa do not have edges. Keep in mind that
#'   edges will look odd if they are bigger than their assocaited node. See the "size" section below
#'   for more information. Default: half of node size.
#' @param node_label_size A numeric value used to determine relative node label size. See the "size"
#'   section below for more information. Default: relative to node size.
#' @param edge_label_size A numeric value used to determine relative edge label size. See the "size"
#'   section below for more information. Default: relative to edge size.
#' @param tree_label_size A numeric value used to determine relative tree label size. See the "size"
#'   section below for more information. Default: relative to subtree size.
#'
#' @param node_color A numeric value used to determine relative node color or a character value
#'   representing a color, either hex code or color name. See the "color" section below for more
#'   information. Default: grey.
#' @param edge_color  A numeric value used to determine relative edge color or a character value
#'   representing a color, either hex code or color name. See the "color" section below for more
#'   information. Default: same as node color.
#' @param node_label_color A numeric value used to determine relative node label color or a
#'   character value representing a color, either hex code or color name. See the "color" section
#'   below for more information.. Default: black.
#' @param edge_label_color A numeric value used to determine relative edge label color or a
#'   character value representing a color, either hex code or color name. See the "color" section
#'   below for more information. Default: black.
#' @param tree_label_color A numeric value used to determine relative subtree label color or a
#'   character value representing a color, either hex code or color name. See the "color" section
#'   below for more information. Default: black.
#'
#' @param node_size_trans See details on transformations. Default: \code{"area"}.
#' @param edge_size_trans See details on transformations. Default: same as \code{node_size_trans}.
#' @param node_label_size_trans See details on transformations. Default: same as
#'   \code{node_size_trans}.
#' @param edge_label_size_trans See details on transformations. Default: same as
#'   \code{edge_size_trans}.
#' @param tree_label_size_trans See details on transformations. Default: \code{"area"}.
#' @param node_color_trans See details on transformations. Default: \code{"linear"}.
#' @param edge_color_trans See details on transformations. Default: same as node color
#'   transformation.
#' @param node_label_color_trans See details on transformations. Default: \code{"linear"}.
#' @param edge_label_color_trans See details on transformations. Default: \code{"linear"}.
#' @param tree_label_color_trans See details on transformations. Default: \code{"linear"}.
#'
#'
#' @section labels:
#'
#'   Inputs will be coerced into characters with \code{\link[base]{as.character}}. If a value for an
#'   individual taxon does not exist, or it is \code{NA}, or \code{""}, that label will not be
#'   printed. Labels can be multiple lines when \code{\n} is used to go to the next line. When
#'   labels are present, the color and size of labels can be used to display numeric information
#'   using \code{*_label_color} and \code{*_label_size} options. By default, labels are moved when
#'   needed to avoid overlap.
#'
#' @section sizes:
#'
#'   The size of nodes, edges, and labels can be mapped to numeric values. Non-numeric inputs will
#'   be coerced into characters with \code{\link[base]{as.numeric}}. If a value for an individual
#'   taxon does not exist, or it is \code{NA}, that element will not be printed. This is useful for
#'   displaying statistics for taxa, such as abundance. Only the relative size of the values given
#'   is used, not the values themselves. The *_size_trans (transformation) options can be used to
#'   change how values are transformed before being linearly mapped to color or size. The
#'   *_size_range options determine that absolute size displayed (in proportion of graph size). The
#'   *_size_interval options can be used to change what range of values will be graphically
#'   represented as the same size as the minimum/maximum *_size_range.
#'
#' @section colors:
#'
#'   The colors of nodes, edges, and labels can be mapped to \code{numeric} values. If the input is
#'   numeric, the relative size of values is mapped to a color ramp defined by the *_color_range
#'   options. If the input is \code{character}, then it must be either a value returned by
#'   \code{\link[grDevices]{colors}} or
#'   \href{https://en.wikipedia.org/wiki/Web_colors#Hex_triplet}{hexadecimal color codes}, or a
#'   mixture of the two. If a value for an individual taxon does not exist, or it is \code{NA}, that
#'   element will not be printed. The *_color_trans (transformation) options can be used to make the
#'   color mapping non-linear. The *_color_range options define the series of colors used to make
#'   the color ramp. The *_color_interval options can be used to change the range of input values
#'   that will be graphically represented as the same color as the minimum/maximum defined by the
#'   *_color_range options.
#'
#' @section transformations:
#'
#'   Before any conditions specified are mapped to an element property (color/size), they can be
#'   transformed to make the mapping non-linear. Any of the transformations listed below can be used
#'   by specifying their name. A customized function can also be supplied to do the transformation.
#'   
#'   \describe{
#'     \item{"linear"}{Proportional to radius/diameter of node}
#'     \item{"area"}{circular area; better perceptual accuracy than \code{"linear"} for sizes}
#'     \item{"log10"}{Log base 10 of radius}
#'     \item{"log2"}{Log base 2 of radius}
#'     \item{"ln"}{Log base e of radius}
#'     \item{"log10 area"}{Log base 10 of circular area}
#'     \item{"log2 area"}{Log base 2 of circular area}
#'     \item{"ln area"}{Log base e of circular area}
#'   }
#'   
