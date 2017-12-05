#' Print something
#'
#' The standard print function for this package. This is a wrapper to make
#' package-wide changes easier.
#'
#' @param ... Something to print
#'
#' @keywords internal
my_print <- function(...) {
  text <- paste0(as.character(list(...)), collapse = "")
  message(text)
}


#' get_edge_parents
#' 
#' @keywords internal
get_edge_parents <-function(graph) {
  igraph::get.edges(graph, 1:igraph::ecount(graph))[,1]
}


#' get_edge_children
#' 
#' @keywords internal
get_edge_children <- function(graph) {
  igraph::get.edges(graph, 1:igraph::ecount(graph))[,2]
}


#' get_node_children
#' 
#' @keywords internal
get_node_children <- function(graph, node) {
  which(igraph::shortest.paths(graph, igraph::V(graph)[node], mode="out") != Inf)
}


#' delete_vetices_and_children
#' 
#' @keywords internal
delete_vetices_and_children <- function(graph, nodes) {
  nodes <- unlist(sapply(nodes, function(x) get_node_children(graph, x)))
  graph <- igraph::delete.vertices(graph, nodes)
  return(graph)
}


#' add_alpha
#' 
#' @keywords internal
add_alpha <- function(col, alpha = 1){
  apply(sapply(col, grDevices::col2rgb) / 255, 2,
        function(x) grDevices::rgb(x[1], x[2], x[3], alpha = alpha))
}


#' get indexes of a unique set of the input
#' 
#' @keywords internal
unique_mapping <- function(input) {
  unique_input <- unique(input)
  vapply(input, function(x) {if (is.na(x)) which(is.na(unique_input)) else which(x == unique_input)}, numeric(1))
}


#' Run a function on unique values of a iterable 
#' 
#' @param input What to pass to \code{func}
#' @param func (\code{function})
#' @param ... passend to \code{func}
#' 
#' @keywords internal
map_unique <- function(input, func, ...) {
  input_class <- class(input)
  unique_input <- unique(input)
  class(unique_input) <- input_class
  func(unique_input, ...)[unique_mapping(input)]
}


#' Converts DNAbin to a named character vector
#' 
#' Converts an object of class DNAbin (as produced by ape) to a named character vector.
#' 
#' @param dna_bin (\code{DNAbin} of length 1) the input.
#' 
#' @keywords internal
DNAbin_to_char <- function(dna_bin) {
  vapply(as.character(dna_bin), paste, character(1), collapse="")
}

#' Test if characters can be converted to numbers
#' 
#' Makes TRUE/FALSE vector
#' 
#' @param input A character vector
#' 
#' @keywords internal
can_be_num <- function(input) {
  suppressWarnings(!is.na(as.numeric(input)))
}

#' Return numeric values in a character
#' 
#' Returns just valid numeric values and ignores others.
#' 
#' @param input
#' 
#' @keywords internal
get_numerics <- function(input) {
  as.numeric(input[can_be_num(input)])
}

#' check for packages
#'
#' check for packages, and stop if not installed.
#' This function was written by Scott Chamerlain, from whom I shamelessly stole
#' it.
#'
#' @param package The name of the package
#'
#' @return `TRUE` if package is present
#'
#' @keywords internal
check_for_pkg <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Please install ", package, call. = FALSE)
  } else {
    invisible(TRUE)
  }
}


#' Print a subset of a character vector
#'
#' Prints the start and end values for a character vector. The number of values
#' printed depend on the width of the screen by default.
#'
#' @param chars (`character`) What to print.
#' @param prefix (`character` of length 1) What to print before
#'   `chars`, on the same line.
#' @param max_chars (`numeric` of length 1) The maximum number of
#'   characters to print.
#' @param type (`"error"`, `"warning"`, `"message"`, `"cat"`, `"print"`, `"silent"``)
#'
#' @return `NULL`
#'
#' @examples
#' taxa:::limited_print(1:100)
#' taxa:::limited_print(1:10000)
#' taxa:::limited_print(1:10000, prefix = "stuff:")
#'
#' @keywords internal
limited_print <- function(chars, prefix = "",
                          max_chars = getOption("width") - nchar(prefix) - 5,
                          type = "message") {
  
  if (length(chars) == 0) {
    cat(prefix)
    return(invisible(NULL))
  }
  
  
  # https://stat.ethz.ch/pipermail/r-help/2006-March/101023.html
  interleave <- function(v1,v2) {
    ord1 <- 2*(1:length(v1))-1
    ord2 <- 2*(1:length(v2))
    c(v1,v2)[order(c(ord1,ord2))]
  }
  
  q = "'"
  interleaved <- interleave(chars[1:(length(chars) / 2)],
                            rev(chars[(length(chars) / 2 + 1):length(chars)]))
  is_greater_than_max <- cumsum(nchar(interleaved) + 2) + 10 > max_chars
  if (all(! is_greater_than_max)) {
    max_printed <- length(chars)
  } else {
    max_printed <- which.max(is_greater_than_max)
  }
  if (max_printed < length(chars)) {
    first_part <-  chars[1:as.integer(max_printed / 2 - 0.5)]
    second_part <-
      chars[as.integer(length(chars) - (max_printed / 2) + 1.5):length(chars)]
    output <- paste0(paste0(collapse = ", ", first_part),
                     " ... ",
                     paste0(collapse = ", ", second_part),
                     "\n")
  } else {
    output <- paste0(paste0(collapse = ", ", chars), "\n")
  }
  output <- paste(prefix, output, collapse = "")
  
  if (type == "error") {
    stop(output)
  } else if (type == "warning") {
    warning(output)
  } else if (type == "message") {
    message(output)
  } else if (type == "cat") {
    cat(output)
  } else if (type == "print") {
    print(output)
  } else if (type != "silent") {
    stop("invalid type option")
  }
  return(invisible(output))
}


#' Get a table from a taxmap object
#' 
#' Get a table from a taxmap object and complain if it does not exist.
#' 
#' @param obj A taxmap object
#' @param dataset The name of the table
#' @param expected_cols The names/indexes of columns expected to exist. If a column is not found, issue a warning.
#' 
#' @return A data.frame
#' 
#' @keywords internal
get_taxmap_table <- function(obj, dataset, expected_cols = NULL) {
  # Check that dataset exists and is a table
  if (! dataset %in% names(obj$data)) {
    stop(paste0('The dataset "', dataset,
                '" is not in the object supplied. Datasets found include:\n  ',
                limited_print(names(obj$data), type = "silent")), call. = FALSE)
  }
  if (! is.data.frame(obj$data[[dataset]])) {
    stop(paste0('The dataset "', dataset,  '" is not a table.'), call. = FALSE)
  }
  
  # Get table
  table <- obj$data[[dataset]]
  
  # Check that all columns exist
  if (! is.null(expected_cols)) {
    invalid_cols <- get_invalid_cols(table, expected_cols)
    if (length(invalid_cols) > 0) {
      warning(paste0('The following ', length(invalid_cols),
                     ' column(s) were not found in dataset "', dataset, '":\n',
                     limited_print(prefix = "  ", invalid_cols, type = "silent")),
              call. = FALSE)
    }
  }

  # Return without printing
  return(invisible(table))
}


#' Return invalid column names/indexes
#' 
#' Return invalid column names/indexes
#' 
#' @param table Table to check
#' @param cols he names/indexes of columns
#' 
#' @keywords internal
get_invalid_cols <- function(table, cols)  {
  if (is.numeric(cols)) { # if column indexes
    invalid_cols <- cols[! cols %in% seq_len(ncol(table))]
  } else { # if column names
    invalid_cols <- cols[! cols %in% colnames(table)]
  }
  return(invalid_cols)
}


#' @keywords internal
rad_to_deg <- function(rad) {
  (rad * 180) / (pi)
}


#' @keywords internal
deg2rad <- function(deg) {
  (deg * pi) / (180)
  }