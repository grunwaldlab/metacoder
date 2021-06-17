#' Print something
#'
#' The standard print function for this package. This is a wrapper to make
#' package-wide changes easier.
#'
#' @param ... Something to print
#' @param verbose If \code{FALSE}, do not print anything.
#'
#' @keywords internal
my_print <- function(..., verbose = TRUE) {
  if (verbose) {
    text <- paste0(as.character(list(...)), collapse = "")
    message(text)
  }
}


#' get_edge_parents
#' 
#' @keywords internal
get_edge_parents <- function(graph) {
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
  
  truncate <- function(x, max_chars = 30, postfix = "[truncated]") {
    if (nchar(x) > max_chars) {
      x <- paste0(substr(x, 0, max_chars - nchar(postfix)), postfix)
    }
    return(x)
  }
  
  q = "'"
  interleaved <- interleave(chars[1:(length(chars) / 2)],
                            rev(chars[(length(chars) / 2 + 1):length(chars)]))
  is_greater_than_max <- cumsum(nchar(interleaved) + 2) + 10 > max_chars
  if (all(! is_greater_than_max)) {
    max_printed <- length(chars)
  } else {
    max_printed <- which.max(is_greater_than_max) - 1
  }
  if (max_printed < length(chars)) {
    if (max_printed < 2) {
      first_part <- truncate(chars[1])
      second_part <- truncate(chars[length(chars)])
    } else {
      first_part <-  chars[1:ceiling(max_printed / 2)]
      second_part <- chars[(length(chars) - floor(max_printed / 2) + 1):length(chars)]
    }
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


#' @keywords internal
rad_to_deg <- function(rad) {
  (rad * 180) / (pi)
}


#' @keywords internal
deg2rad <- function(deg) {
  (deg * pi) / (180)
}


#' Converts decimal numbers to other bases
#'
#' Converts from base 10 to other bases represented by a given set of symbols.
#'
#' @param numbers One or more numbers to convert.
#' @param symbols The set of symbols to use for the new base.
#' @param base The base to convert to.
#' @param min_length The minimum number of symbols in each result.
#'
#' @return character vector
#'
#' @keywords internal
convert_base <- function(numbers, symbols = letters, base = length(symbols),
                         min_length = 0) {
  
  # A modification of the `dec2base` function in the `oro.dicom` package
  #    Copyright (c) 2015, Brandon Whitcher
  convert_one <- function (n)  {
    if (is.na(n)) {
      return(NA_character_)
    }
    max_length <- max(trunc(log(max(n, 1))/log(base)) + 1, min_length)
    power <- rep(1, length(n)) * base^((max_length - 1):0)
    n <- n * rep(1, max_length)
    digits <- floor((n%%(base * power))/power)
    paste(symbols[digits + 1], collapse = "")
  }
  
  vapply(as.integer(numbers), convert_one, character(1))
  
}

#' Return github url
#' 
#' Return github url
#' @keywords internal
repo_url <- function() {
  "https://github.com/grunwaldlab/metacoder"  
}


#' Make a temporary file U's replaced with T
#' 
#' Make a temporary fasta file U's replaced with T without reading in whole file.
#' 
#' @param file_path
#' 
#' @return A path to a temporary file.
#' 
#' @keywords internal
make_fasta_with_u_replaced <- function(file_path) {
  
  # Create temporary file path
  output_path <- tempfile()
  output_con <- file(output_path, open = "w")
  
  # Open file for reading
  input_con <- file(file_path, open = "r")
  
  # Replace U's with T's one line at a time
  while(length(line <- readLines(input_con, n = 1)) > 0) {
    if (! startsWith(line, ">")) {
      line <- gsub(line, pattern = "U", replacement = "T", fixed = TRUE)
      line <- gsub(line, pattern = "u", replacement = "t", fixed = TRUE)
    }
    writeLines(line, output_con)
  }
  
  # Close file connections
  close(output_con)
  close(input_con)
  
  # Return path to temporary file
  return(output_path)
  
}
