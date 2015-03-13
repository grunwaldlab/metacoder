#===================================================================================================
#' Calculate kimura-2 from unaligned sequences
#' 
#' The two input sequences are aligned using Biostrings::pairwiseAlignment and the kimura
#' 2-parameter distance metric is calculated.
#' 
#' This is not set up to deal with IUPAC ambiguity codes. It will treat ambiguous SNPs as 
#' transversions. Use with caution; this is still in the development stage. 
#' 
#' @param x A character vector of DNA sequence. 
#' @param y A character vector of DNA sequence. 
#' @param ... Additional arguments are passed to Biostrings::pairwiseAlignment
#' @return A Kimura2 distance metric as a numeric vector of length 1.
model_kimura2 <- function(x, y, ...) {
  is_transition <- function(x, y) {
    x <- tolower(x)
    y <- tolower(y)
    setequal(c(x,y), c("a", "g")) | setequal(c(x,y), c("c", "t"))
  }
  # Standardize input ------------------------------------------------------------------------------
  x <- unlist(lapply(x, paste, collapse=""))
  y <- unlist(lapply(y, paste, collapse=""))
  # Align sequences --------------------------------------------------------------------------------
  alignment <- Biostrings::pairwiseAlignment(x, y, ...)
  # Calculate distance -----------------------------------------------------------------------------
  x <- seqinr::s2c(tolower(as.character(alignment@pattern)))
  y <- seqinr::s2c(tolower(as.character(alignment@subject)))
  transitions <- sum(unlist(Map(is_transition, x, y))) / length(x)
  equal <- sum(x == y) / length(x)
  tranversions <- 1 - equal - transitions
  -.5 * log(1 - 2 * transitions - tranversions) - .25 * log(1 - 2 * tranversions)
}


#===================================================================================================
#' Calculates a pairwise distance matrix from unaligned sequences
#' 
#' This is still in the developmental stage and currently only uses the kimura 2-parameter metric.
#' @param sequences A list or vector of sequence data
#' @return A distance matrix of class \code{dist}
#' @importFrom proxy dist
pairwise_distance <- function(sequences, seq_name = NULL) {
  sequences <- as.list(sequences)
  if (!is.null(seq_name)) names(sequences) <- seq_name
  dist(sequences, method = model_kimura2)
}