#' Estimate text grob length
#' 
#' Estimate the printed length of `resizingTextGrob` text
#' 
#' @param text \code{character} The text to be printed
#' @param rot The rotation in radians 
#' 
#' @return The estimated length of the printed text as a multiple of its text size (height)
#' 
#' @keywords internal 
text_grob_length <- function(text, rot = 0) {
  do_one <- function(text) {
    as.numeric(grid::widthDetails(grid::textGrob(text, rot = rot * 180 / pi))) / as.numeric(grid::heightDetails(grid::textGrob(text))) * .8
  }
  vapply(text, do_one, numeric(1))
 }

