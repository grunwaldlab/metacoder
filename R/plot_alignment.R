#===================================================================================================
#' Display sequence alignment
#' 
#' Make a plot af a sequence alignment for an overview of alignment structure. 
#' 
#' @param alignment (\code{DNAbin}) A matrix representing a sequence alignment. 
#' 
#' @references ColorBrewer2 was used for the color palette
#' 
#' @return A \code{\link[ggplot2]{ggplot}} object
#' 
#' @examples
#' library(ape)
#' data(woodmouse)
#' plot_alignment(woodmouse)
#' 
#' @export
plot_alignment <- function(alignment) {
  color_key <- c("A" = "#a6cee3", "T" = "#1f78b4", "C" = "#b2df8a", "G" = "#33a02c", "-" = "#DDDDDD")
  alignment <- as.character(alignment)
  molten_alignment <- reshape::melt.matrix(alignment)
  names(molten_alignment) <- c("name", "position", "base")
  molten_alignment$base <-  toupper(molten_alignment$base)
  molten_alignment$color <- color_key[molten_alignment$base]
  molten_alignment$color[is.na(molten_alignment$color)] <- "#FFFFFF"
  ggplot2::ggplot(molten_alignment, ggplot2::aes_string(x = "position", y = "name")) +
    ggplot2::geom_tile(fill = molten_alignment$color) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), 
          panel.background = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          axis.text  = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(), 
          axis.line  = ggplot2::element_blank())
}