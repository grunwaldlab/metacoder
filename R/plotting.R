#===================================================================================================
#' Display sequence alignment
#' 
#' Make a plot pf a sequence alignment for an overview of alignment structure. 
#' NOT FINISHED.
#' 
#' @param alignment (\code{DNAbin}) A matrix representing a sequence alignment. 
#' 
#' @references ColorBrewer2 was used for the color palette
plot_alignment <- function(alignment) {
  color_key <- c("A" = "#a6cee3", "T" = "#1f78b4", "C" = "#b2df8a", "G" = "#33a02c", "-" = "#DDDDDD")
  alignment <- as.character(alignment)
  molten_alignment <- reshape::melt.matrix(alignment)
  names(molten_alignment) <- c("name", "position", "base")
  molten_alignment$base <-  toupper(molten_alignment$base)
  molten_alignment$color <- color_key[molten_alignment$base]
  molten_alignment$color[is.na(molten_alignment$color)] <- "#FFFFFF"
  ggplot(molten_alignment, aes(x = position, y = name)) +
    geom_tile(fill = molten_alignment$color) +
    theme(panel.grid = element_blank(), 
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text  =  element_blank(),
          axis.ticks = element_blank(), 
          axis.line  = element_blank())
}