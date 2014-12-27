#===================================================================================================
#' Adds text grob that scales with viewport size
#' 
#' Creates a new grob class called resizingTextGrob that is supposed to resize automatically
#' It is a thin layer over \code{grid::textGrob}.
#' 
#' Taken from: http://ryouready.wordpress.com/2012/08/01/creating-a-text-grob-that-automatically-adjusts-to-viewport-size/
#' 
#' @import grid
resizingTextGrob <- function(...)
{
  grob(tg=textGrob(...), cl='resizingTextGrob')
}

#===================================================================================================
#' Draws a resizingTextGrob
#' 
#' Called automatically when drawing a grob using grob.draw.
#' 
#' Taken from: http://ryouready.wordpress.com/2012/08/01/creating-a-text-grob-that-automatically-adjusts-to-viewport-size/
#' @import grid
drawDetails.resizingTextGrob <- function(x, recording=TRUE)
{
  grid.draw(x$tg)
}

#===================================================================================================
#' Adjusts text size to viewport
#' 
#' Automatically called before any drawing occures.
#' 
#' Taken from: http://ryouready.wordpress.com/2012/08/01/creating-a-text-grob-that-automatically-adjusts-to-viewport-size/
#' @import grid
#' @import scales
preDrawDetails.resizingTextGrob <- function(x)
{
  cex <- x$tg$gp$size * 4
  h <- convertHeight(unit(1, 'snpc'), 'mm', valueOnly=TRUE)
  pushViewport(viewport(gp = gpar(fontsize = h * cex)))
}


#===================================================================================================
#' Clean up after the drawing.
#' 
#' Clean up after the drawing.
#' 
#' Taken from: http://ryouready.wordpress.com/2012/08/01/creating-a-text-grob-that-automatically-adjusts-to-viewport-size/
#' 
#' @import grid
postDrawDetails.resizingTextGrob <- function(x) {
  popViewport()
}




# #===================================================================================================
# #' Resizing text annotations.
# #'
# #' @section Aesthetics:
# #' \Sexpr[results=rd,stage=build]{ggplot2:::rd_aesthetics("geom", "text")}
# #'
# #' @inheritParams geom_point
# #' @param parse If TRUE, the labels will be parsed into expressions and
# #' displayed as described in ?plotmath
# #' @export
# #' @import ggplot2
# #' @import grid
# geom_resizing_text <- function (mapping = NULL, data = NULL, stat = "identity", position = "identity",
#                        parse = FALSE, ...) {
#   GeomResizingText$new(mapping = mapping, data = data, stat = stat, position = position,
#                parse = parse, ...)
# }
#
# #===================================================================================================
# #' ggplot2 extension 
# #' 
# #' Resizing text annotations that scale with viewport size
# #' 
# #' @exportClass
# #' @import proto
# #' @import ggplot2
# #' @importClassesFrom ggplot2
# GeomResizingText <- proto::proto(Geom, {
#   objname <- "text"
#   draw_groups <- function(., ...) .$draw(...)
#   draw <- function(., data, scales, coordinates, ..., parse = FALSE, na.rm = FALSE) {
#     data <- remove_missing(data, na.rm,
#                            c("x", "y", "label"), name = "geom_text")
#     lab <- data$label
#     if (parse) {
#       lab <- parse(text = lab)
#     }
#     with(coord_transform(coordinates, data, scales),
#          resizingTextGrob(lab, x, y, default.units="native",
#                   hjust=hjust, vjust=vjust, rot=angle,
#                   gp = gpar(col = alpha(colour, alpha), size = size,
#                             fontfamily = family, fontface = fontface, lineheight = lineheight)))
#   }
#   draw_legend <- function(., data, ...) {
#     data <- aesdefaults(data, .$default_aes(), list(...))
#     with(data,
#          textGrob("a", 0.5, 0.5, rot = angle,
#                   gp=gpar(col=alpha(colour, alpha), fontsize = size * .pt))
#     )
#   }
#   default_stat <- function(.) StatIdentity
#   required_aes <- c("x", "y", "label")
#   default_aes <- function(.) aes(colour="black", size=5 , angle=0, hjust=0.5,
#                                  vjust=0.5, alpha = NA, family="", fontface=1, lineheight=1.2)
#   guide_geom <- function(x) "text"
# })