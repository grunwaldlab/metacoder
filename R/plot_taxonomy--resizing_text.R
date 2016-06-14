#' Estimate text grob length
#' 
#' Estimate the printed length of `resizingTextGrob` text
#' 
#' @param text \code{character} The text to be printed
#' 
#' @return The estimated length of the printed text as a multiple of its height
#' 
#' @keywords internal 
text_grob_length <- function(text) {
  do_one <- function(text) {
    text_grob <- grid::textGrob(text)
    as.numeric(grid::widthDetails(text_grob)) / as.numeric(grid::heightDetails(text_grob))
  }
  vapply(text, do_one, numeric(1))
 }




#' Create a list of text grobs
#' 
#' Creates a list of resizing text grobs
#' 
#' @param label Text to display
#' @param x x
#' @param y y
#' @param size The height of text
#' @param color The color of text
#' @param rotation The rotation of text. In radians.
#' @param justification A \code{character} indicating the justification of text. 
#' Is split by \code{"-"}
#' @param x_range The minimum and maximum x value in the plot 
#' @param y_range The minimum and maximum y value in the plot
#' 
#' @return \code{list} of \code{\link{resizingTextGrob}}
#' 
#' @keywords internal
make_text_grobs <- function(label, x, y, size, color, rotation, justification, x_range, y_range) {
  make_text_grob <- function(label, x, y, size, color, rotation, justification) {
    resizingTextGrob(label = label, x = x, y = y,
                     gp = grid::gpar(text_prop = size,
                                     col = color),
                     rot = rotation * 180 / pi,
                     just =  strsplit(justification, split = "-")[[1]])
  }
  square_side_length <- sqrt((x_range[2] - x_range[1]) * (y_range[2] - y_range[1]))
  x <- scales::rescale(x, to = c(0, 1), from = x_range)
  y <- scales::rescale(y, to = c(0, 1), from = y_range)
  size <- scales::rescale(size, to = c(0, 1), from = c(0, square_side_length))
  mapply(make_text_grob, label, x, y, size, color, rotation, justification, SIMPLIFY = FALSE)
}





#===================================================================================================
#' Adds text grob that scales with viewport size
#' 
#' Creates a new grob class called resizingTextGrob that is supposed to resize automatically
#' It is a thin layer over \code{grid::textGrob}.
#' 
#' Taken from: http://ryouready.wordpress.com/2012/08/01/creating-a-text-grob-that-automatically-adjusts-to-viewport-size/
#' 
#' @export
#' @keywords internal
resizingTextGrob <- function(...)
{
  grid::grob(tg=grid::textGrob(...), cl='resizingTextGrob')
}

#===================================================================================================
#' Draws a resizingTextGrob
#' 
#' Called automatically when drawing a grob using grob.draw.
#' 
#' Taken from: http://ryouready.wordpress.com/2012/08/01/creating-a-text-grob-that-automatically-adjusts-to-viewport-size/
#' 
#' @export
#' @keywords internal
drawDetails.resizingTextGrob <- function(x, recording=TRUE)
{
  grid::grid.draw(x$tg)
}

#===================================================================================================
#' Adjusts text size to viewport
#' 
#' Automatically called before any drawing occures.
#' 
#' Taken from: http://ryouready.wordpress.com/2012/08/01/creating-a-text-grob-that-automatically-adjusts-to-viewport-size/
#' 
#' @export
#' @keywords internal
preDrawDetails.resizingTextGrob <- function(x)
{
  cex <- x$tg$gp$text_prop * 4
  h <- grid::convertHeight(grid::unit(1, 'snpc'), 'mm', valueOnly=TRUE)
  grid::pushViewport(grid::viewport(gp = grid::gpar(fontsize = h * cex)))
}


#===================================================================================================
#' Clean up after the drawing.
#' 
#' Clean up after the drawing.
#' 
#' Taken from: http://ryouready.wordpress.com/2012/08/01/creating-a-text-grob-that-automatically-adjusts-to-viewport-size/
#' 
#' @export
#' @keywords internal
postDrawDetails.resizingTextGrob <- function(x) {
  grid::popViewport()
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