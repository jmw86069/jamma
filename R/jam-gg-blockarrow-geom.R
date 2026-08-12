#' Block arrow or rectangle in the plot margins
#'
#' Draw rectangular blocks or block arrows in the plot margins,
#' suitable for annotating volcano plots with hit counts.
#' This is modeled after `ggh4x::geom_rectmargin()` but customized
#' for block arrow visualization.
#'
#' @inheritParams ggplot2::geom_rect
#' @param sides `character` string to control which sides of the
#'    plot the blocks appear. A string containing any letters in
#'    `"trbl"` will set it to top, right, bottom and left respectively.
#' @param direction `character` vector with direction the arrow
#'    should point. It recognizes letters in `'rlud'` as
#'    right, left, up, and down, respectively.
#' @param length `grid::unit()` object that sets the width and
#'    height of the rectangles in the x- and y-directions respectively.
#'    Default is 5% of the plot size, using `'snpc'` units so it
#'    maintains consistent width for vertical and horizontal
#'    block arrows. Change to `'npc'` to use fixed percentage
#'    of the plot size.
#' @param arrow_type `character` indicating the type of drawing:
#'    * `"blockarrow"` for block arrow shapes, default.
#'    * `"rect"` is implemented, a simple rectangle.
#' @param label `character` text label to display on the rectangle.
#'    Default `NULL`, or `NA` displays no label.
#' @param label_size `numeric` size for label text, default 18.
#' @param label_color `character` color for label text. Default NULL
#'    uses a contrasting color with `jamba::setTextContrastColor()`.
#' @param font `integer` to define the font face, default 2 uses
#'    bold font face, use 1 for normal.
#' @param na.rm `logical` not implemented.
#' @param show.legend `logical` or `NA`, not implemented.
#' @param inherit.aes `logical` not implemented.
#'
#' @return A *Layer* ggproto object.
#'
#' @examples
#' #\dontrun{
#' # Simple example with rectangles in margins
#' df <- data.frame(
#'   x = 1:10,
#'   y = rnorm(10),
#'   label = head(LETTERS, 10)
#' )
#' df
#' p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
#'   ggplot2::geom_point() +
#'   geom_blockmargin(inherit.aes=FALSE,
#'     # ggplot2::aes(xmin = 2, xmax = 5, ymin = 0.5, ymax = 1.5),
#'     xmin = 2, xmax = 5, ymin = 0.5, ymax = 1.5,
#'     label = c("Label A", "Label B", rep("", 8)),
#'     fill = c("royalblue", "red", rep(NA, 8)),
#'     alpha = 0.3, sides = c("b", "l", rep("b", 8))
#'   ) +
#'   ggplot2::coord_cartesian(clip = "off")
#'
#' print(p)
#' #}
#'
#' @export
geom_blockmargin <- function(
   mapping = NULL,
   data = NULL,
   stat = "identity",
   position = "identity",
   ...,
   sides = "bl",
   direction = "u",
   length = grid::unit(0.05, "snpc"),
   arrow_type = "blockarrow",
   label = NULL,
   label_size = 18,
   label_color = NULL,
   font = 2,
   na.rm = FALSE,
   show.legend = NA,
   inherit.aes = TRUE
) {
   ggplot2::layer(
      data        = data,
      mapping     = mapping,
      stat        = stat,
      geom        = GeomBlockMargin,
      position    = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params      = list(
         sides        = sides,
         direction    = direction,
         length       = length,
         arrow_type   = arrow_type,
         label        = label,
         label_size   = label_size,
         label_color  = label_color,
         font         = font,
         na.rm        = na.rm,
         ...
      )
   )
}

# ggproto object for GeomBlockMargin
#' @usage NULL
#' @format NULL
#' @keywords internal
#' @noRd
GeomBlockMargin <- ggplot2::ggproto(
   "GeomBlockMargin",
   ggplot2::GeomRect,
   draw_panel = function(
      self, data, panel_params, coord,
      sides = "bl",
      direction = "u",
      length = grid::unit(0.03, "snpc"),
      arrow_type = "rect",
      arrow_width = 0.9,
      arrow_head = 1,
      label = NULL,
      label_size = 10,
      label_color = NULL,
      font = 2
) {

      if (!inherits(length, "unit")) {
         stop("The `length` argument must be a grid unit object.")
      }

      rugs <- list()
      coords <- coord$transform(data, panel_params)
      sides <- coords$sides;
      direction <- coords$direction;

      if (inherits(coord, "CoordFlip")) {
         sides <- chartr("tblr", "rlbt", sides)
         direction <- chartr("udlr", "durl", direction)
      }

      # Rectangles are drawn inside the plot margins (between plot area and axis)
      rug_length <- list(
         min = length,
         max = grid::unit(1, "npc") - length
      )

      gp <- grid::gpar(
         col = coords$colour,
         fill = ggplot2::alpha(coords$fill, coords$alpha),
         lty = coords$linetype,
         lwd = coords$linewidth * ggplot2::.pt
      )
      # jamba::printDebug("coords:");print(coords);# debug

      # Draw rectangles on x-direction (top/bottom)
      # jamba::printDebug("sides:", sides);# debug
      if (!is.null(coords$xmin)) {
         if (!"rect" %in% arrow_type) {
            #
            if (any(grepl("t", sides))) {
               use_xmin <- grid::unit(coords$xmin, "native")
               use_xmax <- grid::unit(coords$xmax, "native")
               # inside margins
               use_ymax <- grid::unit(1, "npc")
               use_ymin <- grid::unit(1, "npc") - rug_length$min
               # # outside margins
               # use_ymax <- grid::unit(1, "npc") + rug_length$min
               # use_ymin <- grid::unit(1, "npc")
            } else if (any(grepl("b", sides))) {
               use_xmin <- grid::unit(coords$xmin, "native")
               use_xmax <- grid::unit(coords$xmax, "native")
               # inside margins
               use_ymax <- rug_length$min
               use_ymin <- grid::unit(0, "npc")
               # # outside margins
               # use_ymax <- grid::unit(0, "npc")
               # use_ymin <- -rug_length$min
            } else if (any(grepl("l", sides))) {
               use_ymin <- grid::unit(coords$ymin, "native")
               use_ymax <- grid::unit(coords$ymax, "native")
               # inside margins
               use_xmin <- grid::unit(0, "npc")
               use_xmax <- rug_length$min
               # # outside margins
               # use_xmin <- -rug_length$min
               # use_xmax <- grid::unit(0, "npc")
            } else if (any(grepl("r", sides))) {
               use_ymin <- grid::unit(coords$ymin, "native")
               use_ymax <- grid::unit(coords$ymax, "native")
               # inside margins
               use_xmin <- grid::unit(1, "npc") - rug_length$min
               use_xmax <- grid::unit(1, "npc")
               # # outside margins
               # use_xmin <- grid::unit(1, "npc")
               # use_xmax <- grid::unit(1, "npc") + rug_length$min
            }
            # print(list(use_xmin=use_xmin, use_xmax=use_xmax, use_ymin=use_ymin, use_ymax=use_ymax));# debug
            xy <- make_blockarrow_coords(
               xmin = use_xmin,
               xmax = use_xmax,
               ymin = use_ymin,
               ymax = use_ymax,
               direction = coords$direction
            )
            rugs$x_b <- grid::pathGrob(
               x = xy$x,
               y = xy$y,
               pathId = xy$id,
               gp = gp
            )
         } else if (grepl("b", sides)) {
            rugs$x_b <- grid::rectGrob(
               x = grid::unit(coords$xmin, "native"),
               y = grid::unit(0, "npc"),
               width = grid::unit(coords$xmax - coords$xmin, "native"),
               height = rug_length$min,
               just = c("left", "bottom"),
               gp = gp
            )
         } else if (grepl("t", sides)) {
            rugs$x_t <- grid::rectGrob(
               x = grid::unit(coords$xmin, "native"),
               y = grid::unit(1, "npc"),
               width = grid::unit(coords$xmax - coords$xmin, "native"),
               height = rug_length$min,
               just = c("left", "top"),
               gp = gp
            )
         }
      } else if (!is.null(coords$ymin)) {
         # Draw rectangles on y-direction (left/right)
         if (grepl("l", sides)) {
            rugs$y_l <- grid::rectGrob(
               x = grid::unit(0, "npc"),
               y = grid::unit(coords$ymax, "native"),
               width = rug_length$min,
               height = grid::unit(coords$ymax - coords$ymin, "native"),
               just = c("left", "top"),
               gp = gp
            )
         }

         if (grepl("r", sides)) {
            rugs$y_r <- grid::rectGrob(
               x = grid::unit(1, "npc"),
               y = grid::unit(coords$ymax, "native"),
               width = rug_length$min,
               height = grid::unit(coords$ymax - coords$ymin, "native"),
               just = c("right", "top"),
               gp = gp
            )
         }
      }

      # Add text labels if provided in the data
      if (
         "label" %in% colnames(data) && any(nchar(data$label) > 0, na.rm = TRUE)
      ) {
         text_grobs <- lapply(seq_len(nrow(coords)), function(i) {
            if (nchar(coords$label[i]) > 0 && !is.na(coords$label[i])) {
               # Calculate text position based on which side
               if (grepl("t", sides[i])) {
                  # For top margin
                  x_pos <- grid::unit(
                     mean(c(coords$xmin[i], coords$xmax[i])),
                     # "native"
                     "npc"
                  )
                  y_pos <- grid::unit(1, "npc") - rug_length$min / 2
                  angle <- 0
               } else if (grepl("b", sides[i])) {
                  # For bottom margin
                  x_pos <- grid::unit(
                     mean(c(coords$xmin[i], coords$xmax[i])),
                     # "native"
                     "npc"
                  )
                  y_pos <- rug_length$min / 2
                  angle <- 0
               } else if (grepl("r", sides[i])) {
                  # For right margin
                  x_pos <- grid::unit(1, "npc") - rug_length$min / 2
                  y_pos <- grid::unit(
                     mean(c(coords$ymin[i], coords$ymax[i])),
                     # "native"
                     "npc"
                  )
                  angle <- 90
               } else if (grepl("l", sides[i])) {
                  # For left margin
                  x_pos <- rug_length$min / 2
                  y_pos <- grid::unit(
                     mean(c(coords$ymin[i], coords$ymax[i])),
                     # "native"
                     "npc"
                  )
                  angle <- 90
               }

               label_color_i <- coords$label_colour[i]
               if (is.na(label_color_i)) {
                  if (length(label_color) == 0) {
                     label_color <- jamba::setTextContrastColor(
                        ggplot2::alpha(colour=coords$fill[i],
                           alpha=coords$alpha[i])
                     )
                  }
                  label_color_i <- label_color
               }
               use_fontface <- c("plain", "bold", "italic", "bold.italic")[font]

               text_grob <- grid::textGrob(
                  label = coords$label[i],
                  x = x_pos,
                  y = y_pos,
                  rot = angle,
                  just = c("center", "center"),
                  gp = grid::gpar(
                     col = label_color_i,
                     fontsize = label_size,
                     fontface = use_fontface
                  )
               )
            }
         })
         if (length(text_grobs) > 0) {
            rugs <- c(rugs, text_grobs)
         }
      }

      grid::gTree(children = do.call(grid::gList, rugs))
   },
   optional_aes = c("x", "y", "xmin", "xmax", "ymin", "ymax", "label", "label_colour"),
   default_aes = ggplot2::aes(
      colour = NA,
      fill = "grey30",
      linewidth = 0.5,
      linetype = 1,
      alpha = 1,
      label = "",
      label_colour = NA
   ),
   draw_key = ggplot2::draw_key_polygon
)
