
#' Define block arrow coordinates horizontal or vertical direction
#'
#' @returns `list` with x, y `numeric` vectors.
#'
#' @param xmin,xmax,ymin,ymax `numeric` or `unit`
#' @param arrow_head `numeric` fraction of the width to define
#'    the length of the triangular arrow head.
#'    Default 0.7 defines the arrow head length to be 70%
#'    the full width.
#' @param arrow_width `numeric` fraction of the width to define
#'    the arrow stem width, default 0.7 defines the stem as
#'    70% the full width.
#' @param direction `character` indicating direction:
#'    * 'u' Up
#'    * 'd' Down
#'    * 'l' Left
#'    * 'r' Right
#' @param do_plot `logical` whether to plot output using base R
#'    graphics.
#' @param add `logical` whether to add lines to an existing
#'    plot, used only when `do_plot=TRUE`.
#' @param ... additional arguments are passed to `plot()` or
#'    `lines()` when `do_plot=TRUE`.
#'
#' @examples
#' make_blockarrow_coords(xmin=0.92, xmax=1, ymin=0.1, ymax=0.98, direction="u")
#' make_blockarrow_coords(xmin=0.92 - 0.2, xmax=1 - 0.2, ymin=0.1, ymax=0.98, direction="u", arrow_head = 0.7, add=TRUE)
#' make_blockarrow_coords(xmin=0.55, xmax=0.98, ymin=0.92, ymax=1, direction="r", add=TRUE)
#' make_blockarrow_coords(xmin=0.02, xmax=0.45, ymin=0.92, ymax=1, direction="l", add=TRUE)
#' make_blockarrow_coords(xmin=0.92, xmax=1, ymin=0, ymax=0.099, direction="d", add=TRUE)
#' @keywords internal
#' @noRd
make_blockarrow_coords <- function(
   xmin = 0,
   xmax = 1,
   ymin = 0,
   ymax = 1,
   arrow_head = 0.7,
   arrow_width = 0.7,
   direction = c("u"),
   do_plot = FALSE,
   add = FALSE,
   ...
) {
   # assume pointing up:
   # head tip
   # head base
   # stem top
   # stem bottom
   #
   # xmin, xmax
   # ymin, ymax
   direction <- match.arg(
      direction,
      choices = c("u", "d", "l", "r"),
      several.ok = TRUE
   )

   xylen <- max(c(length(xmin), length(xmax), length(ymin), length(ymax)))
   if (length(xmin) < xylen) {
      xmin <- rep(xmin, length.out = xylen)
   }
   if (length(xmax) < xylen) {
      xmax <- rep(xmax, length.out = xylen)
   }
   if (length(ymin) < xylen) {
      ymin <- rep(ymin, length.out = xylen)
   }
   if (length(ymax) < xylen) {
      ymax <- rep(ymax, length.out = xylen)
   }
   if (length(direction) < xylen) {
      direction <- rep(direction, length.out = xylen)
   }
   if (length(arrow_head) < xylen) {
      arrow_head <- rep(arrow_head, length.out = xylen)
   }
   if (length(arrow_width) < xylen) {
      arrow_width <- rep(arrow_width, length.out = xylen)
   }

   # debug
   # print(list(
   #    xmin = xmin,
   #    xmax = xmax,
   #    ymin = ymin,
   #    ymax = ymax,
   #    direction = direction
   # )) # debug

   #
   xypaths <- jamba::rbindList(lapply(seq_len(xylen), function(k){
      if (grepl("[ud]", direction[k])) {
         arrow_width_diff <- (1 - arrow_width[k]) / 2
         xdiff <- arrow_width_diff * (xmax[k] - xmin[k])
         xmid <- (xmin[k] + xmax[k]) / 2
         x1 <- xmin[k] + xdiff
         x2 <- xmax[k] - xdiff
         x <- grid::unit.c(
            xmin[k],
            xmid,
            xmax[k],
            x2,
            x2,
            x1,
            x1,
            xmin[k]
         )
         if (grepl("u", direction[k])) {
            yhead <- ymax[k] - arrow_head[k] * (xmax[k] - xmin[k])
            y <- grid::unit.c(
               yhead,
               ymax[k],
               yhead,
               yhead,
               ymin[k],
               ymin[k],
               yhead,
               yhead
            )
         } else if (grepl("d", direction[k])) {
            yhead <- ymin[k] + arrow_head[k] * (xmax[k] - xmin[k])
            y <- grid::unit.c(
               yhead,
               ymin[k],
               yhead,
               yhead,
               ymax[k],
               ymax[k],
               yhead,
               yhead
            )
         }
      } else {
         #####################
         ## right or left
         arrow_width_diff <- (1 - arrow_width[k]) / 2
         ydiff <- arrow_width_diff * (ymax[k] - ymin[k])
         ymid <- (ymin[k] + ymax[k]) / 2
         y1 <- ymin[k] + ydiff
         y2 <- ymax[k] - ydiff
         y <- grid::unit.c(
            ymin[k],
            ymid,
            ymax[k],
            y2,
            y2,
            y1,
            y1,
            ymin[k]
         )
         if (grepl("r", direction[k])) {
            xhead <- xmax[k] - arrow_head[k] * (ymax[k] - ymin[k])
            x <- grid::unit.c(
               xhead,
               xmax[k],
               xhead,
               xhead,
               xmin[k],
               xmin[k],
               xhead,
               xhead
            )
         } else if (grepl("l", direction[k])) {
            xhead <- xmin[k] + arrow_head[k] * (ymax[k] - ymin[k])
            x <- grid::unit.c(
               xhead,
               xmin[k],
               xhead,
               xhead,
               xmax[k],
               xmax[k],
               xhead,
               xhead
            )
         }
      }
      as.data.frame(list(x=I(x), y=I(y), id=k))
   }))

   if (isTRUE(do_plot)) {
      for (i in seq_len(xylen)) {
         idf <- subset(xypaths, id %in% i)
         x <- idf$x;
         y <- idf$y;
         if (isTRUE(add)) {
            lines(x, y, type = "l", asp = 1, xlim = c(0, 1), ylim = c(0, 1), ...)
         } else {
            plot(x, y, type = "l", asp = 1, xlim = c(0, 1), ylim = c(0, 1), ...)
            add <- TRUE;
         }
      }
   }
   # grid::grid.newpage()
   # grid::grid.path(x=x, y=y)
   return(xypaths)
}
