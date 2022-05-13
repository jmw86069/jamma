
#' Add legend in the outer margin
#'
#' Add legend in the outer margin
#'
#' This function is intended to place a legend outside the figure
#' in the outer margin area. The motivating use case is to
#' display a legend below multi-panel base R graphics plots,
#' for example when using `par("mfrow"=c(2,2))`.
#'
#' @param x,y position of the legend, using any method
#'    recognized by `grDevices::xy.coords()`.
#' @param legend character vector, or expression, used to
#'    display each legend label.
#' @param doPar logical indicating whether to apply
#'    `graphics::par()` in order to place the legend in
#'    the outer margin. If `doPar=FALSE` and if no other
#'    method has been done to update the `graphics::par()`
#'    settings, the legend will be placed with the same
#'    mechanism used by `graphics::legend()`.
#' @param bg color used as the background of the legend box.
#' @param box.col color used as the frame border color of the
#'    legend box.
#' @param col `character` vector of colors used to color each
#'    entry in `legend`.
#' @param pt.bg vector of colors used as the background color
#'    when `pch` is any value from 21 to 25. When `pt.bg=NULL`
#'    this value and `col` are adjusted so the outline is
#'    darker, via `jamba::makeColorDarker()`, with the fill
#'    color `pt.bg` using the original `col` value.
#' @param pch integer or character vector of point shapes,
#'    as described in `graphics::points()`.
#' @param ncol integer number of columns to use for legend
#'    entries. By default `ncol` is the smaller of
#'    `length(legend)` and `5`, which is intended for a legend
#'    placed at the bottom of the visualization.
#' @param cex character expansion multiplier, used to size
#'    the legend text overall.
#' @param pt.cex point expansion factor, used to size points
#'    when `pch` is not `NULL`.
#' @param ... additional arguments are passed to
#'    `graphics::legend()`.
#'
#' @family jam utility functions
#'
#' @examples
#' opar <- par(no.readonly=TRUE);
#' on.exit(par(opar));
#'
#' par("mfrow"=c(2,2));
#' for (i in 1:4) {
#'    jamba::nullPlot(plotAreaTitle=paste("Panel", i));
#' }
#' require(jamba);
#' outer_legend("bottom",
#'    legend=c("one", "two", "three"),
#'    col=colorjam::rainbowJam(3));
#'
#' @export
outer_legend <- function
(x="bottom",
 y=NULL,
 doPar=TRUE,
 legend,
 bg="white",
 box.col="grey70",
 col,
 pt.bg=NULL,
 pch=21,
 ncol=min(c(5, length(legend))),
 cex=1.2,
 pt.cex=1.4,
 ...)
{
   ##
   if (length(pch) == 0) {
      pch <- 21;
   }
   if (length(pt.cex) == 0) {
      pt.cex <- 1.4;
   }
   pch <- rep(pch,
      length.out=length(legend));
   if (length(pt.bg) == 0) {
      pt.bg <- col;
      col <- ifelse(pch %in% c(21:25),
         jamba::makeColorDarker(col),
         col);
   }

   if (doPar) {
      opar <- par(no.readonly=TRUE);
      on.exit(par(opar));
      par("mfrow"=c(1,1));
      par(fig=c(0, 1, 0, 1),
         oma=c(0, 0, 0, 0),
         mar=c(0, 0, 0, 0),
         new=TRUE);
      jamba::nullPlot(doBoxes=FALSE);
   }

   legend(x=x,
      y=y,
      legend=legend,
      bg=bg,
      col=col,
      pt.bg=pt.bg,
      pch=pch,
      ncol=ncol,
      cex=cex,
      pt.cex=pt.cex,
      box.col=box.col,
      ...);
}
