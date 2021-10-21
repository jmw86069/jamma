
#' Custom ggplot2 element textbox highlight
#'
#' Custom ggplot2 element textbox highlight
#'
#' This function is used internally by `ggjammaplot()` to colorize
#' ggplot2 facet strip panel background using a named color vector
#' supplied as `colorSub`.
#'
#' @family jam utility functions
#'
#' @export
element_textbox_colorsub <- function
(...,
   colorSub=NULL,
   width=grid::unit(1, "npc"),
   padding=ggplot2::margin(6, 2, 4, 2),
   margin=ggplot2::margin(1, 1, 1, 1),
   linetype=1,
   r=grid::unit(0, "pt"),
   halign=0.5,
   valign=0.5)
{
   #jamba::printDebug("element_textbox_colorsub(): ",
   #   "box.colour: ", paste0("'", box.colour, "'"));
   structure(
      c(ggtext::element_textbox(...,
         width=width,
         padding=padding,
         margin=margin,
         halign=halign,
         valign=valign,
         linetype=linetype,
         r=r),
         list(
            colorSub=colorSub)
      ),
      class=c("element_textbox_colorsub",
         "element_textbox",
         "element_text",
         "element")
   )
}

#' Custom ggplot2 element textbox highlight grob
#'
#' @family jam utility functions
#'
#' @export
element_grob.element_textbox_colorsub <- function
(element,
   label="",
   ...)
{
   #
   if (label %in% names(element$colorSub)) {
      if ("transparent" %in% element$colorSub[label] ||
            jamba::col2alpha(element$colorSub[label]) == 0) {
         element$fill <- "transparent";
         element$colour <- "transparent";
         element$box.colour <- "transparent";
      } else {
         element$fill <- element$colorSub[label];
         element$colour <- jamba::setTextContrastColor(element$fill);
         #element$box.colour <- "#000000";
      }
   } else {
      element$fill <- "transparent";
      element$colour <- "transparent";
      element$box.colour <- "transparent";
   }
   NextMethod()
}

