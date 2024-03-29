% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jam-element_textbox_colorsub.R
\name{element_textbox_colorsub}
\alias{element_textbox_colorsub}
\title{Custom ggplot2 element textbox highlight}
\usage{
element_textbox_colorsub(
  ...,
  colorSub = NULL,
  width = grid::unit(1, "npc"),
  padding = ggplot2::margin(6, 2, 4, 2),
  margin = ggplot2::margin(1, 1, 1, 1),
  linetype = 1,
  r = grid::unit(0, "pt"),
  halign = 0.5,
  valign = 0.5
)
}
\description{
Custom ggplot2 element textbox highlight
}
\details{
This function is used internally by \code{ggjammaplot()} to colorize
ggplot2 facet strip panel background using a named color vector
supplied as \code{colorSub}.
}
\seealso{
Other jam ggtext functions: 
\code{\link{element_grob.element_textbox_colorsub}()}
}
\concept{jam ggtext functions}
