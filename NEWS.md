# jamma version 0.0.3.900

## enhancements

* `jammaplot()` now more consistently handles highlight
point parameters, where each are expected in list form:
`highlightPoints`, `highlightPch`, `highlightCex`, and
`highlightCex`. If `highlightPoints` is sent as a vector,
it is converted to a one-vector list. If `highlightPoints`
is a list, then the other parameters are converted to list
as needed (`as.list()`) and expanded to `length(highlightPoints)`.
This way, each highlight set of points can have custom
color, size, and shape, to help visually discern them.

# jamma version 0.0.2.900

## new functions

* `centerGeneData()` the main tool for calculations by `jammaplot()`,
centers data by subtracting the row mean or median; optionally within
groups of columns; optionally with specified control columns per
group.
* `jammaplot()` the main MA-plot function.
* `points2polygonHull()` takes a set of points and draws a polygon
hull encompassing them. It calls `rgeos::gConvexHull()`. The function
is intended to indicate the location of highlighted points on each
MA-plot panel.

## modifications

* several default settings were changed for `jammaplot()`: larger
`titleCex`; outlierColramp is NULL, in favor of using outlierColor
to change only the first color from `colramp`.
