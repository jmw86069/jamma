# jamma version 0.0.5.900

## changes to existing functions

* `jammaplot()` had a fairly substantial shift in default behavior,
concerning the x-axis coordinate when using `centerGroups` to
center a subset of samples. Previous default was to use one global
x-axis coordinate across all samples, regardless of the `centerGroups`.
New default behavior is to calculate an x-axis value for each
`centerGroups` so within each group the data is self-contained.
This behavior can be controlled with `groupX=TRUE` which applies
grouping to the x-axis, and `groupX=FALSE` for the previous
behavior. I admit, this change has pros and cons. On the one hand,
there is inherent benefit for keeping each
subset grouping independent, with the possible downside that
comparisons across groupings may not be visible.


# jamma version 0.0.4.900

## enhancements

* Bumped the version number to cover a few smaller updates.

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
