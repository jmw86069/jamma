# jamma version 0.0.12.900

## changes

* `jammaplot()` argument `useRank=TRUE` will produce the
rough equivalent of a non-parametric MA-plot, where
the y-axis unit is either the actual change in rank,
or the fractional change in rank adjusted by the
number of rows.

# jamma version 0.0.11.900

## new functions

* `outer_legend()` is used to place a color legend outside
multi-panel base R graphics, for example when using
`par("mfrow"=c(2,2))`. By default it places the legend
at the bottom, centered in the plot device. It
of course may move to `jamba` but that package
is already getting overloaded.

## changes

* `jammaplot()` has transitioned from using `groupSuffix`
to add a text suffix to each panel title, to using
`subtitle` to display text in the bottom-left corner
of each panel. In most cases, this change results in
better use of whitespace on each plot, and helps
keep the title label smaller.
* `jammaplot()` allows a list color ramps in the
argument `colramp`, which applies a different color ramp
in each panel.
* `jammaplot()` changed how highlightPoints are handled,
making it more robust to different types of input, vector
or list. Fix error when highlightPoints had no list name.

# jamma version 0.0.10.900

## new functions

* `jammacalc()` is used to calculate data ready for
MA-plots, to separate this logic from the `jammaplot()`
function used for visualization.
* The `jammacalc()` function slightly
improves upon the `jammaplot()` calculations by
simplifying the centering by calling `centerGeneData_new()`
(or `centerGeneData()`) instead of performing similar
calculations itself. The results should be notably
faster especially for larger data matrices.
* The `jammacalc()` function also changes the result
of `groupedX=TRUE`, which now correctly uses each group
median value instead of a median across the group medians.
It allows each group of MA-plot panels defined by
`centerGroups` to be completely independent of other groups,
which has benefits and discrepancies that depend upon how
the `centerGroups` values are defined relative to the intended
downstream statistical comparisons. For example, we recommend
centerGroups broadly include the sample groups that will
eventually be compared, specifically so the distributions are
directly compared in the MA-plot visualization. However,
certain per-group QC checks are ideally performed by using
`centerGroups` for each sample group. Overall, be aware of
the effects of each scenario, and adjust the expectations
and assumptions accordingly.
* `jammaplot()` new argument `ma_method` where
`ma_method="old"` calls the previous code for MA-plot calculations,
and `ma_method="jammacalc"` calls `jammacalc()` with the
newer equivalent calculations.

## changes

* Added `"colorjam"` to R package dependencies.
* Removed R package dependencies on `"rgeos"` and `"sp"`.
* `points2polygonHull()` now uses `grDevices::chull()`
instead of `rgeos::gConvexHull()`. The function accepts
numeric `matrix`, `data.frame`, `tibble`, `tbl`, and `SpatialPoints`
input, and can output either a two-column numeric `matrix`,
or `SpatialPolygons` as output.
* the default argument to `centerGeneData()` was changed from
`"indata"` to `"x"` to be consistent with Tidyverse syntax.

## future changes to centerGeneData()

* `centerGeneData()` will soon be replaced by `centerGeneData_new()`
after a period of testing. The functions should be equivalent.
The new function uses `jamba::rowGroupMeans()` to produce the
summary values used for centering. The change will consolidate the
logic of calculating row group values into `jamba::rowGroupMeans()`,
as opposed to duplicating only certain parts of that functionality
inside `centerGeneData()`. The new function also surfaces outlier
detection and removal during centering in one step.
The potential downside, and need for a period of testing, is that
the full input matrix is centered versus the expanded row group
matrix in one step. Briefly there will be two matrices with identical
dimensions, which will produce a higher memory profile.
For extremely large matrices, this change may become problematic,
and may require an alternative strategy during the centering step.

* `jammaplot()` will be refactored to streamline certain calculations
using the newer `centerGeneData_new()`.

# jamma version 0.0.9.900

## changes to jammaplot()

* now creates a default `groupSuffix` based upon the
`centerGroups` if any, the goal is to ensure the center groups are
by default included in each panel label.
* Columns with missing data (e.g. all NA values) are now handled
by plotting a blank plot panel with background color
`outlierColor`. The MAD value is assigned `Inf`.
* now uses `maintitle` as the title for the overall
set of plots, displayed in the top outer margin of each page. New
argument `maintitleCex` controls title text size.
* the default input argument is changed from `"object"`
to `"x"`, to be consistent with other tidyverse-friendly workflows.
* attribute names are changed from "mvaMADs", "mvaMADoutliers",
"mvaMADfactors" to "MADs", "MADoutliers", "MADfactors". The names
using "mva" are deprecated and will be retired in a future version
of jamma.
* outlier MAD label cex is now `0.9*titleCex`, making the MAD
label slightly smaller than the panel label.
* Outlier plot panels are now highlighted using `>= outlierMAD`,
consistent with the formal outlier detection. Previously it was
erroneously using `> outlierMAD`.
* new argument `subtitleBoxColor` which defines the background color
for optional subtitles printed at the bottom of each plot panel.
Used to color-code `centerGroups`, if `groupSuffix=""` and subtitle
contains `centerGroups` as labels.

## Notes

* Very likely, the base R plotting may be replaced or augmented by
using ggplot2. Past reason for using base R were due to panel layout,
and wanting detailed control over smoothScatter per-panel density
plotting.

* In future, the `centerGroups` value may be displayed inside
the `subtitle` label box at the bottom-left corner of each panel,
color-coded to help indicate the centering group.
* In future, `controlSamples` may be indicated by appending an
asterisk "*" at the end of the main label for each panel.

# jamma version 0.0.8.900

## changes

* `drawLabels()` and `coordPresets()` were moved to the `jamba` package,
consistent with the need for wider re-use among Jam packages.
* Added dependency on the version of jamba that includes `drawLabels()`.

# jamma version 0.0.7.900

## changes

* Based upon feedback regarding difficulty installing the "rgeos" package,
changed `jammaplot()` to use `grDevices::chull()` when
`doHighlightPolygon=TRUE`. The output should be identical, without the
need to install the rgeos package.
* changed `drawLabels()` now uses label height to define the effects of
`boxCexAdjust` on the label box size, making the x- and y-adjustment
more consistent when using long labels; it also uses `jamba::getPlotAspect()`
to correct for aspect ratio. Note it is probably far less useful to supply
x and y values for `boxCexAdjust`.
* made adjustments based upon the single-line height and not the total
string height, affecting multi-line labels.

# jamma version 0.0.6.900

## new functions

* `drawLabels()` which draws a colored box with label inside, positioned
by default at the center of the x and y coordinates provided. This function
will be tested then moved to the "jamba" package.

## changes to existing functions

* `jammaplot()` by default will use `drawLabels()` to print a label in
each panel, instead of the previous method which used `legend()` and was
not properly centered. In future, labels may detect whether to enforce
word-wrap in order to keep labels from overlapping adjacent panels.

## other changes

* pkgdown was added to the site documentation.

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
