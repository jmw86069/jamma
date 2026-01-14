# jamma 0.0.38.900

* Removed dependencies: crayon, ggtext

## Updates to existing functions

* `jammaplot()`

   * Fixed bug with `apply_transform_limit` when input data contained NA.
   * Removed deprecated arguments: filterNeg, groupSuffix

# jamma 0.0.37.900

* Added Enhances: 'SummarizedExperiment', 'Biobase', to support additional
input object types, but are not "Suggests" so they avoid some CRAN checks.
* Added Imports: 'cli' for improved commandline messaging.

## changes to existing functions

* `jammaplot()`, `ggjammaplot()`

   * Added support for `ExpressionSet` input class, which includes
   NanoStringGeoMxSet for example.
   * Fixed the warning regarding `pmax(margins, ...)`.
   * Arguments `centerGroups`, `subtitle`, `titleBoxColor`, `subtitleBoxColor`
   now accept an unnamed vector if `length(centerGroups)` equals `ncol(x)`,
   making it easier to use.
   It is preferred to use colnames from `SummarizedExperiment`
   or `ExpressionSet` which pulls from colData/pData.

* `get_se_assaydata()` and `get_se_colData()` now support `ExpressionSet`.


# jamma 0.0.36.900

## changes

* Added 'withr' to dependencies.
* Bumped 'jamba' requirement to CRAN version 1.0.2, woot!
* `jammaplot()`

   * 'SummarizedExperiment' input now allows args `centerGroups`,
   `titleBoxColor`, `subtitle`, `subtitleBoxColor` to match `colData(x)`.
   * ylim,xlim can no longer have length 1 to avoid errors in smooth scatter.
   * Calls to `par()` use `withr::local_par()` or `withr::with_par()`
   to revert changes.

* `ggjammaplot()`

   * 'SummarizedExperiment' input now allows args `centerGroups`,
   `titleBoxColor`, `subtitle`, `subtitleBoxColor` to match `colData(x)`.
   * `subtitle` is displayed consistent with `jammaplot()`, bottom-left
   corner, and not as a 'ggplot2' subtitle.
   * MAD factor is displayed in the bottom-right corner.

# jamma 0.0.35.900

## changes

* `jammaplot()` now uses `jamba::shadowText()` for `displayMAD=TRUE`,
and the control asterisk labels, to improve visibility. It might be
an improvement or distraction.

## bug fixes

* Fixed `jammaplot()` for data with NA values when checking
`apply_transform_limit`. Now `apply_transform_limit` can be NULL,
skipping this step, otherwise it only checks non-NA values.


# jamma 0.0.34.900

* Added dependency `ggh4x` to handle facet fill color. In future, it
may provide improved faceting options, groups, multiple layers, etc.

## bug fixes

* `ggjammaplot()`

   * fixed error with missing `colnames(x)` or `rownames(x)`

## updates

* Bumped jamba version dependency to 0.0.104.900 to support SparseMatrix
objects as with `SingleCellExperient` and `Seurat` matrix data.
* `centerGeneData()`

   * Updated to add `includeAttributes=FALSE` with `jamba::rowGroupMeans()`.
   * Added tests to all basic use cases, and to cover SparseMatrix input.

* `jammaplot()`

   * Removed argument `ma_method`, it was deprecated.
   * Added argument `apply_transform_limit=40` to apply `jamba::log2signed()`
   when any value is above this threshold, only when `useRank=FALSE`.

* `ggjammaplot()`

   * Added argument `apply_transform_limit=40` to apply `jamba::log2signed()`
   when any value is above this threshold, only when `useRank=FALSE`.
   * Now uses `ggh4x::facet_wrap2()` to apply color to facet strips.

## internal functions removed

* `element_textbox_colorsub()`, `element_grob.element_textbox_colorsub()`
were removed, in favor of using `ggh4x` functions for ggplot2 facet strip
color fill.

# jamma 0.0.33.900

## changes to existing functions

* `jammaplot()`

   * A warning was shown `"Warning in min(x, na.rm = na.rm) :
  no non-missing arguments to min; returning Inf"` which was caused by
  `xlim` not being defined for each plot panel, thereby causing
  `ablinesV` to throw an error because no values fit within the undefined
  x-axis range. Now when `xlim` argument is not provided, it uses
  the range of x-axis values from `jammacalc()` for each panel.
  * This change makes panel x-axis ranges more consistent by default.
  They should have been consistent already, but in some cases `NA` values
  caused the range to differ slightly from panel to panel. In most cases,
  it caused no visible effect in x-axis labels.
  * This change also allows samples within `centerGroups` groupings to
  have distinct x-axis ranges, since each individual group can represent
  potentially very different values depending upon the experiment.

# jamma 0.0.32.900

## bug fixes

All instances of `if (class(x) %in% c("a", "b"))` were corrected to use
format compatible with R version 4 which considers this syntax an error
and not a warning. New format is `if (any(c("a", "b") %in% class(x)))`.

* `jammaplot()`, `update_list_elements()`, `centerGeneData()`, `ggjammaplot()`

   * Fixed `if` statement logic to remove length=2 error in all scenarios.

## other changes

* `jammaplot()`

   * argument help docs were updated for clarity.
   * change argument default to `filterNeg=FALSE`, although this argument
   is deprecated and only used when `ma_method="old"`. The change was
   made to reflect the current default behavior which uses `noise_floor=0`
   and `noise_floor_value=NA`, which replaces all values at or below 0
   with `NA`, causing these values not to contribute to the data
   centering and MA-plot output.


# jamma 0.0.31.900

## changes to existing functions

* `jammaplot()`

   * New argument `outer_margins` to control the whitespace around
   the outside of the MA-plot panels, only when `doPar=TRUE`.
   * The y-axis labels are only shown on the first plot each row,
   and plot panels are condensed side-by-side.
   * The default `margins` has much smaller `left` value, so plots
   can be placed closer side-by-side.
   * Far very many MA-plot panels, the plot panels are much more visibly
   condensed, especially when `useRank=TRUE`, which previously increased
   each panel `left` margin by 2 lines; now it only increases the `left`
   outer margin.


# jamma 0.0.30.900

## changes to existing functions

* `jammacalc()`

   * now inherits arguments `controlFloor`, `naControlAction`, `naControlFloor`
   consistent with `centerGeneData()`.
   * now calls `matrix_to_column_rank()` instead of calculating internally.

* `jammaplot()` and `ggjammaplot()`

   * gained arguments `controlFloor`, `naControlAction`, `naControlFloor`
   consistent with `centerGeneData()`. These arguments are passed
   to `jammacalc()`.

# jamma 0.0.29.900

## new functions

* `matrix_to_column_rank()`

   * converts `numeric` matrix to column rank data, as used when
   `useRank=TRUE`.
   * This function becomes its own standalone function so it can be
   called when necessary in other areas.

## updates to existing functions

* `jammacalc()`

   * calls `matrix_to_column_rank()` instead of calculating internally.

# jamma 0.0.28.900

## updates to existing functions

* `centerGeneData()`

   * Note that default behavior is unchanged.
   * New argument `controlFloor` to impose a minimum control summary value
   when the summary value is below `controlFloor`.
   * New argument `naControlAction` used when all control sample values
   are `NA`:
   
      1. `naControlAction="na"`: leave values centered values as `NA`
      (default behavior)
      2. `naControlAction="row"`: center versus remaining non-NA row values
      3. `naControlAction="floor"`: center versus numeric floor
      4. `naControlAction="min"`: center versus the minimum observed value
   
   * New argument `naControlFloor=0` used only when `naControlAction="floor"`.
   When zero, values are effectively un-centered. However, sometimes
   the range of signal is 20-35, in which case it may be practically more
   useful to use `naControlFloor=20` so data is centered relative to
   the minimum range of detection.

* Added proper prefix to `ggplot2::after_stat()` for `ggjammaplot()`.

## testthis unit tests

* Began the process of covering package functions with unit tests.
Currently: `centerGeneData()`.


# jamma 0.0.27.900

## updates to existing functions

* `ggjammaplot()` was updated:

   * Point density now shares the same range for all panels,
   including outlier panels. Previously, outlier panels generated
   their own color gradient, due to their specific numeric range
   of point density values calculated by ggplot on the fly.
   The new code will define the numeric range using non-outlier
   samples, then re-uses the same numeric range for outlier panels
   when there are outlier panels. This step appears to add some lag
   to rendering, although there is already more lag than with
   base R graphics.
   * New argument `detail_factor` is a single argument to adjust the
   level of detail in the point density per facet panel.

* `jammacalc()`

   * The calculation of rank data using `rank(..., na.last=FALSE)` was
   assigning ties to repeated `numeric` values, but was assigning
   unique (arbitrary, by order they occur) ranks to `NA` values.
   The new code secretly converts `NA` to the lowest `numeric` value
   so they are assigned tied rank values.

# jamma 0.0.26.900

* `volcano_plot()` updates:

   * Fixed bug where `lfc_colname` was used even when only `fold_colname`
   was appropriate for the data. It now properly uses the `lfc_values`
   internal vector for relevant steps.
   * `hi_points` now accepts a `list` input, to color code different subsets
   of points different colors
   * new argument `hi_colors` coincident with `hi_points` having `list`
   input: `length(hi_colors)` should match `length(hi_points)` so colors
   will be applied to each vector in the `list`. If no colors are assigned,
   then bright categorical colors will be defined by `colorjam::group2colors()`.
   * when `hi_points` is a `list`, a color key is placed at the bottom
   of the plot, making the bottom margin 2 lines taller, and shifting other
   captions up by 2 lines.

# jamma 0.0.25.900


Main change: this version:

1. Automatic indication of `controlSamples` in each MA-plot panel.
2. Ability to supply custom names for each MA-plot panel.

## updates to existing functions

* `jammaplot()`

   * `whichSamples` accepts `character` or `numeric` and will subset
   itself to match available `names(x)` or `nsamples` accordingly.
   * New argument `sample_labels` to provide custom labels for each
   sample without having to rename `colnames(x)` in the process.
   An example has been added.
   * New argument `controlIndicator`:
   
      * `"titlestar"` - asterisk appended to the title of each plot
      * `"labelstar"` - asterisk labeled inside the top corner of each plot
      * `"none"` - no indicator is drawn
   
   * `par("mar")` is now only updated when `doPar=TRUE`.
   * New argument `panel_hook_function` to supply a custom function to be
   run after drawing each plot panel.
   An example has been added, showing how to draw a box around each figure.
   * y-axis labels are drawn using `format(..., big.mark=",")` as
   an improvement for `useRank=TRUE`.
   * Margins have been adjusted so the title box default above each plot
   has more room when there is no `maintitle` that would otherwise
   have added space to the top of each overall figure.

* `find_colname()` is now exported.
* added `log2fold_to_fold()` and `fold_to_log2fold()` functions as required
by `volcano_plot()`.

# jamma 0.0.24.900

## bug fixes

* `volcano_plot()` was throwing an error due to the new use of
`utils::modifyList()`, ultimately caused by `update_function_params()`.
That function was updated, resolving the error.
* `jammaplot()` was updated to pass `panelWidth="minimum"` to
`jamba::drawLabels()` for the title box when using base R
graphics plots. This change ensures the box is at least the plot panel
width, and will be larger if the label is wider than each panel.


# jamma 0.0.23.900

## changes to existing functions

* `jammaplot()` title box labels are now placed *above* plot panels!

   * This changes is being tested, but it seems like a sensible default
   not to obscure points in each plot panel.
   * two new arguments: `titleAdjPreset`, `subtitleAdjPreset`
   * passed to `jamba::coordPresets()` which defines coordinate positions
   in base R graphics context.
   The two main arguments are `preset` and `adjPreset`, where `preset`
   defines the coordinate "edge" or "corner" position, typically something
   like "top" or "bottomleft"; and `adjPreset` defines relative placement
   of the label. By default, `adjPreset` chooses the opposite orientation,
   to keep labels inside the plot panel. However `preset="top"` defines
   the top edge of the plot, and `adjPreset="top"` places the label
   above the top edge of the plot panel.
   * The calculation of top outer margin is dependent upon the number of
   text lines in `maintitle`, and was increased to accomodate the new
   default placement of title box labels above the top panel. It does
   not account for multi-line title box labels, in that case, supply
   `maintitle` with empty extra lines at the end to add spacing.

* `jammaplot()` argument `margins=c(3.5, 2, 0.3, 0.2)` is typically good
for non-ranked MA-plots. When `useRank=TRUE` the numbers typically have
4 or more digits, and need more space. In this case by default,
`margins + c(0, 1, 0, 0)` is used.

## bug fixes / improvements

* `jammaplot()` correctly represents `highlightPch` point shapes in
the highlight point legend at the bottom, via `outer_legend()`.
* `outer_legend()` when `pt.bg=NULL` it will adjust `col` slightly
darker for point shapes ranging from 21 to 25, and will assign `pt.bg=col`
so the fill color is correct.


# jamma 0.0.22.900

## new functions: volcano_plot()

The `volcano_plot()` function has been in use for a long time,
but not in R package form. The sticking point has been the
block arrows drawn in plot margins (outside the plot panel)
which describe the number of points that exceed statistical
thresholds.
Colleagues seem to like having these block arrows
with summary numbers, and they've appeared in some published papers.
That said, sizing block arrows consistently for all plot aspect
ratios has been tricky, and migrating everything to ggplot2
is on the short to-do list. Unclear how block arrows would work
in ggplot2.

* `volcano_plot()` - Draw volcano plot, work in progress
* `blockArrowMargin()` - draws block arrows in plot margins
* `gradient_rect()` - draw rectangle with gradient color fill
* `logAxis()` - draw log-transformed `-log10(pvalue)` axis labels
* `find_colname()` - find colname matching a character string
* `update_function_params()` - update function arguments in place
* `update_list_elements()` - utility called by `update_function_params()`.

## immediate change to existing functions

* `update_list_elements()` appears to have identical purpose and
functionality as `utils::modifyList()`

   * `update_function_params()` was updated to call `utils::modifyList()`
   * `update_list_elements()` is deprecated while this change is tested.


# jamma 0.0.21.900

## notable changes to default values

`jammaplot()` argument default `filterFloorReplacement=NA`:

* This change corrects an issue where the row mean was calculated
using values that may contain a large proportion of zero `0`,
therefore pushing the x-axis position toward zero. The effect
often also pushed the y-axis difference-from-mean to slightly
larger magnitude, since the baseline value was typically closer
to zero.
* The effect of the change in argument defaults is that the x-axis
position of points displayed will reflect the mean/median of
values above zero, and y-axis will reflect the difference from
mean/median of values above zero, thus displaying mean and variation
of those measurements that detected a signal, without being influenced
by the proportion of those points that detected a signal.
* Specifics:

   * `noise_floor` and `noise_floor_value` will replace deprecated
   arguments `filterFloor` and `filterFloorReplacement`, respectively.
   But the main change is to the default values.
   * previous default: `filterFloor=0` and `filterFloorReplacement=filterFloor`
   set any value at or below `0` to `0`.
   * new default: `noise_floor=0` and `noise_floor_value=NA`
   sets any value at or below `0` to `NA`.
   * The change mostly affects data with large number of zero `0` values.
   * Points below `noise_floor` will no longer contribute toward
   the sample MAD factor, which should be correct as the default.
   * The x-axis summary value, by default is the mean unless `useMedian=TRUE`
   in which case it uses the median, is calculated using the remaining
   measurements, above the `noise_floor`, so by default values above zero
   will be displayed in each plot panel.
   * The assumption is that a measurement `0` is an absence of measurement,
   and therefore should not contribute to the mean for that row of data,
   nor to the variability (or absence of variability) for the sample column.
   * The same assumption should hold true for values below zero, where
   upstream processing should in theory have generated negative values
   as a result of background subtraction/adjustment by appropriate methods.
   Negative values are therefore also considered "absence of measurement"
   and the specific numeric value has no quantitative meaning in this
   function.
   
* All that said, the previous behavior can be used with argument
`noise_floor_value=noise_floor`.
* Finally, some case can be made that the points filtered may optionally
be displayed, in order to indicate their presence in the data. This option
may be enabled in future if it seems necessary.

   * It is unclear how to plot points whose values were zero.
   * If for example the non-zero mean x-axis position for a row is 12,
   does it imply the points at for this row would use x=12 and y=-12?
   Or should these points use x=0 and y=0?
   * It seems more sensible not to display points whose measurements have been
   filtered from analysis - at least for now.


## changes to make `jammaplot()` and `ggjammaplot()` consistent:

* `ggjammaplot()` was modified for consistency with `jammaplot()`:

   * `colramp` argument added
   * `outlierColor` can be a color function, for example `colramp="inferno"`
   and `outlierColor="inferno_r"` will use the reverse color ramp for
   a nice visual effect.
   * new argument `fillBackground=TRUE` optionally fills each
   plot panel with the first color ramp color, affecting outlier panels.
   Previously each plot panel has a small region around the outside
   that was not filled, due to ggplot2 axis range expansion slightly
   beyond the 2d density coordinates.
   * MAD outlier panels were not properly recognized, the MAD values
   were displayed but the background was not highlighted with `outlierColor`.
   Issue has been resolved.
   * `noise_floor=0` is the new default value, consistent with `jammaplot()`.
   Any value **at or below** is set to `noise_floor_value` which by default
   is NA, but can be changed to `noise_floor`.
   * The MAD factor text color uses a color that contrasts with the base
   color gradient used in each plot panel, using white text versus
   dark background for example. It calls `jamba::setTextContrastColor()`.
   Outlier text is either red or gold, as appropriate.

* `jammaplot()` was modified for consistency with `ggjammaplot()`:

   * New arguments `noise_floor` and `noise_floor_value`, which replace
   * Deprecated arguments `filterFloor` and `filterFloorValue`.
   * `filterFloorReplacement=NA` is the new default value, which may
   be a substantial change for data that contains a large proportion
   of missing values.
   * The MAD factor text color uses a contrasting color, as with
   `ggjammaplot()`.
   * `outlierColor="lemonchiffon"` default, consistent with `ggjammaplot()`.
   It sounds yummy.


# jamma 0.0.20.900

## bug fixes

* `jammaplot()` was not properly passing argument `assay_name` to
downstream function `get_se_assaydata()` when the input was
`SummarizedExperiment`. It was by default always using the
first assay stored in the SE object.


# jamma 0.0.19.900

## updates to package dependencies

Unfortunately the `gridtext` R package issue
https://github.com/wilkelab/gridtext/issues/22 results in removing
this package and `ggtext` as dependencies from `jamma`. The installation
requires specific recent versions of GCC compiler that I have trouble
resolving on my own linux servers. As a result, they are moved to
`Suggests` and the corresponding `jamma` functions will handle
them as optional.


# jamma 0.0.18.900

Bumped dependency on `jamba` to version `0.0.66.900` to
avoid error with `rowMedians()` that should use `matrixStats::rowMedians()`.

Note: The `pkgdown::build_site()` failed due to SSL certificate expiration,
apparently from the crandb.r-pkg.org site. The workaround is to set this
option beforehand, which prevents `pkgdown::build_news()` from trying
to access that URL.
`options("pkgdown.internet"=FALSE)`

## new functions

* `ggjammaplot()` - the ggplot2 equivalent to `jammaplot()`. See examples,
there are a variety of example figures, using straight MA-plots, raw
and normalized data, using subset of samples, and optionally highlighted
points.

   Still to-do:
   * when `blankPlotPos` is used, the plot itself should be empty, without grid lines
   * add `subtitle` to the bottom-left corner of each plot panel
   * use `geom_text_repel()` to label `highlightPoints` - should probably be optional.
   * more testing for large number of plot panels, for example test the
   aspect ratio for density and pixel size calculations

* `element_textbox_colorsub()` - a custom ggplot2 element that
enables colored facet strip background colors using name-value
pairs.

   * The technique was based upon a post by Claus O. Wilke,
   and is the closest to a "true ggplot2" methodology that I found.
   Other alternatives involved using the `gtable` package, and
   modifying `grob` grid graphical objects directly.
   See stackoverflow questions/60332202/conditionally-fill-ggtext-text-boxes-in-facet-wrap
   * Note that this function very likely will move into the `colorjam`
   package to become the default of `colorjam::theme_jam()`
   * This function has not been tested with `facet_grid()` but in
   principle it should work.

* `element_grob.element_textbox_colorsub()` - required for the ggplot2
workaround with facet strip label background colors.


## new dependency

* `ggtext` package is a dependency for `ggjammaplot()`, however it
is a lightweight and useful package, whose dependencies are largely
also in common with `ggplot2`.



# jamma 0.0.17.900

## changes to existing functions

* `jammaplot()`, `jammanorm()` and `jammacalc()` were
updated to handle `useMedian` by default, but to accept
`useMean` for backward compatibility.


# jamma 0.0.16.900

## changes to existing functions

* `centerGeneData()` has been changed to what was previously
`centerGeneData_new()`, and the old function was renamed
to `centerGeneData_v1()`. The new function is designed to
be backward compatible.
* `jammaplot()` now accepts `SummarizedExperiment` objects
as input, used alongside new argument `assay_name` which
determines the matrix to use for MA-plots. The `assay_name`
input can be an integer index, or by default it will use
the first assay matrix in `x`. Alternatively to use the
last assay matrix (if normalized data is stored in the
last assay position for example) supply `assay_name=Inf`
and the default will use the last assay in the list.
* `jammacalc()` now calls `centerGeneData()`.
* `jammanorm()` help text was expanded.

Added various entries to `TODO.md` for future work.


# jamma 0.0.15.900

## bug fixes

* `jammaplot()` fixed issue when applying `filterFloor`
and `useRank=TRUE`, it was converting to rank before
applying filter floor. An edge case but important for
count data, where using `filterFloor=0` and
`filterFloorReplacement=NA` is useful, and using
`useRank=TRUE` is informative especially when using
methods like DESeq2/edgeR downstream with count data.

## minor updates

* `jammacalc()` new argument `useRank` so the conversion
to rank happens inside this function, as it should.
* Added some package prefix to functions such as
`jamba::nameVector()`, `jamba::rmNA().
* Removed `require(jamba)` and `require(matrixStats)` in
favor of using proper package prefixing, and check for
`matrixStats` using `rownames(installed.packages())`.
* Cleaned up some R code, removing sections that were deprecated.
* Edited the wording for `jammaplot()` function parameters.


## changes to existing functions

* `jammaplot()` cleaned up logic on handling `colramp` and
`colrampOutlier` arguments.
* `jammaplot()` default `ma_method="jammacalc"` which is equivalent
but slightly streamlines the internal calculations.


# jamma 0.0.14.900

## new function

* `jammanorm()` a basic normalization method that uses the output
of `jammaplot()` or `jammacalc()`, and normalizes data such that
the controlGenes are centered at y=0. In words, it makes the
controlGenes have mean zero change from average across samples.
It optionally takes pre-defined controlGenes (housekeeper genes),
and optionally applies an expression minimum threshold, to ensure
genes are expressed above noise.


# jamma 0.0.13.900

## changes

* `jammaplot()` argument default changed to `useMean=TRUE`
to plot the mean centered value, instead of the median.


# jamma 0.0.12.900

## changes

* `jammaplot()` argument `useRank=TRUE` will produce the
rough equivalent of a non-parametric MA-plot, where
the y-axis unit is either the actual change in rank,
or the fractional change in rank adjusted by the
number of rows.

# jamma 0.0.11.900

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

# jamma 0.0.10.900

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

# jamma 0.0.9.900

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

# jamma 0.0.8.900

## changes

* `drawLabels()` and `coordPresets()` were moved to the `jamba` package,
consistent with the need for wider re-use among Jam packages.
* Added dependency on the version of jamba that includes `drawLabels()`.

# jamma 0.0.7.900

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

# jamma 0.0.6.900

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

# jamma 0.0.5.900

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


# jamma 0.0.4.900

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

# jamma 0.0.2.900

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
