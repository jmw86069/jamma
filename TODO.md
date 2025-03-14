# TODO for jamma

## 14mar2025

* Add visual tests, as used in venndir.
* Consider multi-page output with `ggjammaplot()`.
Similarly, consider padding empty cells when `ncol` and `nrow`
are provided but do not fill all available rows/columns.
* DONE. Consider improving visual detail with ggplot2 smooth scatter output.
* Add control asterisk to `ggjammaplot()`.

## 07mar2025

* DONE. Fix bug in `jammaplot()` when NA values are present,
causing an error with `any(x > apply_transform_limit)`.

## 21feb2025

* Prepare for publishing a brief report.

   * Useful references on MA-plots:
   
      * Futshick et al 2018
      https://www.sciencedirect.com/science/article/pii/B9780128027615000079
      Molecular Pathology textbook
      https://doi.org/10.1016/B978-0-12-802761-5.00007-9

## 03feb2025

* Better idea for categorical coloring to indicate factor values

   * Allow one or more `colnames(colData(se))` to be used with
   `sample_color_list` which will add one tile per factor value
   in the corner of each panel. Potential replacement for `subtitle`
   and `subtitleBoxColor`.

## 24jan2025

* DONE. Display MAD factor using `shadowText()` or equivalent.
* DONE. Display control asterisk using `shadowText()` or equivalent.

## 07oct2024

* DONE. Fix `ggjammaplot()` bug where ggplot2 facet strip background color
no longer applies `titleBoxColor` properly. Unclear when this stopped
working, probably something in the ggplot2 3.5.x update.

   * Fix was to use `ggh4x::facet_wrap2()` to apply color and fill to
   ggplot facet strip rectangles.

* Enable `subtitle` arguments for `ggjammaplot()` consistent with
`jammaplot()`.


## 17sep2024

* Add support for `SingleCellExperiment`/`Seurat` objects.

   * Todo: Ensure input to `jammaplot()` and `ggjammaplot()`
   work for Seurat objects.
   * Todo: Ensure all steps are properly compatible with alternate
   matrix objects, mostly from Matrix package.
   * Adopt similar mechanism used in `jamses::heatmap_se()`
   * Use `rowRanges()` when necessary.

* Update `centerGeneData()` to handle `Matrix` objects, not just `matrix`

   * for example, SingleCellExperiment data may provide `dgCMatrix`
   * `inherits(x, "sparseMatrix")` is valid for Matrix sparse matrices,
   and `sparseMatrixStats` may provide equivalent methods to `matrixStats`
   if installed.
   * Check if matrix is sparse, check if `sparseMatrixStats` is installed,
   if so use `sparseMatrixStats::rowMedians()`, otherwise `apply(x, 1, median)`

* Consider `subtitle` ability to use `colData()` colnames.
* Consider `colorSub` to define colors used by arguments
`titleBoxColor` and `subtitleBoxColor`.
* Add asterisk to corner of `ggjammaplot()` output

## 14aug2024

* Consider other data centering options:

   * Use case:
   
      * Compare GroupA to GroupB and GroupC.
      * Goal is to show up-regulated genes in GroupA.
      * Center versus GroupA shows everything down in GroupB/GroupC.
      * Centering by GroupB/GroupC uses the mean of GroupB/GroupC.
      * Proposal is to center by group "min" or group "max" instead
      of group "mean".
      * Possible to accomplish with `rowStatsFunc` argument?
      * If the goal is a heatmap, `jamses::heatmap_se()` could accomplish
      it entirely with custom code:
      
         * Define custom `rowStatsFunc` function that recognized colnames,
         assigned them to groups, returned the appropriate group stat
         (min, max, mean, median).
         * It also needs to label the centering accordingly.
         Entirely specific to `heatmap_se()`. Or is it?

   * Consider including attributes which describe the centering.
   
      * Method: mean/median
      * controlSamples
      * centerGroups
      * Optional labels? `control_label`, `centerby_label`

   * Consider method which accepts `SummarizedExperiment` data, and arguments:
   
      * `assay_name` - if one value, one `matrix` is returned, otherwise `list`
      * `centerby_colnames` - to define `centerGroups`
      * `controlSamples` - same as usual, to define controls for centering
      * Output includes attribures including consistent centering label,
      to be used by MA-plots, `jamses::heatmap_se()` and other tools.


## 09jul2024

* Consider RLE plot - which is a flattened form of MA-plot

   * Motivation: single cell data with thousands of columns, they should
   not become independent panels.
   * Could be useful for bulk 'omics data by providing a quick global view.
   * "Problem samples" could be expanded to show the full MA-plot.
   * For reference see https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html
   the `plotRLExpr()` function shows example before and after normalization.
   * Weakness of this style is that it doesn't show non-linear MA-plot,
   however the box-whisker/median with range is effective at showing
   when one sample has much different variability than others.

* Consider adding signed-significance plot, perhaps `ssplot()`

   * Goal is to provide two `data.frame` objects that have identifiers (gene?)
   that can be aligned to one another.
   * Each data has significance column (adjusted P-value, P-value, Q-value, FDR)
   and a directional column (fold change, log fold change, ratio, etc.).
   * The `-log10(significance) * sign(fold)` is used for each axis.
   * Statistical thresholds are drawn for each dataset, usually P-value
   threshold. Shown as dotted/dashed lines on respective axes.
   * Consider option to impose fold change threshold for each data also?
   * Consider option to use log fold change on each axis.

## 30may2024

* `volcano_plot()`

   * return `data.frame` equivalent to the data plotted, make it easy
   to find hits, and highlighted points.

## 13mar2024

* Consider how to utilize `colData()` sample annotations.

   * Go with "easiest thing that works" since in most cases, MA-plots
   are data-agnostic QC for technical quality assessment.
   However, centering within known sample annotations could be convenient.
   * For example: colorized `subtitle` boxes, or augmented labeling.
   * Situation is that `SummarizedExperiment` input offers more sample
   annotations than are easily added via `subtitle`.
   * Could mimic arguments in `jamses::heatmap_se()`: `centerby_colnames`,
   `colData_colnames`.
   * Could apply sort order to samples by argument `colData_byCols`,
   pass to: `mixedSortDF(colData(se), byCols=colData_byCols)`
   * New argument `sample_color_list`, optional color assignment for `colData`
   columns.

## 27sep2023

* develop better method to organize MA-plot panels

   * for example columns organized by batch
   * generally some assignment of sample to panel number

## 01aug2023

* `centerGeneData()`

   * Enable `SummarizedExperiment` input, either as a separate function,
   or as option for this function. It should require `assay_name`,
   and return a `numeric` matrix.

* Figure out a reasonable way to include experiment design factors in
`jammaplot()` and `ggjammaplot()` panel labels.

   * For example `"BAH013760_101_NEG_MS1"` is not a helpful identifier,
   it would be much more interesting to include `"Male, Dex-Treated, Time 3"`.
   * Bonus points for including each factor on its own line, using
   categorical colors for each value to make it easier to scan by eye.

* Fix missing strip text colors in `ggjammaplot()` when multiple colors
are supplied for `titleBoxColor`.

   * Somehow the ggplot2 `element_grob()` dispatch is not calling
   the custom functon `element_grob.element_textbox_colorsub()`.

## 31jul2023

* DONE. Fix error `ggjammaplot()` when input data `x` does not
contain colnames, or rownames. Either case produces an error.


## 19jul2023

* Some method to detect non-linear (non-horizontal) MA-plot signal

   * Example case is one sample whose values are nearly all zero,
   the MA-plot looks like a diagonal line aiming down from left to right.
   After log-ratio normalization (for example with `jammanorm()`),
   the signal is shift up so the mean signal is at y=0, however
   the diagonal is still there.
   * Simplest possible idea is to fit a line to points above
   the minimum threshold defined by `mad_row_min`, and test the
   slope.
   * It is unclear how much to trust a linear fit, and assuming the linear
   fit is reasonably good, what to do with the slope afterward.
   
      * How much deviation from slope=0 is tolerable?
      * Should the range of acceptable slopes be data-defined?
   
   * Another idea: Could `jammanorm()` only normalize using rows
   where the **raw** values are above the noise threshold `mad_row_min`?
   Currently the process uses the x-axis calculated value (row mean/median)
   to apply that threshold. Logically it seems reasonable that
   measurements below a "noise floor" are noise, and therefore
   should not be considered during normalization.
   
      * For samples where signal is completely "below noise"
      it would result in data not being normalized, which is
      preferable since there is no reliable signal. In this
      case the usual MA outlier strategy would be effective.
      * For samples with reasonable signal with slope deviating from zero,
      it should have no effect (no negative effect).

## 27jun2023

* `jammaplot()` and `ggjammaplot()` - controlSamples

   * improve the indication that a sample (panel) was used as a control
   during data centering. Currently displays an asterisk "*" in the top-left,
   however it is not visible in all plots, especially when there are outlier
   samples. Also, there is no legend indicating the meaning of the asterisk.
   * Consider using asterisk only for samples excluded during centering,
   therefore by default no asterisk is necessary.
   * Consider using `ComplexHeatmap::Legend()` to create a custom legend.
   It would be used to display `"* - samples excluded during centering"`
   and optional highlight points.
   * Consider argument to display/hide the color legend, useful when the
   color legend might be too large to display comfortably.
   Instead provide method for users to create their own legend.

* `jammaplot()` and `ggjammaplot()` - MAD factor label box

   * Consider displaying the MAD factor inside a label box, to improve
   legibility when the text overlaps the point density. This change
   could be optional, controlled with a new argument.
   * Consider using `grid` functions for `drawLabels()` to allow much more
   control over label placement, for example avoiding overlaps
   between `subtitle` and `MAD factor` labels by shifting one label by
   the height of the other label.
   This suggestion might affect `jamba::drawLabels()`.

## 30may2023

* `volcano_plot()`

   * method to highlight points (rows) and assign colors should
   be consistent with `jammaplot()`, with similar color legend drawn.
   * consider colorized smooth scatter in panels that match the
   statistical thresholds used.
   
      * Idea is to render smooth scatter in panels, where color density
      is restricted to points inside these regions, then render these
      panels exactly at these borders. Effectively like running the plot
      three times (red/blue/grey) then copy/pasting only the relevant
      regions for use in the final plot. In fact, that could be one
      potential implementation strategy. It could be faster and more
      reliable than subsetting points, where density at the borders
      might be mis-calculated.
      * top-left: down-significant (blue)
      * top-right: up-significant (red)
      * top-center: no change, met significance (grey)
      * bottom-left, bottom-center, bottom-right: no statistical change (grey)
   
   * Consider option to use scatter points instead of smooth scatter,
   with colorized points in each relevant region.
   
      * Primary driver would be for volcano plots with many relatively
      few points to display, the density plot would be too obscure.
      E.g. NanoString, SomaLogic, with 100 to 1,000 points.


## 25may2023

* Adjust outer margin for base plots

   * display y-axis labels only on the first plot each row
   * shrink the margin between plot panels

## 16may2023

* `volcano_plot()`

   * allow `SummarizedExperiment` input
   * allow optional labeling of genes
   * allow optional `ggplot2` output, which would make labeling genes
   much more feasible.

## 22mar2023

* Empty control group updates

   * Date centering by `centerGeneData()` should be mirrored in MA-plots
   for the same situation where control samples are entirely NA,
   causing remaining values to become `NA` and therefore not be visible
   on MA-plot panels.
   


## 12mar2023

* `centerGeneData()`

   * situation occurs when centering versus `controlSamples` results in
   rows where all `controlSamples` have `NA` values, thereby causing
   all centered values to become `NA`.
   * goal is to offer alternatives where appropriate:
   
      1. `"na"`: leave values NA (current behavior, default)
      2. `"row"`: center versus row non-NA values
      3. `"floor"`: center versus numeric floor
      4. `"min"`: center versus the minimum observed value


## 22nov2022

* `jammaplot_se()` - customized `jammaplot()` for SummarizedExperiment input

   * analogous to `jamses::heatmap_se()` so it should share argument style
   * `normgroup` default uses column(s): `colData(SE)[[normgroup]]
   * `centerby_colname` default uses column(s): `colData(SE)[[centerGroups]]
   * `sample_color_list` optional input for colorization

## 28jun2022

* `ggjammaplot()` issues:

  * `bw_factor` was behaving exactly opposite as expected
  * COMPLETE: color gradient by `ggplot2::stat_density_2d()` is inconsistent
  between outlier and non-outlier point ranges.

* `jammaplot()` and `ggjammaplot()` may need one argument to adjust
detail of density plots overall. Let it handle passing `binpi` and `bwpi`
to `jamba::plotSmoothScatter()`, or using ggplot2 geom arguments.

## 31may2022

* COMPLETE: Related to indicating `controlSamples` below, some plot hook to add
annotations to each panel while plotting for R base graphics. This
feature is likely already possible with `ggjammaplot()` with custom
ggplot2 layers.

   * COMPLETE: new argument `plot_hook_function` allows full customization
   after each panel is rendered, for `jammaplot()` in base R graphics.
   * Note: `ggjammaplot()` is fairly slow for this purpose, rendering
   may overall be slower than base R graphics, though this difference
   could be from lack of particular optimizations.

* Visual differences in `jammaplot()` and `ggjammaplot()` assignment of
colors to point density. Unclear whether this difference is due to
cropping of points outside the visual range, as happens with base R
graphics `plotSmoothScatter()` and `applyRangeCeiling=TRUE`.
* `ggjammaplot()` does not honor the order of samples when applying
`facet_wrap()`, it should convert the facet column to factor with
levels equal to the order of columns to be plotted.
* Consider `subtitle` being able to use one or more `colnames` in
`colData(se)` when the input data is `SummarizedExperiment`.

   * When multiple columns are provided, include them as multiple lines,
   each colored using `colorSub`?

* Consider using `shadowText()` to display the MAD values, otherwise
it is not visible in some panels.

## 26may2022

* COMPLETE: There should be some indication for samples that are `controlSamples`
during centering.

   * Perhaps an asterisk in the topleft corner inside each plot panel?

* Allow `sample_color_list` input as alternative to `colorSub`?
Or auto-detect `list` input and try to match values in the `list`.
* COMPLETE: Update the usage of `drawLabel()` to size the title box at least
the width of each plot panel.
* FIXED: `volcano_plot()` throwing an error:

   * ```"Error in utils::modifyList(default_params, new_values):
   is.list(x) is not TRUE```
   * `update_function_params(function_name = "volcano_plot", param_name = "color_set", 
    new_values = color_set) at jam-volcano-plot.R#374`

## 11may2022

* COMPLETE: `jammaplot()` should have some ability to provide column labels,
in place of using `colnames(x)` which may be a super-long text string.
* `jammaplot()` consider using `jamba::adjustAxisLabelMargins()` for
panel margins by default, making `margins` optional for custom use.
This change would ensure each margin fully displays text labels,
closer to how ggplot2 works.

Potential bugs:

* COMPLETE: `jammaplot()` highlightPch point shapes are not honored in the legend.
* COMPLETE: `jammaplot()` draws plot labels after `highlightPoints`, which can obscure
the `highlightPoints`. Ideally, draw the labels then the highlighted points.
This bug may not be evident with `ggjammaplot()`.

   * Partially completed by moving title box labels outside the plot. It
   is still possible with subtitle box labels, but will leave as-is for now.

* `ggjammaplot()` does not display `subtitle` box in the bottom-left.
Consider using `ggtext` or `ggplot2::geom_label()`.
* `jammaplot()` and `ggjammaplot()` should somehow indicate which samples
were used as controlSamples. Perhaps asterisk "*" in the title?
* The MAD factor label is difficult to read when it overlaps points,
particularly highlight points.
* When only one sample per centerGroups, hide MAD factor `"NaN"`.

## outlier detection in jammaplot by leave-one-out

Basic workflow:

* within each `centerGroups` grouping
* iterate each sample to leave one sample out
* call `jammacalc()` to calculate MADfactors for all replicates,
focusing on the omitted sample.
* determine if there are robust sample outliers

## volcano_plot()

* COMPLETE. Add new function.
* Add ggplot2 variation on this function.
* Add `volcano_sestats()` that takes `sestats` input from
`jamses::se_contrast_stats()`.

## ggjammaplot()

* make empty plot panels completely blank, relevant when using `blankPlotPos`
* set `xlab`, `ylab` using summary, difference labels
* add subtitleBox to bottomleft corner
* investigate using `element_text()` instead of `ggtext::element_textbox()`
* allow `geom_text_repel()` to label highlighted points
* consider selectable x- and y-ranges, to highlight a box and points within it

* Needed a workaround to build pkgdown site, see: https://github.com/r-lib/pkgdown/issues/1157

## New functions

* `jammaplotDispEsts()` wrapper for `DESeq2::plotDispEsts()`, although
this function needs the newer colorized `plotSmoothScatterG2()`.


## compatibility with SummarizedExperiment


## Bug fixes

* `jammaplot()` argument `ablineV` is not functioning properly.
* `jammaplot()` highlight legend is not using `highlightPch`.

## Vignettes

1. Guides for MA-plots.

* Basic guide to MA-plots for gene expression data.

  * check data normality, apply appropriate transform, normalize data
  * check within-group variability
  * highlight points
  * common patterns and what they mean:
  
    * shifted up/down
    * skewed up/down
    * low signal-to-noise, no signal
    * the 45-degree lines
    * batch effects
    * technical versus biological controls
    
  * non-parametric (rank-based) MA-plot

* When to do use `centerGroups`, `controlSamples`.
* How to interpret global-, group-, and technical replicate-centered data.
* How to detect sample outliers using a MAD factor threshold.
* How to guide data normalization by MA-plot review.


## useful functions

* `log2fold_to_fold()` and `fold_to_log2fold()` - to interconvert:

   * normal space fold changes, which could be represented as positive and
   negative values (e.g. 2-fold and -2-fold) or as ratios
   (e.g. 2-fold and 0.5-fold); and
   * log2 fold change values
   * Consider `fold_sign()` function that takes either as input and
   returns either `1` or `-1`, and by default never returns `0` since
   1-fold change cannot easily be multiplied by its `fold_sign()` without
   causing it to become zero. More thought to be had as to the workflow.


## Interface with other R packages

### DESeq2

* `plotMA()` is used for the `DESeq2::DESeqResults` object,
highlighting points that meet `alpha < 0.1` by default.
* Add new method to run `jammaplot()` on the count data,
optionally transformed using their recommended approach,
or use `useRank=TRUE`.


## Bug fixes

* Version 0.0.10.900 fixed a bug where rowGroupMeans() was used to
center values, but used default `na.rm=FALSE`, which caused groups
with missing values not to display a centered value for other non-NA
samples. The new version should provide two options:

    1. Hide or display values where groups with n>1 members have
    only one non-NA value. The reasoning is that the MA-plot is
    intended to show "difference from mean" and has far less utility
    when a fraction of rows has only one non-NA value in the group, 
    thus emphasizing the display of "zero difference" in those
    samples. Hiding these values therefore helps display the
    variability where data is available to judge the difference
    from mean.
    2. Hide groups with any NA values.



## Refactoring ideas

* Consider using `jamba::rowGroupMeans()` within `centerGeneData()`
which would allow optional outlier detection, which could substantially
improve the quality of MA-plot panels.
* `jammagg()` or something similarly named, to produce a ggplot2 object.
* One goal for ggplot2 output is to produce
interactive visualizations with plotly; however, if that is not
sufficient for the desired properties, e.g. useful hover text over
the smooth scatter density, then a custom interactive plot function
may be needed, for example direct calls to plotly instead of using
the ggplotly wrapper function.


### Update to handle SummarizedExperiment objects

* centerGroups should be able to take colnames from colData(x)
to define groupings.
* `jammaplot.SE()` could be specific to SummarizedExperiment objects.
It would require the names(assays(SE)) to define the data matrix to use.
It could even recognize sample groups from `colData(SE)`, or
from internal design matrix used for statistical testing.


### Highlighted points

* Consider allowing highlighted points per panel, for example,
showing gene hits only in the panels relevant to that comparison.
However, there is a workaround, using `whichSamples` which will
create the MA-plot data, but only display the samples of interest.


### Color legend for highlighted points

* Add optional color legend at bottom of the figure labeling the
colors highlighted on the plot.

### Enhance the returned data

* Currently returns a list of 2-column matrices
* Return data.frame with properly labeled colnames
* Add columns:
    * indicating the centerGroups value if defined
    * whether the sample is a control sample within its centerGroup
    * highlight label, highlight point color
    * title box color

User should be able to take the returned data.frame list, and make
a tall data.frame sufficient for ggplot2, or sufficient to answer
the question "What is that outlier datapoint?"


### Interactive plots using plotly

Future idea:

* Convert plotSmoothScatter to generate a heatmap-ready density
plot for each panel.
* For each cell in the panel, optionally update the label so it
lists the point label(s) represented in that cell; alternatively
simply label with the number of points represented.

### Accessory function to select points

* A common workflow request stems from the question "What are those
points?" A wrapper function was written to allow clicking twice inside
an active plot device to define a rectangle, then returning the
rownames of the corresponding points.
* The steps above could be used to supply coordinates from something
like R-shiny, to request the rownames of the points specified.
