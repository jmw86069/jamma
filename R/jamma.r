##
## jamma.R
##
## Targeted R functions as required to create MA-plots Jam-style
##
## Functions included:
## jamma
## centerGeneData
## points2polygonHull
##
## Requires: KernSmooth,sp,rgeos
## Suggests: matrixStats
##


#' jamma: MA-plots for omics data
#'
#' The jamma package creates MA-plots for omics data, with
#' important customizations to handle specific controls, and
#' experimental subsets.
#'
#' @section MA-plots functions:
#'    jammaplot,
#'    centerGeneData
#' @section Additional plot functions:
#'    points2polygonHull
#'
#' @docType package
#' @name jamma
NULL


#' Produce MA-plot of omics data.
#'
#' Produce MA-plot of omics data.
#'
#' `jammaplot` takes a numeric matrix, typically of gene expression data,
#' and produces an MA-plot (Bland-Altman plot), also known as a
#' median-difference plot. One panel is created for each column of
#' data. Within each panel, the x-axis represents the mean or median
#' expression of each row; the y-axis represents the difference from
#' mean or median for that column.
#'
#' By default, the plot uses `jamba::plotSmoothScatter()`, with optional
#' highlighted points draw using `points()`.
#'
#' The function will determine an appropriate layout of plot panels,
#' which can be overridden using `ncol` and `nrow` to specify the
#' number of columns and rows of plot panels, respectively. For now,
#' this function uses base R graphics instead of ggplot2, in order
#' to accomodate some custom features.
#'
#' This function uses "useRaster=TRUE" by default, which causes
#' `jamba::plotSmoothScatter()` to render a rasterized image as opposed
#' to a composite of colored rectangles. This process substantially
#' reduces the render time in all cases, and reduces the image size
#' when saving as PDF or SVG.
#'
#' ## Notable features
#'
#' ### Highlighting points
#'
#' Specific points can be highlighted with argument `highlightPoints`
#' which can be a vector or named list of vectors, containing `rownames(x)`.
#' When using a list, point colors are assigned to each element in the
#' list in order, using the argument `highlightColor`.
#'
#' ### Centering by control samples
#'
#' Typical MA-plots are "global-centered", which calculates the
#' mean/median across all columns in `x`, and this value is subtracted
#' from each individual value per row.
#'
#' By specifying `controlSamples`
#' the mean/median is calculated using only the `colnames(x)` which match
#' `controlSamples`, thus representing "difference from control."
#'
#' It may also be useful to center data by known high-quality samples,
#' so the effect of potential outlier samples is avoided.
#'
#' ### Centering within subgroups
#'
#' By specifying `centerGroups` as a vector of group names,
#' the centering is calculated within each group of `colnames(x)`.
#' In this way, subsets of samples can be treated independently in
#' the MA-plots. A good example might be producing MA-plots for
#' "kidney" samples, and "muscle" samples, which may have
#' fundamentally different signal distributions. A good rule
#' of thumb is to apply `centerGroups` to represent separate
#' groups of samples where you do not intend to apply direct
#' statistical comparisons across those samples, without at
#' least applying a two-way contrast, a fold change of fold
#' changes.
#'
#' Another informative technique is to center by sample group,
#' for example `centerGroups=sample_group`.
#' This technique produces MA-plots that depict the
#' "difference from group" for each sample replicate of a sample
#' group, and is very useful for identifying sample replicates
#' with markedly higher variability to its sample group than
#' others. In general, the variability within sample group
#' should be substantially lower than variability across
#' sample groups. Use `displayMAD=TRUE` and `outlierMAD=2`
#' as a recommended starting point for this technique.
#'
#' ### Applying a noise floor
#'
#' The argument `filterFloor` provides a numeric lower threshold,
#' where individual values **at or below** this threshold are
#' set to a defined value. The default defined value is the
#' floor itself, which has the effect of removing information from
#' points that are already below the noise threshold and therefore
#' are unreliable for this purpose.
#'
#' Another useful alternative
#' is to define `filterFloor=0` and `filterFloorReplacement=NA`
#' so that values of zero `0` are set to `NA` and are not included
#' in the MA-plot calculations, and are not represented as points
#' in each MA-plot panel. For data with many sparse missing values
#' represented as zero, these options can be very helpful because
#' each MA-plot panel will only represent actual measurements,
#' compared to only those sample which also have actual
#' measurements for those rows.
#'
#' ### Customizing the panel layout
#'
#' Panels are drawn using the order of `colnames(x)` by row,
#' from left-to-right, then top-to-bottom.
#' The argument `blankPlotPos` is intended to insert an empty panel
#' at a particular panel position, to help customize the alignment
#' of sample panels.
#' This option is typically used with `ncol` and `nrow` to define
#' a fixed layout of panel columns and rows. `blankPlotPos` refers
#' to panels numbered as drawn per row of panels,
#'
#' ### Identifying potential sample outliers
#'
#' Use argument `displayMAD=TRUE` to display the per-sample MAD factor
#' relative to its `centerGroups` value, if provided. The MAD value
#' for each MA-plot panel is calculated using rows whose mean
#' is at or above `outlierRowMin`. The median MAD value is calculated
#' for each `centerGroups` grouping when `groupedMAD=TRUE`, by default.
#' Finally, each MA-plot panel MAD factor is the ratio of its MAD value
#' to the relevant median MAD value. MA-plot panels with MAD factor
#' above `outlierMAD` are considered outliers, and the color ramp
#' uses `outlierColramp` or `outlierColor` as a visual cue.
#'
#' Putative outlier samples should usually not be determined
#' when:
#'
#'  * `controlSamples` are defined to include only a subset
#' of sample groups,
#' * `centerGroups` is not defined, or represents more than one
#' set of sample groups that are not intended to be statistically
#' compared directly to one another.
#'
#' Putative outlier samples may be defined when:
#'
#' * `centerGroups` represents a set of sample groups that are
#' intended to be involved in direct comparisons
#' * `centerGroups` represents each sample group
#'
#' Potential sample outliers may be identified by setting a threshold
#' with `outlierMAD`, by default 5xMAD. For a sample to be considered
#' an outlier, its median difference from mean/median needs to be
#' five times higher than the median across samples.
#'
#' We typically recommend an `outlierMAD=2` when centering
#' by sample groups, or when centering within experiment subsets.
#' For one sample to have 2xMAD factor, its variance needs
#' to be uniquely twice as high as the majority of other samples, which
#' is typically symptomatic of possible technical failure.
#'
#' There are exceptions to this suggested guideline, which includes
#' scenarios where a batch effect may be involved.
#'
#'
#' To do:
#' \itemize{
#'    \item{Accept other object types as input, including Bioconductor
#'       classes: ExpressionSet, SummarizedExperiment, MultiExperimentSet.}
#'    \item{Make it efficient to convey group information, for example
#'       define titleBoxColor with group colors, allow centerByGroup==TRUE
#'       which would re-use known sample group information.'}
#'    \item{Adjust the suffix to indicate when \code{centerGroups} are being
#'       used. For example indicate 'sampleID vs groupA' instead of
#'       'sampleID vs median'.}
#' }
#'
#' @param x `numeric` `matrix` typically containing log-normal
#'    values, with measurement rows, and sample columns. For example,
#'    with gene or protein expression data, the genes or proteins are
#'    represented in rows, and biological samples are represented
#'    in columns.
#' @param colramp one of several inputs recognized by
#'    `jamba::getColorRamp()`. It typically recognizes either the name of
#'    a color ramp from RColorBrewer, the name of functions from the
#'    `viridis` package such as `viridis::viridis()`, or single R colors, or
#'    a vector of R colors. When a single color is supplied, a gradient
#'    is created from white to that color, where the default base color
#'    can be customized with `defaultBaseColor="black"` for example.
#' @param colrampOutlier one of several inputs recognized by
#'    `jamba::getColorRamp()` to define a specific color ramp for
#'    MA-plot outlier panels, used when `outlierMAD` is defined.
#'    When `colrampOutlier` is `NULL` the `outlierColor` is used.
#' @param outlierColor `character` string representing one R color,
#'    used when `colrampOutlier` is `NULL` and when `outlierMAD` is
#'    defined. This color is used for MA-plot outlier panels by
#'    substituting the first color from the `colramp` color ramp,
#'    to act as a visual cue that the panel represents an outlier.
#' @param applyRangeCeiling `logical` passed to
#'    `jamba::plotSmoothScatter()` which determines how to handle points
#'    outside the plot x-axis and y-axis range: `applyRangeCeiling=TRUE`
#'    will place points at the border of the plot, which is helpful
#'    to indicate that there are more points outside the viewing range;
#'    `applyRangeCeiling=FALSE` will crop and remove points outside
#'    the viewing range, which is helpful for example when a large
#'    number of points are at zero and overwhelm the point density.
#'    When there are a large proportion of values at zero, it
#'    can be helpful to apply `xlim=c(0.01, 20)` and
#'    `applyRangeCeiling=FALSE`.
#' @param whichSamples `NULL` or `integer` vector, representing
#'    an index subset of samples to include in the MA-plots. When
#'    `whichSamples` represents a subset of samples in `x`, the
#'    MA-plot calculations are performed on all samples, then only
#'    samples in `whichSamples` are displayed. This argument keeps
#'    the MA-plot calculations consistent even when viewing only
#'    one or a subset of samples in more detail.
#' @param maintitle `character` string with the main title for
#'    the overall plot, printed at the top of each page in the
#'    top outer margin.
#' @param maintitleCex `numeric` cex character expansion used to
#'    resize the `maintitle`.
#' @param subtitle `NULL` or `character` vector to be drawn at
#'    the bottom left corner of each plot panel, the location
#'    is defined by `subtitlePreset`.
#' @param subtitlePreset character value describing where to position the
#'    subtitle, using terms valid in `jamba::coordPresets()`. The default
#'    `subtitlePreset="bottomleft"` places the subtitle at the
#'    bottom left corner of each plot panel.
#' @param titleBoxColor,subtitleBoxColor `character` vector of
#'    R colors applied to title text, or subtitle text, respectively.
#'    When `doTitleBox=TRUE`
#'    one or no value is supplied, it defines colors using
#'    `jamba::setTextContrastColor()` to use a contrasting color.
#' @param titleColor `character` vector of colors applied to title text
#'    in each MA-plot panel. When `doTitleBox=TRUE` and `titleColor`
#'    contains only one or no value, the title color is defined by
#'    `jamba::setTextContrastColor()` along with `titleBoxColor`.
#' @param doTitleBox `logical` indicating whether to draw plot
#'    titles using a colored box. When `doTitleBox=TRUE` the
#'    `jamba::drawLabels()` is called to display a label box at
#'    the top of each plot panel, with `drawBox=TRUE`. When
#'    `doTitleBox=FALSE`, `jamba::drawLabels()` is called with
#'    `drawBox=FALSE`.
#' @param titleFont `integer` font compatible with `par("font")`.
#'    Values are recycled across panels, so each panel can use a custom
#'    value if needed.
#' @param titlePreset `character` value describing where to position the
#'    subtitle, using terms valid in `jamba::coordPresets()`.
#'    The default `titlePreset="top"` centers the label at the
#'    top of each panel.
#' @param xlab,ylab `character` x- and y-axis labels, respectively.
#'    The default values are blank `""` because there are a wide variety
#'    of possible labels, and the labels take up more space
#'    than is often useful for most MA-plots.
#' @param xlabline,ylabline `numeric` number indicating the text line
#'    distance from the edge of plot border to place `xlab` and `ylab`
#'    text, as used by `graphics::title()`.
#' @param groupSuffix `character` text appended to each MA-plot
#'    panel title. This argument is deprecated in favor of
#'    using `subtitle` to place additional text at the bottom left
#'    corner of each MA-plot panel.
#' @param highlightPoints `NULL`, or `character` vector, or a
#'    `list` of `character` vectors indicating `rownames(x)` to
#'    highlight in each MA-plot panel. When `NULL`, no points are
#'    highlighted; when `character` vector, points are highlighted in
#'    all MA-plot panels; when `list` of `character` vectors, each
#'    `character` vector in the list is highlighted using a unique
#'    color in `highlightColor`. Points are drawn using
#'    `graphics::points()` and colored using `highlightColor`,
#'    which can be time-consuming for a large number of highlight
#'    points.
#' @param highlightColor `character` vector used when `highlightPoints`
#'    is defined. It is recycled to `length(highlightPoints)` and
#'    is applied either to
#' @param highlightCex `numeric` value recycled to `length(highlightPoints)`
#'    indicating the highlight point size.
#' @param doHighlightPolygon `logical` indicating whether to draw
#'    a shaded polygon encompassing `highlightPoints`, using
#'    `highlightColor`.The polygon is defined by `grDevices::chull()`
#'    via the function `points2polygonHull()`.
#' @param highlightPolygonAlpha `numeric` value indicating alpha
#'    transparency, where `0` is fully transparent, and `1` is completely
#'    not transparent.
#' @param doHighlightLegend `logical` indicating whether to print a
#'    color legend when `highlightPoints` is defined. The legend is
#'    displayed in the bottom outer margin of the page using
#'    `outer_legend()`, and the page is adjusted to add bottom
#'    outer margin.
#' @param smoothPtCol `color` used to draw points when `nrpoints` is
#'    non-zero, which draws points in the extremities of the
#'    smooth scatter plot. See `jamba::plotSmoothScatter()`.
#'    The effect can also be achieved by adjusting `transFactor` to
#'    a lower value, which increases the visual contrast of individual
#'    points in the point density.
#' @param margins `numeric` vector of margins compatible with
#'    `graphics::par("mar")`. Default values
#'    are applied, but provided here for convenience.
#' @param useRaster `logical` indicating whether to draw the
#'    smooth scatter plot using raster logic, `useRaster=TRUE` is
#'    passed to `jamba::plotSmoothScatter()`. The default `TRUE`
#'    creates a much smaller plot object by rendering each plot
#'    panel as a single raster image instead of rendering individual
#'    colored rectangles.
#' @param ncol,nrow `integer` number of MA-plot panel columns and rows
#'    passed to `graphics::par("mfrow")` when `doPar=TRUE`. When only one
#'    value is supplied, `nrow` or `ncol`, the other value is defined
#'    by `ncol(x)` and `blankPlotPos` so all panels can be contained on
#'    one page. When `nrow` and `ncol` are defined such that multiple
#'    pages are produced, each page will be annotated with `maintitle`
#'    and `doHighlightLegend` as relevant.
#' @param doPar `logical` indicating whether to apply
#'    `graphics::par("mfrow")` to define MA-plot panel rows and columns.
#'    When `doPar=FALSE` each plot panel is
#'    rendered without adjusting the `graphics::par("mfrow")` setting.
#' @param las `integer` value `1` or `2` indicating whether axis labels
#'    should be parallel or perpendicular to the axes, respectively.
#' @param ylim,xlim `NULL` or `numeric` vector `length=2` indicating
#'    the y-axis and x-axis ranges, respectively. The values are useful
#'    to define consistent dimensions across all panels. The
#'    default `ylim=c(-4,4)` represents 16-fold up and down range in
#'    normal space, and is typically a reasonable starting point for
#'    most purposes. Even if numeric values are all between
#'    `-1.5` and `1.5`, it is still recommended to keep a range in
#'    context of `c(-4, 4)`, to indicate that the observed values
#'    are lower than typically observed. The `c(-4, 4)` may be adjusted
#'    relative to the typical ranges expected for the data.
#'    It is sometimes helpful to define `xlim` slightly above zero for
#'    datasets that have an extremely large proportion of zeros, in order
#'    to reduce the visual effect of having that much point density at
#'    zero, for example with `xlim=c(0.001, 20)` and
#'    `applyRangeCeiling=FALSE`.
#' @param controlSamples `character` vector of `colnames(x)` to define
#'    control samples during the data centering step. Default MA-plots are
#'    calculated relative to the overall mean, which uses all `colnames(x)`
#'    as `controlSamples`. However it can be helpful and informative
#'    to view MA-plots relative to a set of known control samples.
#' @param centerGroups vector of groups, with length `ncol(x)`
#'    indicating optional groupings to use during data centering.
#'    When `centerGroups` is defined, each subset of data is centered
#'    independently. It is useful to center within batches, or within
#'    subsets of samples that are not intended to be compared across
#'    one another. It is also informative to center by each sample
#'    group in order to view the variability among sample group
#'    replicates, which should be much lower than variability across
#'    sample groups. See `centerGeneData()` for more specific examples.
#' @param useMean `logical` indicating whether to center data using
#'    `mean` or `median` values. The default `useMean=TRUE` is chosen
#'    because it visually represents data as seen by typical parametric
#'    analysis methods downstream. When outlier sample are observed,
#'    it may be more useful to apply `useMean=FALSE` to use the `median`
#'    which is less prone to outlier effects. That said, if a particular
#'    sample is an outlier, another alternative is to define
#'    `controlSamples` which excludes the outlier sample(s), and
#'    therefore data centering is applied using only the non-outlier
#'    samples as the reference.
#' @param customFunc `NULL` or `function` used instead of `mean` or
#'    `median` during the data centering step to generate a row summary
#'    statistic. It should take `matrix` input, and return a `numeric` vector
#'    output summarizing each row in `x`, to be subtracted from each
#'    `numeric` value by row in `x`. It is intended to
#'    provide custom row statistics, for example geometric mean, or other
#'    row summary function.
#' @param filterNeg `logical` deprecated argument, use `filterFloor`
#'    instead. The `filterNeg` indicates whether to change all negative
#'    values to zero `0` before proceeding with data centering.
#'    Negative values are often the result of measurements
#'    being below a noise threshold in upstream data processing,
#'    and therefore the magnitude of negative value is usually
#'    either not informative, or not on similar scale as positive
#'    values.
#'    When `filterNeg=TRUE`, negative values are set to zero, and
#'    can result in a characteristic 45 degree angle line originating at
#'    `x=0` extending to the right.
#' @param filterNA,filterNAreplacement `logical` and `vector` respectively.
#'    When `filterNA=TRUE`, all `NA` values are replaced with
#'    `filterNAreplacement`, which can be helpful to handle `NA` values
#'    as zero `0` for example. In reality, `NA` values should probably
#'    be left as-is, so subsequent data centering does not use these values,
#'    and so the MA-plot panel does not draw a point when no measurement
#'    exists.
#' @param filterFloor,filterFloorReplacement `numeric` or `NULL` indicating
#'    a numeric floor, where values in `x` **at or below** `filterFloor` are
#'    replaced with `filterFloorReplacement`. Note that this argument
#'    can be used to replace zero `0` with `NA` in the event that
#'    zeros do not represent measurements. One can typically tell whether
#'    input data includes zero `0` values by the presence of characteristic
#'    45-degree angle lines originating from `x=0` angled to the right.
#'    The default values replace any values **at or below zero** with zero,
#'    which also applies a numeric floor to negative values.
#'    For some platform data technologies, it can be useful to define a
#'    `filterFloor` roughly equivalent to its noise threshold. For example
#'    quantitative PCR sometimes uses `log_expression = (40 - Ct)`, where
#'    `Ct` values above `35` are considered to be noise. That noise threshold
#'    implies that any `expression` values `5` or lower are roughly
#'    equivalent noise, so applying `filterFloor=5` is appropriate.
#' @param transFactor `numeric` adjustment to the visual density of
#'    points by `jamba::plotSmoothScatter()`. The value is based upon
#'    `graphics::smoothScatter()` argument `transformation` which uses
#'    `function(x)x^0.24`. The `transFactor` is equivalent to the
#'    exponential in the form: `function(x)x^transFactor`. Lower values
#'    increase the visual density more intense, higher values make the
#'    visual density less intense.
#' @param nrpoints `integer` or `NULL` indicating the number of points
#'    to display on the extremity of the smooth scatter density,
#'    passed to `jamba::plotSmoothScatter()`.
#' @param smoothScatterFunc `function` used for the smooth scatter plot,
#'    default `jamba::plotSmoothScatter()`. Note that a custom function
#'    may not recognize `nrpoints` or `transformation`.
#' @param ablineH,ablineV `numeric` vector indicating horizontal and
#'    vertical lines to draw in each MA-plot panel, respectively.
#'    These values are passed to `graphics::abline()`.
#' @param doTxtplot `logical` not yet implemented, indicating whether
#'    to plot results using colored text output.
#' @param blankPlotPos `NULL` or `integer` vector indicating
#'    plot panel positions to be drawn blank, skipped. Blank panel
#'    positions are intended to help customize the visual alignment
#'    of MA-plot panels. The mechanism is similar to `ggplot2::facet_wrap()`
#'    except that blank positions can be manually defined by what makes
#'    sense to the experiment design.
#' @param displayMAD `logical` indicating whether to display each
#'    MA-plot panel MAD factor (median absolute deviation). A MAD
#'    value for each panel is calculated by taking the median absolute
#'    deviation from zero across all points, using points whose mean
#'    value is equal or greater than `outlierRowMin`. The overall MAD
#'    is defined by the median MAD from the MA-plot panels. The MAD
#'    factor is defined as the ratio of each MA-plot panel MAD value
#'    to the overall MAD value, and therefore most MAD factor values
#'    should be roughly `1`. The overall MAD value is defined by the
#'    median across all samples when `groupedMAD=FALSE`, or defined
#'    within each `centerGroup` when `groupedMAD=TRUE`. A value
#'    with MAD factor 2 is interpreted as a sample whose median
#'    deviation from zero is twice as high as the typical sample,
#'    which is a reasonably indication that this sample has twice
#'    the inherent level of noise compared to other samples. Note
#'    that MAD values should be interpreted within sample processing
#'    batches if relevant, or within logical experimental units --
#'    roughly interpreted to mean sets of samples within which direct
#'    statistical comparisons are intended to be applied. For example,
#'    gene expression data that include brain and liver samples would
#'    probably use `centerGroups` for brain and liver to be centered
#'    separately, therefore the MAD factors should be separately
#'    calculated for brain and for liver.
#' @param groupedMAD logical indicating how the MAD calculation
#'    should be performed: `groupedMAD=TRUE` will calculate the
#'    median MAD and corresponsing MAD factor within each
#'    `centerGroups` grouping; `groupedMAD=FALSE` will calculate
#'    one overall median MAD, and corresponding MAD factor values
#'    will be performed across all samples.
#' @param outlierMAD `numeric` threshold above which a MA-plot panel
#'    MAD factor is considered an outlier. When a MA-plot panel is
#'    considered an outlier, the `outlierColramp` or `outlierColor`
#'    is applied to the panel color ramp to display a visual
#'    indication.
#' @param outlierRowMin `numeric` value indicating the minimum mean
#'    value as displayed on the MA-plot panel x-axis, in order for
#'    the row to be included in MAD calculations. This argument is
#'    intended to prevent measurements whose mean value is below
#'    a noise threshold from being included, therefore only including
#'    points whose mean measurement is above noise and represents
#'    "typical" variability.
#' @param fillBackground `logical` passed to `jamba::plotSmoothScatter()`
#'    indicating whether to fill the plot panel with the background color,
#'    using the first value in the color ramp for each MA-plot panel.
#' @param ma_method character string indicating whether to perform
#'    MA-plot calculations using the old method `"old"`; or `"jammacalc"`
#'    which uses the function `jammacalc()`.
#' @param doPlot `logical` indicating whether to create plots. When
#'    `doPlot=FALSE` only the MA-plot panel data is returned.
#' @param useRank `logical` indicating whether to create column-wide
#'    ranks, then create MA-plots using the rank data. When `useRank=TRUE`
#'    the y-axis represents the rank difference from mean, and the
#'    x-axis represents the mean rank. Using `useRank=TRUE` is a good
#'    method to evaluate whether data can be normalized, or whether
#'    data across samples is inherently noisy.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional parameters sent to downstream functions,
#'    \code{\link{jamba::plotSmoothScatter}}, \code{\link{centerGeneData}}.
#'
#' @return List of data used by jamma() to produce an MA-plot, sufficient
#'          to use as input to produce another MA-plot.
#'
#' @seealso `jamba::plotSmoothScatter()` for enhanced
#'    `graphics::smoothScatter()`, `centerGeneData()`,
#'    `jamba::imageDefault()`, `jamba::nullPlot()`
#'
#' @aliases MAplot, MVAplotMED, MVAplotMEDsmooth, MVAplot, plotMA, plotMVA,
#'    ma.plot, plot.ma, mva.plot
#'
#' @family jam plot functions
#'
#' @examples
#' # Note the example data requires the affydata Bioconductor package
#' if (suppressPackageStartupMessages(require(affydata))) {
#'    data(Dilution);
#'    edata <- log2(1+exprs(Dilution));
#'    jammaplot(edata);
#' }
#'
#' @export
jammaplot <- function
(x,
 maintitle=NULL,
 titleBoxColor="#DDBB9977",
 subtitleBoxColor=titleBoxColor,
 centerGroups=NULL,
 controlSamples=colnames(x),
 useMean=TRUE,
 ylim=c(-4,4),
 xlim=NULL,
 highlightPoints=NULL,
 outlierMAD=5,
 outlierRowMin=5,
 displayMAD=FALSE,
 groupedMAD=TRUE,
 colramp=c("white", "lightblue", "blue", "navy", "orange", "orangered2"),
 colrampOutlier=NULL,
 outlierColor="palegoldenrod",
 whichSamples=NULL,
 maintitleCex=1.8,
 subtitle=NULL,
 subtitlePreset="bottomleft",
 titleCexFactor=1,
 titleCex=NULL,
 doTitleBox=TRUE,
 titleColor="black",
 titleFont=2,
 titlePreset="top",
 xlab="",
 xlabline=2,
 ylab="",
 ylabline=1.5,
 groupSuffix=NULL,
 highlightPch=21,
 highlightCex=1.5,
 highlightColor="#00AAAA66",
 doHighlightPolygon=FALSE,
 highlightPolygonAlpha=0.3,
 doHighlightLegend=TRUE,
 smoothPtCol="#00000055",
 margins=c(3.5, 2, 0.3, 0.2),
 useRaster=TRUE,
 ncol=NULL,
 nrow=NULL,
 doPar=TRUE,
 las=2,
 groupedX=TRUE,
 customFunc=NULL,
 filterNA=TRUE,
 filterNAreplacement=NA,
 filterNeg=TRUE,
 filterFloor=0,
 filterFloorReplacement=filterFloor,
 transFactor=0.18,
 nrpoints=0,
 smoothScatterFunc=jamba::plotSmoothScatter,
 applyRangeCeiling=TRUE,
 doTxtplot=FALSE,
 ablineV=0,
 ablineH=c(-2,0,2),
 blankPlotPos=NULL,
 fillBackground=TRUE,
 useRank=FALSE,
 ma_method=c("jammacalc", "old"),
 doPlot=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to provide an MA plot, based upon the MVAplotMED from
   ## the defunct Agi4x44PreProcess package,  except in this case using
   ## smoothScatter since 20k to 40k datapoints is too many to show without
   ## them being occluded.
   ##
   ## You must manually specify the number columns and rows, and set doPar=TRUE
   ## for this function to define multi-panel plotting, e.g.
   ## "nrow=5, ncol=4, doPar=TRUE"
   ## Otherwise, you can manually set these parameters with
   ## "par('mfrow'=c(5,4));"
   ##
   ## This function also uses "useRaster=TRUE" by default which causes the
   ## imageDefault() function used by smoothScatter() to produce a temporary
   ## raster image which is then resized with interpolation to produce a
   ## properly-blended representation of the full image. This process also
   ## substantially reduces the image size when saving as a vector output
   ## format like PDF or SVG.
   ##
   ## To highlight certain points, e.g. marked outliers, or statistical hits,
   ## supply rownames into highlightPoints as a vector, or a list of vectors.
   ## If using a list of vectors, also supply a list of colors in
   ## highlightColor, so each subset is colored uniquely.
   ##
   ## ctrlSamples is deprecated in favor of controlSamples.
   ##
   ## controlSamples indicates which samples should be used for centering
   ## data, by default all samples. However, centering by control samples
   ## visually indicates changes from those control samples.
   ##
   ## centerGroups is used to apply the centering to subsets of columns
   ## of data, for example if two very different cell types should be separately
   ## centered for closer scrutiny of the within-cell variability.
   ## One can take this further, and center each sample group, which has
   ## the helpful effect of displaying variability of replicates per
   ## sample group. This method is quite useful for reviewing whether any
   ## replicates are notably different than counterparts in the same
   ## sample group.
   ##
   ## Note: the "filterNAandNeg" will floor log2 values at 0, which has the
   ## effect of adding the "45 degree lines" originating from c(0,0).
   ##
   ## blankPlotPos is intended for specific plot layouts, to keep the plots
   ## organized by some fixed criteria, e.g. certain sample groups per
   ## row or column, but where there may not be consistent replicates per
   ## group.
   ## For ncol and nrow, blankPlotPos indicates which panels should be
   ## left empty, where each panel is numbered by row.  For example,
   ## if ncol=4, then the panels are numbered 1-4, 5-8, 9-12, etc.
   ##
   ## outlierMAD is only viable when centerEachGroup=TRUE; in which case
   ## outlierMAD is a threshold for the MAD factor of outliers.  Each MA
   ## plot panel is given a MAD score, using probes whose x-values are
   ## greater than outlierRowMin;, the median of the absolute value is
   ## used as the MAD.
   ## These MAD values are used to calculate the global MAD, and individual
   ## MADs are divided by global MAD to determine a MAD factor.
   ## The outlierMAD is the threshold above which a MAD factor is considered
   ## an outlier.
   ##
   #if (!suppressPackageStartupMessages(require(jamba))) {
   #   stop("jammaplot() requires the jamba package.");
   #}
   #if (suppressPackageStartupMessages(require(matrixStats))) {
   #   doMS <- TRUE;
   if ("matrixStats" %in% rownames(installed.packages())) {
      doMS <- TRUE;
      rowMedians <- matrixStats::rowMedians;
      colMedians <- matrixStats::colMedians;
   } else {
      doMS <- FALSE;
      rowMedians <- function(x, na.rm=TRUE) {
         apply(x, 1, median, na.rm=na.rm);
      }
      colMedians <- function(x, na.rm=TRUE) {
         apply(x, 2, median, na.rm=na.rm);
      }
   }
   ma_method <- match.arg(ma_method);
   if (doTxtplot) {
      blankPlotPos <- NULL;
      smoothScatterFunc <- function(...){
         plotTextSmoothScatter(height=20,
            width=80,
            doLegend=FALSE,
            ...);
      }
   }

   transformation <- function(x){
      x^transFactor;
   }
   if ("list" %in% class(x)) {
      nsamples <- length(x);
      x_names <- names(x);
   } else {
      nsamples <- ncol(x);
      x_names <- colnames(x);
   }
   if (length(x_names) == 0) {
      x_names <- jamba::makeNames(rep("V", ncol(x)));
      if ("list" %in% class(x)) {
         names(x) <- x_names;
      } else {
         colnames(x) <- x_names;
      }
   }

   ## if colrampOutlier is a single color, use it to replace
   ## the first color of colramp
   if (length(outlierColor) == 0) {
      outlierColor <- "palegoldenrod";
   }
   if (length(colrampOutlier) == 0) {
      colrampOutlier <- outlierColor;
   }
   if (!is.list(colramp)) {
      colramp <- list(colramp);
   }
   if (!is.list(colrampOutlier)) {
      colrampOutlier <- list(colrampOutlier);
   }
   colramp <- rep(colramp,
      length.out=length(x_names));
   names(colramp) <- x_names;
   colramp <- lapply(colramp, jamba::getColorRamp, n=NULL);

   colrampOutlier <- rep(colrampOutlier,
      length.out=length(x_names));
   names(colrampOutlier) <- x_names;
   # apply outlier colors to each color ramp
   for (i in seq_along(x_names)) {
      icolramp <- colramp[[i]];
      icolrampOutlier <- colrampOutlier[[i]];
      if (!is.function(icolrampOutlier)) {
         if (length(icolrampOutlier) == 1 && jamba::isColor(icolrampOutlier)) {
            ## substitute single color with the first color in colramp
            icolrampOutlier <- jamba::getColorRamp(
               c(icolrampOutlier,
                  tail(icolramp(101), -1)),
               n=NULL,
               ...);
         } else {
            ## For multi-color colrampOutlier, use it to define
            ## the complete outlier color gradient
            icolrampOutlier <- jamba::getColorRamp(icolrampOutlier,
               n=NULL,
               ...);
         }
         colrampOutlier[[i]] <- icolrampOutlier;
      }
   }

   ## If groupSuffix is supplied, make sure its length
   ## is length(x_names)
   if (length(groupSuffix) > 0) {
      groupSuffix <- rep(groupSuffix,
         length.out=length(x_names));
   }

   ## highlightPoints
   ## Todo: check values in highlightPoints versus rownames(x)
   if (length(highlightPoints) == 0) {
      highlightColor <- NULL;
   } else {
      if (!is.list(highlightPoints)) {
         if (length(highlightPoints) == length(highlightColor)) {
            if (length(names(highlightPoints)) == 0) {
               names(highlightPoints) <- makeNames(highlightPoints);
            }
            highlightPoints <- as.list(highlightPoints);
         } else {
            highlightPoints <- list(highlighted=highlightPoints);
         }
      } else {
         if (length(names(highlightPoints)) == 0) {
            names(highlightPoints) <- makeNames(
               rep("highlight",
                  length.out=length(highlightPoints)));
         }
      }
      if (!is.list(highlightColor)) {
         highlightColor <- as.list(highlightColor);
      }
      highlightColor <- rep(highlightColor,
         length.out=length(highlightPoints));
      if (length(names(highlightColor)) > 0 &&
            all(names(highlightColor) %in% names(highlightPoints))) {
         highlightColor <- highlightColor[names(highlightPoints)];
      } else {
         names(highlightColor) <- names(highlightPoints);
      }

      if (!is.list(highlightPch)) {
         highlightPch <- as.list(highlightPch);
      }
      highlightPch <- rep(highlightPch,
         length.out=length(highlightPoints));
      if (length(names(highlightPch)) > 0 &&
            all(names(highlightPch) %in% names(highlightPoints))) {
         highlightPch <- highlightPch[names(highlightPoints)];
      } else {
         names(highlightPch) <- names(highlightPoints);
      }
      if (!class(highlightCex) %in% c("list")) {
         highlightCex <- as.list(highlightCex);
      }
      highlightCex <- rep(highlightCex,
         length.out=length(highlightPoints));
      if (length(names(highlightCex)) > 0 &&
            all(names(highlightCex) %in% names(highlightPoints))) {
         highlightCex <- highlightCex[names(highlightPoints)];
      } else {
         names(highlightCex) <- names(highlightPoints);
      }
   }
   ## By default, when subtitle is not supplied, and
   ## centerGroups contains multiple unique values,
   ## use centerGroups in subtitle
   if (length(subtitle) == 0 &&
         length(unique(centerGroups)) > 1) {
      subtitle <- rep(centerGroups,
         length.out=length(x_names));
      names(subtitle) <- x_names;
   }

   ## Define a new par for panel layout
   if (doPar) {
      oPar <- par(no.readonly=TRUE);
      on.exit(par(oPar));
      if (is.null(ncol)) {
         if (is.null(nrow)) {
            parMfrow <- jamba::decideMfrow(nsamples);
            nrow <- parMfrow[1];
            ncol <- parMfrow[2];
         } else {
            ncol <- ceiling(nsamples / nrow);
         }
      } else if (is.null(nrow)) {
         nrow <- ceiling(nsamples / ncol);
      }
      par(mfrow=c(nrow, ncol));
      if (length(maintitle) > 0) {
         maintitle_nlines <- length(unlist(strsplit(maintitle, "\n")));
         par("oma"=pmax(par("oma"),
            c(0, 0, 1.5+1.5*maintitle_nlines, 0)));
      }
      if (length(highlightPoints) > 0 && doHighlightLegend) {
         par("oma"=pmax(par("oma"),
            c(3, 0, 0, 0)));
      }
   }
   if (length(titleCex) == 0) {
      titleCex <- 1;
   }
   titleCex <- rep(titleCex, length.out=nsamples);
   titleFont <- rep(titleFont, length.out=nsamples);
   titleBoxColor <- rep(titleBoxColor, length.out=nsamples);

   ## Adjust titleColor to have contrast from the titleBoxColor
   if (length(titleColor) <= 1) {
      if (doTitleBox) {
         titleColor <- jamba::setTextContrastColor(titleBoxColor);
      } else {
         if (length(titleColor) == 0) {
            titleColor <- "black";
         }
         titleColor <- rep(titleColor, length.out=nsamples);
      }
   }

   par("mar"=margins);
   nARR <- nsamples;
   if (is.null(whichSamples)) {
      whichSamples <- 1:nARR;
   } else {
      if (is.numeric(whichSamples)) {
         whichSamples <- x_names[whichSamples]
      }
      whichSamples <- match(whichSamples, x_names);
   }
   names(whichSamples) <- x_names[whichSamples];
   gaveMVA <- FALSE;
   if ("data.frame" %in% class(x)) {
      x <- as.matrix(x);
   }
   if (class(x) %in% "list" &&
      all(c("x", "y") %in% colnames(x[[1]]))) {
      ## Re-use MA-plot data as provided
      gaveMVA <- TRUE;
      if (is.null(xlim)) {
         xlim <- range(sapply(x, function(iObj){
            range(c(range(iObj[,"x"])));
         }));
         #xlim <- c(0, xlim1);
      }
      mvaDatas <- x;
   } else {
      ## Prepare MA-plot data
      if (ma_method %in% "jammacalc") {
         ## Newer method using jammacalc()
      } else {
         if (verbose) {
            jamba::printDebug("jammaplot(): ",
               "Processing matrix input type, dim(x)", dim(x));
         }
         if (filterNA) {
            x[is.na(x)] <- filterNAreplacement;
         }
         if (filterNeg) {
            x[x <= 0] <- 0;
         }
         if (!is.null(filterFloor)) {
            x[!is.na(x) & x <= filterFloor] <- filterFloorReplacement;
         }

         # apply rank here
         if (useRank) {
            if (is.list(x)) {
               stop("Cannot combine useRank=TRUE with list input, it requires numeric matrix input.");
            }
            ## Convert x to rank per column
            if (verbose) {
               jamba::printDebug("jammaplot(): ",
                  "Converting input matrix to per-column rank values.");
            }
            x1 <- apply(x, 2, rank, na.last=FALSE);
            if (any(is.na(x))) {
               x1[is.na(x)] <- NA;
            }
            x <- x1;
            if (all(abs(ylim) == 4)) {
               ylim <- c(-1,1) * round(nrow(x) * .40);
            }
         }

         if (length(controlSamples) == 0) {
            controlSamples <- x_names;
         } else {
            controlSamples <- intersect(controlSamples, x_names);
            if (length(controlSamples) == 0) {
               controlSamples <- x_names;
            }
         }
         ## Calculate the summary value per row of input data x,
         ## used as the x-axis value for each panel
         if (groupedX && length(unique(centerGroups)) > 1) {
            yM <- centerGeneData(x,
               controlSamples=controlSamples,
               centerGroups=centerGroups,
               returnValues=FALSE,
               returnGroupedValues=TRUE);
            centerGroups <- jamba::nameVector(centerGroups, x_names);
            ## TODO: Apply rowGroupsMeans() for custom y-values per group
            if (useMean) {
               y <- rowMeans(yM, na.rm=TRUE);
            } else {
               y <- rowMedians(yM, na.rm=TRUE);
            }
         } else {
            if (length(customFunc) > 0) {
               y <- customFunc(x[,controlSamples,drop=FALSE], na.rm=TRUE);
            } else if (useMean) {
               y <- rowMeans(x[,controlSamples,drop=FALSE], na.rm=TRUE);
            } else {
               y <- rowMedians(x[,controlSamples,drop=FALSE], na.rm=TRUE);
            }
         }
         if (verbose) {
            jamba::printDebug("jammaplot(): ",
               "Processing matrix input type, dim(x)", dim(x));
            jamba::printDebug("jammaplot(): ",
               "Processing matrix input type, dim(y):", dim(y));
         }
         if (is.null(xlim)) {
            ## Calculate x range using the range of data overall
            xlim <- range(c(range(x, na.rm=TRUE), range(y, na.rm=TRUE)));
            #xlimMax <- mean(c(max(object, na.rm=TRUE), max(y, na.rm=TRUE)));
            #xlim <- c(filterFloor, xlimMax);
         }
         mvaDatas <- lapply(jamba::nameVector(whichSamples), function(i){NA});
      }
   }

   if (ma_method %in% "jammacalc") {
      ## Newer method using jammacalc()
      if (verbose) {
         jamba::printDebug("jammaplot(): ",
            "Calling jammacalc().",
            "useMedian:", !useMean);
      }
      mvaDatas <- jammacalc(x=x,
         na.rm=TRUE,
         useMedian=!useMean,
         controlSamples=controlSamples,
         centerGroups=centerGroups,
         groupedX=groupedX,
         grouped_mad=groupedMAD,
         whichSamples=whichSamples,
         noise_floor=filterFloor,
         noise_floor_value=filterFloorReplacement,
         naValue=filterNAreplacement,
         centerFunc=centerGeneData_new,
         returnType="ma_list",
         mad_row_min=outlierRowMin,
         useRank=useRank,
         verbose=verbose,
         ...);
      # for useRank=TRUE adjust default ylim
      if (useRank && all(abs(ylim) == 4)) {
         ylim <- c(-1,1) * round(nrow(x) * .40);
      }
   } else {
      if (length(centerGroups) > 0) {
         objectCtr <- centerGeneData(x,
            centerGroups=centerGroups,
            mean=useMean,
            needsLog=FALSE,
            returnGroupedValues=FALSE,
            returnValues=TRUE,
            controlSamples=controlSamples,
            verbose=verbose);
         if (verbose) {
            jamba::printDebug("jammaplot(): ",
               "objectCtr:");
            print(head(objectCtr));
            jamba::printDebug("head(y): ",
               round(head(y), digits=2));
         }
      }

      ## If not given MVA data, then calculate the MVA data upfront
      ## so we can test for outliers beforehand
      if (!gaveMVA) {
         if (verbose) {
            jamba::printDebug("jammaplot(): ",
               "Calculating MVA data.");
         }
         for (i in whichSamples) {
            if (groupedX && length(unique(centerGroups)) > 1) {
               ## Optionally, if centerGroups are supplied,
               ## we center within each controlGroup,
               ## and plot the difference within each controlGroup
               iCol <- x_names[i];
               iGroup <- centerGroups[iCol];
               yUse <- yM[,iGroup];
               #x <- objectCtr[,i] + y;
               xi <- objectCtr[,i] + yUse;
               if (verbose) {
                  jamba::printDebug("jammaplot(): ",
                     "Applying unique x per centerGroups");
               }
            } else {
               ## Typical values use difference from median on y-axis and
               ## average of current with median on x-axis
               ##
               ## New Jam values will retain the y-axis,
               ## but use the median on the x-axis
               #mvaValues <- c((x + y)/2, (x - y));
               xi <- x[,i];
               yUse <- y;
            }
            mvaValues <- c(yUse, (xi - yUse));
            mvaData <- matrix(data=mvaValues,
               nrow=length(xi),
               byrow=FALSE,
               dimnames=list(genes=rownames(x),
                  MAcoords=c("x", "y")));
            mvaDatas[[i]] <- mvaData;
            names(mvaDatas)[i] <- x_names[i];
         }
      }
   }

   ## Optionally calculate the MAD per panel
   ## assumptions: that zero should be the median, all deviations are absolute values since
   ## we are comparing with zero, and we are not subtracting an actual median value
   if (ma_method %in% "jammacalc") {
      mvaMADs <- attr(mvaDatas, "MADs");
      mvaMADfactors <- attr(mvaDatas, "MADfactors");
   } else {
      if (verbose) {
         jamba::printDebug("jammaplot():",
            "Calculating MA data MAD factors.");
      }
      mvaMADs <- sapply(jamba::nameVector(whichSamples, x_names[whichSamples]), function(i){
         if (verbose) {
            jamba::printDebug("   i:",
               i);
         }
         mvaData <- mvaDatas[[x_names[i]]];
         iWhich <- (!is.na(mvaData[,"x"]) & abs(mvaData[,"x"]) >= outlierRowMin);
         jamba::rmNA(
            median(abs(mvaData[iWhich,"y"]),
               na.rm=TRUE),
            naValue=Inf);
      });

      ## Calculate MAD factor, a ratio to the median
      if (groupedMAD && length(centerGroups) > 0) {
         ## if grouped, use the median of each center group
         mvaMADfactors <- unlist(unname(
            tapply(mvaMADs, centerGroups, function(i){
               i / median(i, na.rm=TRUE);
            })))[names(mvaMADs)];
      } else {
         ## if not grouped, use the overall median across all samples
         mvaMAD <- median(mvaMADs, na.rm=TRUE);
         mvaMADfactors <- mvaMADs/mvaMAD;
      }

      #mvaMADthreshold <- outlierMAD * mvaMAD;
      names(mvaMADfactors) <- names(whichSamples);
      names(mvaMADs) <- names(whichSamples);
   }
   mvaMADoutliers <- x_names[whichSamples][which(mvaMADfactors >= outlierMAD)];
   if (verbose) {
      jamba::printDebug("mvaMADs:");
      print(format(digits=2, mvaMADs));
      jamba::printDebug("mvaMADfactors:");
      print(format(digits=2, mvaMADfactors));
      jamba::printDebug("mvaMADoutliers:", mvaMADoutliers);
   }
   attr(mvaDatas, "mvaMADs") <- mvaMADs;
   attr(mvaDatas, "mvaMADoutliers") <- mvaMADoutliers;
   attr(mvaDatas, "mvaMADfactors") <- mvaMADfactors;
   attr(mvaDatas, "MADs") <- mvaMADs;
   attr(mvaDatas, "MADoutliers") <- mvaMADoutliers;
   attr(mvaDatas, "MADfactors") <- mvaMADfactors;
   attr(mvaDatas, "outlierMAD") <- outlierMAD;
   attr(mvaDatas, "outlierRowMin") <- outlierRowMin;

   ## check_panel_page() checks if the panel number will cause a new
   ## page to be created, if maintitle is supplied then it will be
   ## printed atop each page
   check_panel_page <- function
   (iPanelNumber,
    maintitle,
    doHighlightLegend=FALSE,
    highlightColor=NULL,
    verbose=FALSE,
    ...) {
      if (length(maintitle) == 0 &&
            !(doHighlightLegend && length(highlightColor) > 0)) {
         if (verbose) {
            jamba::printDebug("check_panel_page(): ",
               "No extra title nor color legend is required.");
         }
         return();
      }
      max_panels <- prod(par("mfrow"));
      if (iPanelNumber > 1 && ((iPanelNumber-1) %% max_panels) == 0) {
         if (length(maintitle) > 0) {
            ## Print maintitle
            maintitle_txt <- paste(unlist(strsplit(maintitle, "\n")),
               collapse="\n");
            title(outer=TRUE,
               main=maintitle_txt,
               cex.main=maintitleCex);
         } else {
            if (verbose) {
               jamba::printDebug("check_panel_page(): ",
                  "No title was supplied.");
            }
         }
         if (doHighlightLegend && length(highlightColor) > 0) {
            ## print color legend
            if (verbose) {
               jamba::printDebug("highlightColor:");
               print(highlightColor);
            }
            outer_legend(x="bottom",
               legend=names(highlightColor),
               col=jamba::makeColorDarker(unlist(highlightColor)),
               pt.bg=unlist(highlightColor));
         } else {
            if (verbose) {
               jamba::printDebug("check_panel_page(): ",
                  "No highlightColor was supplied.");
            }
         }
      } else {
         if (verbose) {
            jamba::printDebug("check_panel_page(): ",
               "Panel did not force a page break.");
         }
      }
   }

   ## iPanelNumber keeps track of the numbered panels as they are plotted,
   ## so we can insert blank panels at the specified positions.
   if (doPlot) {
      iPanelNumber <- 0;
      for(i in whichSamples) {
         #mvaData <- mvaDatas[[i]];
         mvaData <- mvaDatas[[x_names[i]]];
         iPanelNumber <- iPanelNumber + 1;
         check_panel_page(iPanelNumber,
            maintitle,
            doHighlightLegend=doHighlightLegend,
            highlightColor=highlightColor,
            verbose=verbose,
            ...);
         if (verbose) {
            jamba::printDebug("iPanelNumber: ", iPanelNumber,
               fgText=c("orange", "lightgreen"));
         }
         if (length(blankPlotPos) > 0) {
            while(iPanelNumber %in% blankPlotPos) {
               jamba::nullPlot(doBoxes=FALSE);
               if (verbose) {
                  jamba::printDebug("   Inserted a blank panel at position: ",
                     iPanelNumber);
               }
               iPanelNumber <- iPanelNumber + 1;
               check_panel_page(iPanelNumber,
                  maintitle,
                  doHighlightLegend=doHighlightLegend,
                  highlightColor=highlightColor,
                  verbose=verbose,
                  ...);
               if (verbose) {
                  jamba::printDebug("new iPanelNumber: ",
                     iPanelNumber);
               }
            }
         }
         par("mar"=margins);
         groupName <- paste(c(x_names[i],
            groupSuffix[i]),
            collapse="");
         titleText <- x_names[i];

         ## Calculate the MAD, i.e. the median absolute deviation from zero
         mvaMAD <- mvaMADs[whichSamples[i]];
         #mvaMAD <- median(abs(mvaData[,"y"]), na.rm=TRUE);

         ## Determine color ramp per panel
         if (length(outlierMAD) > 0 && x_names[i] %in% mvaMADoutliers) {
            colrampUse <- colrampOutlier[[i]];
         } else {
            colrampUse <- colramp[[i]];
         }

         if (all(is.na(mvaData[,"y"]))) {
            jamba::nullPlot(doBoxes=FALSE,
               xlim=xlim,
               ylim=ylim
            );
            jamba::usrBox(fill=outlierColor)
            box();
         } else {
            mva <- smoothScatterFunc(mvaData[,c("x","y")],
               colramp=colrampUse(21),
               xlab="",
               ylab="",
               las=las,
               xlim=xlim,
               ylim=ylim,
               transformation=transformation,
               #col=smoothPtCol,
               useRaster=useRaster,
               nrpoints=nrpoints,
               fillBackground=fillBackground,
               applyRangeCeiling=applyRangeCeiling,
               ...);
         }
         ## Add axis labels
         title(xlab=xlab,
            line=xlabline);
         title(ylab=ylab,
            line=ylabline);

         ## Add some axis lines across the plot for easy visual reference
         if (!is.null(ablineH)) {
            if ("list" %in% class(ablineH)) {
               if (all(names(ablineH) %in% x_names)) {
                  h <- ablineH[[x_names[i]]];
               } else {
                  h <- ablineH[[i]];
               }
            } else {
               h <- unique(unlist(ablineH));
            }
            h <- h[h >= ylim[1] & h <= ylim[2]];
            if (length(h) > 0) {
               abline(h=h,
                  col="#44444488",
                  lty="dashed",
                  lwd=1,
                  ...);
            }
         }
         ## Add some axis lines across the plot for easy visual reference
         if (!is.null(ablineV)) {
            if ("list" %in% class(ablineV)) {
               if (all(names(ablineV) %in% x_names)) {
                  v <- ablineV[[x_names[i]]];
               } else {
                  v <- ablineV[[i]];
               }
            } else {
               v <- unique(unlist(ablineV));
            }
            v <- v[v >= xlim[1] & v <= xlim[2]];
            if (length(v) > 0) {
               abline(v=v,
                  col="#44444488",
                  lty="dashed",
                  lwd=1,
                  ...);
            }
         }

         ## Optionally highlight a subset of points
         if (length(highlightPoints) > 0) {
            hp1 <- lapply(seq_along(highlightPoints), function(highI){
               highP <- highlightPoints[[highI]];
               hiData <- mvaData[which(rownames(mvaData) %in% highP),,drop=FALSE];

               ## Make sure to restrict the y-values to fit within the axis limits,
               ## consistent with smoothScatterFunc().
               yValues <- jamba::noiseFloor(hiData[,"y"],
                  minimum=min(ylim),
                  ceiling=max(ylim),
                  ...);

               ## Optionally draw a polygon hull around highlighted points.
               if (doHighlightPolygon && length(yValues) > 2) {
                  hiM <- cbind(x=hiData[,"x"],
                     y=yValues);
                  hiHull <- hiM[grDevices::chull(hiM),,drop=FALSE];
                  #hiHull <- points2polygonHull(data.frame(x=hiData[,"x"],
                  #   y=yValues),
                  #   returnClass="matrix");
                  if (length(highlightPolygonAlpha) == 0) {
                     highlightPolyAlpha <- jamba::col2alpha(
                        highlightColor[[highI]]) / 2;
                  }
                  highlightPolyColor <- jamba::alpha2col(
                     highlightColor[[highI]],
                     alpha=highlightPolyAlpha);
                  polygon(hiHull,
                     add=TRUE,
                     col=highlightPolyColor,
                     border=highlightPolyColor);
               }
               ## Draw the highlighted points
               points(x=hiData[,"x"],
                  y=yValues,
                  pch=highlightPch[[highI]],
                  bg=highlightColor[[highI]],
                  col=jamba::makeColorDarker(darkFactor=1.2,
                     sFactor=1.5,
                     highlightColor[[highI]]),
                  xpd=FALSE,
                  cex=highlightCex[[highI]]);
            });
         }
         ## Optionally draw title boxes as labels per plot
         if (verbose) {
            jamba::printDebug("jammaplot(): ",
               "titleFont:",
               titleFont);
         }
         titleBoxTextColor <- titleColor[i];
         #titleBoxTextColor <- ifelse(colnames(object)[i] %in% mvaMADoutliers,
         #   blendColors(rep(c(titleColor[i], "red"), c(3,1))),
         #   titleColor[i]);
         if (doTitleBox) {
            if (verbose) {
               jamba::printDebug("jammaplot(): ",
                  "   titleBoxColor:",
                  titleBoxColor[i],
                  bgText=list(NA, titleBoxColor[i]),
                  fgText=list("orange", titleColor[i]));
            }
            parXpd <- par("xpd");
            par("xpd"=NA);
            ## Use drawLabels() for more control over the text box
            jamba::drawLabels(preset=titlePreset,
               txt=paste0(titleText, groupSuffix[i]),
               boxColor=titleBoxColor[i],
               boxBorderColor=jamba::makeColorDarker(titleBoxColor[i]),
               labelCol=titleColor[i],
               labelCex=titleCex[i],
               drawBox=doTitleBox,
               font=titleFont[i],
               ...);
            par("xpd"=parXpd);
         } else {
            titleLine <- margins[3] - 1.5;
            title(main="",
               sub=subtitle,
               cex.sub=titleCex[i],
               cex.main=titleCex[i],
               line=titleLine,
               col.main=titleBoxTextColor,
               font.main=titleFont[i],
               ...);
            title(main=paste(titleText, groupSuffix[i]),
               cex.main=titleCex[i],
               line=titleLine-1,
               col.main=titleBoxTextColor,
               font.main=titleFont[i],
               ...);
         }
         ## Subtitle with centered label
         if (length(subtitle) > 0) {
            ## Use drawLabels() for more control over the text box
            jamba::drawLabels(preset=subtitlePreset,
               txt=subtitle[i],
               boxColor=subtitleBoxColor[i],
               boxBorderColor=jamba::makeColorDarker(subtitleBoxColor[i]),
               labelCol=jamba::setTextContrastColor(subtitleBoxColor[i]),
               labelCex=titleCex[i]*0.9,
               drawBox=TRUE,
               font=titleFont[i],
               ...);
         }

         ## Optionally print the MAD factor in the bottom right corner
         if (length(outlierMAD) > 0 && displayMAD == 1) {
            jamba::drawLabels(preset="bottomright",
               labelCex=titleCex[i]*0.9,
               font=titleFont[i],
               txt=paste0("MAD x",
                  format(digits=2,
                     mvaMADfactors[whichSamples[i]])
               ),
               drawBox=FALSE,
               labelCol=ifelse(mvaMADfactors[whichSamples[i]] >= outlierMAD,
                  "red3", "grey40")
               );
         } else if (length(outlierMAD) > 0 && displayMAD == 2) {
            jamba::drawLabels(preset="bottomright",
               labelCex=titleCex[i]*0.9,
               font=titleFont[i],
               txt=paste0("MAD:",
                  format(digits=2,
                     mvaMADs[whichSamples[i]])),
               drawBox=FALSE,
               labelCol=ifelse(mvaMADfactors[whichSamples[i]] >= outlierMAD,
                  "red3", "grey40")
            );
         }
         #mvaData;
      }
      ## End of the per-panel MVA plot loop
      check_panel_page(iPanelNumber=prod(par("mfrow")) + 1,
         maintitle,
         doHighlightLegend=doHighlightLegend,
         highlightColor=highlightColor,
         verbose=verbose,
         ...);
   }

   invisible(mvaDatas);
}


#' Center gene data
#'
#' Performs row centering on a matrix of data in log space
#'
#' This function is a relatively simple wrapper function which subtracts
#' the row median (or row mean when mean=FALSE) from each row. The function
#' allows defining a subset of columns to be used in determining the
#' row control value via controlSamples, Similarly, columns can be grouped,
#' where columns are centered versus their relevant control samples within
#' each group of columns.
#'
#' @param x numeric matrix typically containing measurements (genes)
#'    as rows, and samples as columns.
#' @param floor optional numeric floor, below which values are set to the
#'    floor value, useful when one wants to avoid centering to values which
#'    may be below a noise threshold which might otherwise result in
#'    artificially inflated log fold changes.
#' @param controlSamples optional character vector of colnames(x) to be
#'    used as controls when centering data. In the event that centerGroups
#'    is also defined, the controlSamples are only used within each
#'    group of colnames. When this value is NULL or NA, or when any group
#'    of colnames defined by centerGroups contains no controlSamples, then
#'    all samples are used for centering. This relationship is clearly
#'    described in the attribute named "centerGroups".
#' @param centerGroups optional character vector named by colnames, whose
#'    values are group names. Alternatively, a list of vectors of colnames,
#'    where each list element contains colnames in explicit groups.
#' @param needsLod logical, indicating whether to perform log2 transformation
#'    of data prior to centering. If NULL, then if any value is above 40,
#'    it sets needsLog=TRUE and uses log2(x) for centering.
#' @param mean logical indicating whether to use row means, or row medians.
#'    If the matrixStats package is available, it uses
#'    `matrixStats::rowMedians()` for calculations, otherwise
#'    falling back to apply(x, 1, median) which is notably slower for
#'    large data matrices.
#' @param returnGroupedValues logical indicating whether to append columns
#'    which contain the control mean or median values used during centering.
#' @param showGroups logical indicating whether to print the sample centring
#'    relationship to screen during processing. Note this information is
#'    also contained in attribute "centerGroups".
#' @param scale character values indicating whether to scale data by row, or
#'    perform no row scaling. Scaling is dependent upon whether median or mean
#'    values are used in centering. If mean values are used, scaling is
#'    accomplished by dividing row values by the standard deviation. If median
#'    is used, then scaling divides row values by the MAD which is derived
#'    using the median instead of the mean.
#'
#' @return numeric matrix output, with the same row and column dimensions
#' as input data. If returnGroupedValues=TRUE, the additional columns contain
#' the row median or mean values, dependent upon mean=FALSE or mean=TRUE,
#' respectively.
#' An attribute \code{centerGroups} is included, which describes the specific
#' relationship between each colname, and associated control sample colnames,
#' and optional centerGroups grouping of colnames. When columns are grouped
#' and centered to specific control samples, is it important to keep this
#' information during downstream scrutiny of results.
#'
#' @family jam matrix functions
#'
#' @examples
#' x <- matrix(1:100, ncol=10);
#' colnames(x) <- letters[1:10];
#' # basic centering
#' centerGeneData(x);
#'
#' # grouped centering
#' centerGeneData(x,
#'    centerGroups=rep(c("A","B"), c(5,5)));
#'
#' # centering versus specific control columns
#' centerGeneData(x,
#'    controlSamples=letters[c(1:3)]);
#'
#' # grouped centering versus specific control columns
#' centerGeneData(x,
#'    centerGroups=rep(c("A","B"), c(5,5)),
#'    controlSamples=letters[c(1:3, 6:8)]);
#'
#' # confirm the centerGroups and controlSamples
#' x_ctr <- centerGeneData(x,
#'    centerGroups=rep(c("A","B"), c(5,5)),
#'    controlSamples=letters[c(1:3, 6:8)],
#'    showGroups=TRUE);
#' attr(x_ctr, "centerDF");
#'
#' @export
centerGeneData <- function
(x,
 floor=NA,
 controlSamples=NA,
 centerGroups=NULL,
 needsLog=NULL,
 mean=FALSE,
 returnGroupedValues=FALSE,
 returnValues=TRUE,
 groupPrefix="group",
 showGroups=FALSE,
 scale=c("none", "row"),
 verbose=FALSE,
 ...)
{
   ## Purpose is to provide quik data-centering (subtracting average from each row, usually
   ## in log2 space).
   ##
   ## If mean=FALSE, then median is used for centering, sometimes better than
   ## using mean because it is less sensitive to outliers
   ##
   ## New: if data contains text and numeric columns, center the numeric columns only
   ##
   ## centerGroups subdivides samples into groups.  Each group will
   ## be separately centered independent of the other groups.  Likewise, controlSamples
   ## are only used within each group, if the samples exist in each group, otherwise
   ## the group mean will be used for centering.
   ##
   ## centerGroups can be a list of vectors, containing colnames(x),
   ## or it can be a vector of group names, in order of colnames(x).
   ##
   ## If centerGroups is a named vector, and all colnames(x) are contained in named(centerGroups)
   ## then centerGroups[colnames(x)] will be used, to ensure they are both
   ## in consistent order, and to allow using a subset of x accordingly.
   ##
   ## showGroups=TRUE will by default print a quick data.frame summary showing
   ## samples, centerGroups, and controlSamples.
   ##
   ## returnGroupedValues=TRUE will append a column with the group values,
   ## either median or mean, as used in centering.
   ##
   ## scale is optionally intended to allow centered-scaled data, usually
   ## mean divided by standard deviation

   scale <- match.arg(scale);
   indata <- x;
   if (!returnValues && !returnGroupedValues) {
      stop("One must be TRUE: returnValues or returnGroupedValues.");
   }

   ## Determine whether the "matrixStats" R package is available
   if ("matrixStats" %in% rownames(installed.packages())) {
      rowMedians <- matrixStats::rowMedians;
      colMedians <- matrixStats::colMedians;
      rowSds <- matrixStats::rowSds;
      rowMads <- matrixStats::rowMads;
   } else {
      rowMedians <- function(x, na.rm=TRUE) {
         apply(x, 1, median, na.rm=na.rm);
      }
      colMedians <- function(x, na.rm=TRUE) {
         apply(x, 2, median, na.rm=na.rm);
      }
      rowSds <- function(x, na.rm=TRUE) {
         apply(x, 2, sd, na.rm=na.rm);
      }
      rowMads <- function(x, na.rm=TRUE) {
         apply(x, 2, mad, na.rm=na.rm);
      }
   }

   ## Separate character columns from numeric columns
   hasCharColumns <- FALSE;
   if (!any(class(indata) %in% c("matrix", "numeric"))) {
      colClass <- sapply(1:ncol(indata), function(i){class(indata[,i])});
      colClassChar <- which(!colClass %in% c("numeric", "integer"));
      if (length(colClassChar) > 0) {
         hasCharColumns <- TRUE;
         indataChar <- indata[,colClassChar, drop=FALSE];
         indata <- as.matrix(indata[,-colClassChar, drop=FALSE]);
      }
   }

   if (length(controlSamples) == 0 ||
       all(is.na(controlSamples))) {
      controls <- rep(TRUE, ncol(indata));
   } else {
      if (any(class(controlSamples) %in% c("integer", "numeric")) &&
            any(controlSamples <= ncol(indata))) {
         controls <- (seq_len(ncol(indata)) %in% controlSamples);
      } else if (length(colnames(indata)) > 0 &&
            any(controlSamples %in% colnames(indata))) {
         controls <- (colnames(indata) %in% controlSamples);
      } else {
         controls <- rep(TRUE, ncol(indata));
      }
      if (!any(controls)) {
         controls <- rep(TRUE, ncol(indata));
      }
   }
   if (length(colnames(indata)) > 0) {
      c_names <- colnames(indata);
      controlCols <- c_names[controls];
   } else {
      controlCols <- controls;
      c_names <- seq_len(ncol(indata));
   }

   ## We will try to auto-detect whether to log2 transform the data
   if (length(needsLog) == 0) {
      needsLog <- max(indata, na.rm=TRUE) > 100;
   }
   if (needsLog) {
      if (verbose) {
         jamba::printDebug("centerGeneData(): ",
            "applying log2 transform:",
            "log2(1 + x)");
      }
      indata <- log2(1 + indata);
   }

   ## Create summary data.frame
   centerDF <- as.data.frame(
      jamba::rmNULL(
         list(sample=c_names,
            centerGroups=centerGroups,
            centerControl=controls)));

   if (verbose) {
      jamba::printDebug("centerGeneData(): ",
         "dim(indata):",
         dim(indata));
      jamba::printDebug("centerGeneData(): ",
         "dim(centerDF):",
         dim(centerDF));
   }
   rownames(centerDF) <- centerDF[,"sample"];

   ## optionally show summary for visual confirmation
   if (showGroups) {
      print(centerDF);
   }

   ## Optionally center groups separately
   if (length(centerGroups) > 0) {
      if (!jamba::igrepHas("list", class(centerGroups))) {
         if (length(names(centerGroups)) > 0) {
            if (all(c_names %in% names(centerGroups))) {
               centerGroups <- centerGroups[c_names];
            } else {
               stop("colnames(indata) must be contained in names(centerGroups) when centerGroups has names.");
            }
         }
         centerGroups <- split(c_names, centerGroups);
      }
      ## Now iterate through the list
      centeredSubsets <- lapply(nameVectorN(centerGroups), function(iGroupN){
         iGroup <- centerGroups[[iGroupN]];
         iGroupCols <- colnames(indata[,iGroup,drop=FALSE]);
         iControls <- intersect(iGroupCols, controlCols);
         iM <- centerGeneData(indata[,iGroup,drop=FALSE],
            controlSamples=iControls,
            floor=floor,
            mean=mean,
            returnGroupedValues=returnGroupedValues,
            returnValues=returnValues,
            groupPrefix=iGroupN,
            needsLog=FALSE,
            showGroups=FALSE,
            centerGroups=NULL,
            scale=scale,
            ...);
         iM;
      });
      centeredData <- do.call(cbind, centeredSubsets);
      ## Correct rare cases when colnames are duplicated
      if (length(tcount(colnames(centeredData), minCount=2)) > 0) {
         colnames(centeredData) <- gsub("_v0$", "",
            makeNames(colnames(centeredData), startN=0));
      }
      centeredData <- centeredData[,unique(c(intersect(colnames(indata),
         colnames(centeredData)), colnames(centeredData))), drop=FALSE];
      attr(centeredData, "centerGroups") <- centerGroups;
      #return(centeredData);
   } else if (mean) {
      ## Note: Switched to using sweep() because it is much faster than apply()
      indataMeans <- rowMeans(indata[,controls, drop=FALSE],
         na.rm=TRUE);
      if (returnGroupedValues && !returnValues) {
         centeredData <- matrix(indataMeans,
            ncol=1,
            dimnames=list(rownames(indata), groupPrefix));
      } else {
         centeredData <- sweep(indata, 1, indataMeans);
         if (scale %in% "row") {
            indataSds <- rowSds(centeredData[,controls,drop=FALSE], na.rm=TRUE);
            ## Make sure no sd values are zero
            indataSds[indataSds == 0] <- mean(indataSds[indataSds != 0]);
            indataSds[indataSds == 0] <- 1;
            centeredData <- sweep(centeredData, 1, indataSds, "/")
         }
         if (returnGroupedValues) {
            centeredData <- cbind(centeredData, "mean"=indataMeans);
            if (length(groupPrefix) > 0 && nchar(groupPrefix) > 0) {
               n1 <- ncol(centeredData);
               colnames(centeredData)[n1] <- paste(
                  c(groupPrefix, colnames(centeredData)[n1]),
                  collapse="_");
            }
            if (scale %in% "row") {
               centeredData <- cbind(centeredData, "SD"=indataSds);
               if (length(groupPrefix) > 0 && nchar(groupPrefix) > 0) {
                  n1 <- ncol(centeredData);
                  colnames(centeredData)[n1] <- paste(
                     c(groupPrefix, colnames(centeredData)[n1]),
                     collapse="_");
               }
            }
         }
      }
   } else {
      ## Center by median
      indataMedians <- rowMedians(indata[,controls,drop=FALSE],
         na.rm=TRUE);
      if (returnGroupedValues && !returnValues) {
         centeredData <- matrix(indataMedians,
            ncol=1,
            dimnames=list(rownames(indata), groupPrefix));
      } else {
         centeredData <- sweep(indata, 1, indataMedians);
         if (scale %in% "row") {
            indataMads <- rowMads(centeredData[,controls, drop=FALSE],
               na.rm=TRUE);
            ## Make sure no sd values are zero
            indataMads[indataMads == 0] <- mean(indataMads[indataMads != 0]);
            indataMads[indataMads == 0] <- 1;
            centeredData <- sweep(centeredData, 1, indataMads, "/")
         }
         if (returnGroupedValues) {
            centeredData <- cbind(centeredData, "median"=indataMedians);
            if (length(groupPrefix) > 0 && nchar(groupPrefix) > 0) {
               n1 <- ncol(centeredData);
               colnames(centeredData)[n1] <- paste(
                  c(groupPrefix, colnames(centeredData)[n1]),
                  collapse="_");
            }
            if (scale %in% "row") {
               centeredData <- cbind(centeredData, "MAD"=indataMads);
               if (length(groupPrefix) > 0 && nchar(groupPrefix) > 0) {
                  n1 <- ncol(centeredData);
                  colnames(centeredData)[n1] <- paste(
                     c(groupPrefix, colnames(centeredData)[n1]),
                     collapse="_");
               }
            }
         }
      }
   }
   if (hasCharColumns) {
      if (verbose) {
         jamba::printDebug("centerGeneData(): ",
            "Keeping non-numeric columns as-is.");
      }
      centeredData <- cbind(indataChar, centeredData)
   }
   attr(centeredData, "centerDF") <- centerDF;
   return(centeredData);
}


#' define a polygon hull around points
#'
#' Define a polygon hull to encompass a set of data points
#'
#' A polygon hull can be useful as a visual indicator of the location
#' of a set of points. This function creates a polygon around all points,
#' and not just the most dense region of points, using
#' `grDevices::chull()`.
#'
#' The output is slightly modified to duplicate the first coordinate,
#' in order to "close" the resulting polygon, which is consistent
#' with the output of `rgeos::gConvexHull()`.
#'
#' @param x two-column numeric matrix of data points, or an object of
#'    class `sp::SpatialPoints`.
#' @param returnClass character string indicating the class to return,
#'    either `"matrix"` to return a two-column numeric matrix, or
#'    `"SpatialPolygons"` to return `sp::SpatialPolygons`.
#' @param verbose logical whether to print verbose output.
#'
#' @return
#' A two-column numeric matrix of polygon coordinates when
#' `returnClass="matrix"` by default, or when
#' `returnClass="SpatialPolygons"` it returns an object
#' of class `sp::SpatialPolygons`.
#'
#' @family jam plot functions
#'
#' @examples
#' set.seed(123);
#' xy <- cbind(x=1:10,
#'    y=sample(1:10, 10));
#' xyhull <- points2polygonHull(xy);
#' plot(xy, pch=20, cex=1, col="purple4",
#'    main="Polygon hull around points");
#' polygon(xyhull, border="purple2", col="transparent");
#'
#' @export
points2polygonHull <- function
(x,
 returnClass=c("matrix", "SpatialPolygons"),
 ...)
{
   ## Purpose is to take coordinates of a set of points, and return
   ## a polygon which includes all points inside it.
   ## If there are 3 or fewer points, the original points are returned.
   ##
   ## It will accept SpatialPoints or a two-column matrix,data.frame,tbl
   ## of x,y coordinates as input.
   ##
   ## It will return either a matrix of x,y coordinates, or SpatialPolygons
   ## containing one polygon (suitable for plotting directly) as output.
   ##
   returnClass <- match.arg(returnClass);
   x_class <- class(x);
   if (jamba::igrepHas("matrix|data.*frame|tbl", x_class)) {
      if (ncol(x) != 2) {
         stop(paste0("points2polygonHull() requires x with class",
            " matrix, data.frame, tibble, tbl",
            " to contain 2 columns."));
      }
   } else if (jamba::igrepHas("SpatialPoints", x_class)) {
      x <- SpatialPoints(x)@coords;
   }

   ## Determine convex hull coordinates
   x_hull <- x[grDevices::chull(x[,1:2]),,drop=FALSE];
   if (nrow(x_hull) > 2) {
      ## Add final row to "close" the polygon
      x_hull <- rbind(x_hull, x_hull[1,,drop=FALSE]);
   }
   if ("SpatialPolygons" %in% returnClass) {
      spp <- SpatialPolygons(list(
         Polygons(list(
            Polygon(x_hull)),
            "hull"
         )
      ));
      return(spp);
   }
   return(x_hull);
}

#' Center gene data (modified)
#'
#' Performs per-row centering on a numeric matrix
#'
#' This function centers data by subtracting the median or
#' mean for each row.
#'
#' Optionally columns can be grouped using the
#' argument `centerGroups`, which will center
#' data within each group independently.
#'
#' Data can be centered relative to specific control columns
#' using the argument `controlSamples`. When used with
#' `centerGroups`, each group of columns defined by
#' `centerGroups` that does not contain a corresponding
#' value in `controlSamples` will be centered using the
#' entire group.
#'
#' Confirm the `centerGroups` and `controlSamples` are
#' correct using the attribute `"center_df"` of the results,
#' see examples below.
#'
#' Note: This function assumes input data is log2-transformed,
#' or appropriately transformed to fit the assumption of
#' normality. This assumption is necessary for two reasons:
#'
#' 1. The group value (mean or median) is correct only when
#' data is transformed so the mean or median is not affected
#' by skewed data. Alternatively, `rowStatsFunc` can be
#' used to specify a custom group summary function.
#' 2. The centering subtracts the group value from each column
#' value.
#'
#' @param x numeric matrix of input data. See assumptions,
#'    that data is assumed to be log2-transformed, or otherwise
#'    appropriately transformed.
#' @param centerGroups character vector of group names, or
#'    `NULL` if there are no groups.
#' @param na.rm logical indicating whether NA values should be
#'    ignored for summary statistics. This argument is passed
#'    to the corresponding row stats function.
#' @param controlSamples character vector of values in `colnames(x)`
#'    which defines the columns to use when calculating group
#'    summary values.
#' @param useMedian logical indicating whether to use group median
#'    values when calculating summary statistics (`TRUE`), or
#'    group means (`FALSE`). In either case, when `rowStatsFunc`
#'    is provided, it is used instead.
#' @param rmOutliers logical indicating whether to perform outlier
#'    detection and removal prior to row group stats. This
#'    argument is passed to `jamba::rowGroupMeans()`. Note that
#'    outliers are only removed during the row group summary step,
#'    and not in the centered data.
#' @param madFactor numeric value passed to `jamba::rowGroupMeans()`,
#'    indicating the MAD factor threshold to use when `rmOutliers=TRUE`.
#'    The MAD of each row group is computed, the overall group median
#'    MAD is used to define 1x MAD factor, and any MAD more than
#'    `madFactor` times the group median MAD is considered an outlier
#'    and is removed. The remaining data is used to compute row
#'    group values.
#' @param rowStatsFunc optional function used to calculate row group
#'    summary values. This function should take a numeric matrix as
#'    input, and return a one-column numeric matrix as output, or
#'    a numeric vector with length `nrow(x)`. The function should
#'    also accept `na.rm` as an argument.
#' @param returnGroupedValues logical indicating whether to include
#'    the numeric matrix of row group values used during centering,
#'    returned in the attributes with name `"x_group"`.
#' @param returnGroups logical indicating whether to return the
#'    centering summary data.frame in attributes with name "center_df".
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to `jamba::rowGroupMeans()`.
#'
#' @family jam matrix functions
#'
#' @examples
#' x <- matrix(1:100, ncol=10);
#' colnames(x) <- letters[1:10];
#' # basic centering
#' centerGeneData_new(x);
#'
#' # grouped centering
#' centerGeneData_new(x,
#'    centerGroups=rep(c("A","B"), c(5,5)));
#'
#' # centering versus specific control columns
#' centerGeneData_new(x,
#'    controlSamples=letters[c(1:3)]);
#'
#' # grouped centering versus specific control columns
#' centerGeneData_new(x,
#'    centerGroups=rep(c("A","B"), c(5,5)),
#'    controlSamples=letters[c(1:3, 6:8)]);
#'
#' # confirm the centerGroups and controlSamples
#' x_ctr <- centerGeneData_new(x,
#'    centerGroups=rep(c("A","B"), c(5,5)),
#'    controlSamples=letters[c(1:3, 6:8)],
#'    returnGroups=TRUE);
#' attr(x_ctr, "center_df");
#'
#' @export
centerGeneData_new <- function
(x,
 centerGroups=NULL,
 na.rm=TRUE,
 controlSamples=NULL,
 useMedian=TRUE,
 rmOutliers=FALSE,
 madFactor=5,
 rowStatsFunc=NULL,
 returnGroupedValues=FALSE,
 returnGroups=FALSE,
 verbose=FALSE,
 ...)
{
   ## This function is a refactor of centerGeneData() to consolidate
   ## some logic into rowGroupMeans()
   if (length(x) == 0 || ncol(x) == 0 || nrow(x) == 0) {
      return(x);
   }

   ## Ensure that x has colnames
   if (length(colnames(x)) == 0) {
      colnames(x) <- seq_len(ncol(x));
   }

   ## Process controlSamples
   if (any(class(controlSamples) %in% c("numeric", "integer"))) {
      ## Numeric controlSamples are used to subset colnames(x)
      if (!"integer" %in% class(controlSamples) &&
            !all(controlSamples == as.integer(controlSamples))) {
         stop("controlSamples must be integer or numeric integer values, decimal values were detected.");
      }
      controlSamples <- colnames(x)[seq_len(ncol(x)) %in% controlSamples];
   } else {
      controlSamples <- intersect(controlSamples, colnames(x));
   }
   if (length(controlSamples) == 0) {
      if (verbose) {
         jamba::printDebug("centerGeneData_new(): ",
            "controlSamples using all colnames(x)");
      }
      controlSamples <- colnames(x);
   }

   ## Process centerGroups
   if (length(centerGroups) == 0) {
      centerGroups <- rep("Group", ncol(x));
   } else if (length(centerGroups) == 1) {
      ## Note that NA is converted to "NA"
      centerGroups <- rep(
         jamba::rmNA(centerGroups,
            naValue="NA"),
         length.out=ncol(x));
   }
   if (length(centerGroups) < ncol(x)) {
      stop("length(centerGroups) must equal ncol(x), or be length 0 or 1 to indicate no groups.");
   }
   names(centerGroups) <- colnames(x);

   ## Confirm all centerGroups contains controlSamples,
   ## use all samples when any centerGroup has no controlSamples
   samples_l <- split(colnames(x),
      centerGroups);
   controls_l <- lapply(samples_l, function(i){
      if (!any(i %in% controlSamples)) {
         i;
      } else {
         intersect(i, controlSamples);
      }
   })
   controls_v <- unname(unlist(controls_l));
   center_df <- data.frame(sample=colnames(x),
      centerGroups=centerGroups,
      controlSamples=colnames(x) %in% controls_v);

   ## Calculate row summary values
   x_group <- jamba::rowGroupMeans(x[,controls_v, drop=FALSE],
      na.rm=na.rm,
      groups=centerGroups[controls_v],
      useMedian=useMedian,
      rmOutliers=rmOutliers,
      madFactor=madFactor,
      rowStatsFunc=rowStatsFunc,
      verbose=verbose,
      ...);

   ## Now produce centered values by subtracting the group summary values
   x_centered <- (x - x_group[,centerGroups[colnames(x)], drop=FALSE]);

   if (returnGroupedValues) {
      attr(x_centered, "x_group") <- x_group;
   }
   if (returnGroups) {
      attr(x_centered, "center_df") <- center_df;
   }

   return(x_centered);
}

#' Calculate MA-plot data
#'
#' Calculate MA-plot data
#'
#' This function takes a numeric matrix as input, and calculates
#' data sufficient to produce MA-plots. The default output is a
#' list of two-column numeric matrices with `"x"` and `"y"` coordinates,
#' representing the group median and difference from median,
#' respectively.
#'
#' The mean value can be used by setting `useMedian=FALSE`.
#'
#' Samples can be grouped using the argument `centerGroups`.
#' In this case the y-axis value will be "difference from
#' group median."
#'
#' Control samples can be specified for centering using the
#' argument `controlSamples`. In this case, the y-axis value will
#' be "difference from control median".
#'
#' The sample grouping, and control samples can be combined,
#' in which case the y-axis values will be "difference from
#' the control median within the centering group."
#'
#' @family jam matrix functions
#'
#' @param x numeric matrix typically containing log-normal measurements,
#'    with measurement rows, and sample columns.
#' @param na.rm logical indicating whether to ignore NA values
#'    during numeric summary functions.
#' @param controlSamples character vector containing values in
#'    `colnames(x)` to define control samples used during centering.
#'    These values are passed to `centerGeneData()`.
#' @param centerGroups character vector with length equal to `ncol(x)`
#'    which defines the group for each column in `x`. Data will
#'    be centered within each group.
#' @param groupedX logical indicating how to calculate the x-axis
#'    value when `centerGroups` contains multiple groups. When
#'    groupedX=TRUE, the mean of each group median is used, which
#'    has the effect of representing each group equally. When
#'    groupedX=FALSE, the median across all columns is used, which
#'    can have the effect of preferring sample groups with a larger
#'    number of columns.
#' @param useMedian logical indicating whether to use the median
#'    values when calculating the x-axis and during data centering.
#'    The median naturally reduces the effect of outlier points on
#'    the resulting MA-plots., when compared to using the mean.
#'    When useMedian=FALSE, the mean value is used.
#' @param useMean (deprecated) logical indicating whether to use the
#'    mean instead of the median value. This argument is being removed
#'    in order to improve consistency with other Jam package functions.
#' @param whichSamples character vector containing `colnames(x)`, or
#'    integer vector referencing column numbers in `x`. This argument
#'    specifies which columns to return, but does not change the columns
#'    used to define the group centering values. For example, the
#'    group medians are calculated using all the data, but only the
#'    samples in `whichSamples` are centered to produce MA-plot data.
#' @param noise_floor numeric value indicating the minimum numeric value
#'    allowed in the input matrix `x`. When `NULL` or `-Inf` no noise
#'    floor is applied. It is common to set `noise_floor=0` to limit
#'    MA-plot data to use values zero and above.
#' @param noise_floor_value single value used to replace numeric values
#'    at or below `noise_floor` when `noise_floor` is not NULL. By default,
#'    `noise_floor_value=noise_floor` which means values at or below
#'    the noise floor are set to the floor. Another useful option is
#'    `noise_floor_value=NA` which has the effect of removing the point
#'    from the MA-plot altogether. This option is recommended for sparse
#'    data matrices where the presence of values at or below zero are
#'    indicative of missing data (zero-inflated data) and does not
#'    automatically reflect an actual value of zero.
#' @param naValue single value used to replace any `NA` values in
#'    the input matrix `x`. This argument can be useful to replace
#'    `NA` values with something like zero.
#' @param mad_row_min numeric value defining the minimum group
#'    value, corresponding to the x-axis position on the MA-plot,
#'    required for a row to be included in the MAD calculation.
#'    This threshold is useful to filter outlier data below a noise
#'    threshold, so that the MAD calculation will include only the
#'    data above that value. For example, with count data, it is
#'    useful to filter out counts below roughly 8, where Poisson
#'    noise is a more dominant component than real count data.
#'    Remember that count data should already be log2-transformed,
#'    so the threshold should also be identically transformed,
#'    for example using `log2(1 + 8)` to set a minimum count
#'    threshold of at least 8.
#' @param grouped_mad logical indicating whether the MAD value
#'    should be calculated per group when `centerGroups` is
#'    supplied, from which the MAD factor values are derived.
#'    When `TRUE` it has the effect of highlighting outliers
#'    within each group using the variability in that group.
#'    When `FALSE` the overall MAD is calculated, and a
#'    particularly high variability group may have all its
#'    group members labeled with a high MAD factor.
#' @param centerFunc function used for centering data, by default
#'    one of the functions `centerGeneData()` or `centerGeneData_new()`.
#'    This argument will be removed in the near future and is mainly
#'    intended to allow testing the two centering functions.
#' @param returnType character string indicating the format of data
#'    to return: `"ma_list"` is a list of MA-plot two-column
#'    numeric matrices with colnames `c("x","y")`; "tidy"
#'    returns a tall `data.frame` suitable for use in ggplot2.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#'
#' @export
jammacalc <- function
(x,
 na.rm=TRUE,
 controlSamples=NULL,
 centerGroups=NULL,
 groupedX=TRUE,
 useMedian=TRUE,
 useMean=NULL,
 whichSamples=NULL,
 noise_floor=-Inf,
 noise_floor_value=noise_floor,
 naValue=NA,
 mad_row_min=0,
 grouped_mad=TRUE,
 centerFunc=centerGeneData_new,
 useRank=FALSE,
 returnType=c("ma_list", "tidy"),
 verbose=FALSE,
 ...)
{
   ## Purpose is to separate jammaplot() from the underlying
   ## math that produces the data for MA-plots.
   if (length(x) == 0 || ncol(x) == 0) {
      return(NULL);
   }
   returnType <- match.arg(returnType);
   if (length(useMean) > 0) {
      useMedian <- !useMean;
   }

   ## Validate whichSamples
   if (length(whichSamples) == 0) {
      whichSamples <- seq_len(ncol(x));
   } else if (any(c("numeric","integer") %in% class(whichSamples))) {
      whichSamples <- whichSamples[whichSamples %in% seq_len(ncol(x))];
      if (length(whichSamples) == 0) {
         whichSamples <- seq_len(ncol(x));
      }
   } else {
      whichSamples <- whichSamples[whichSamples %in% colnames(x)];
   }
   if (verbose) {
      jamba::printDebug("jammacalc(): ",
         "whichSamples:",
         whichSamples);
   }

   if (!suppressPackageStartupMessages(require(matrixStats))) {
      rowMedians <- function(x, na.rm=TRUE) {
         apply(x, 1, median, na.rm=na.rm);
      }
      colMedians <- function(x, na.rm=TRUE) {
         apply(x, 2, median, na.rm=na.rm);
      }
   }

   ## Apply noise_floor if needed
   noise_floor <- head(noise_floor, 1);
   noise_floor_value <- head(noise_floor_value, 1);
   if (length(noise_floor) == 1 && !is.infinite(noise_floor)) {
      if (all(noise_floor %in% noise_floor_value)) {
         if (jamba::rmNA(naValue=0, any(x < noise_floor))) {
            if (verbose) {
               jamba::printDebug("jammacalc(): ",
                  "flooring values below ",
                  noise_floor,
                  " to ",
                  noise_floor_value);
            }
            x[x < noise_floor] <- noise_floor_value;
         }
      } else if (jamba::rmNA(naValue=0, any(x <= noise_floor))) {
         if (verbose) {
            jamba::printDebug("jammacalc(): ",
               "flooring values at or below ",
               noise_floor,
               " to ",
               head(noise_floor_value, 10));
         }
         x[x <= noise_floor] <- noise_floor_value;
      }
   }
   ## Replace NA if needed
   if (!all(naValue %in% c(NA))) {
      if (any(is.na(x))) {
         if (verbose) {
            jamba::printDebug("jammacalc(): ",
               "replacing NA with ",
               naValue);
         }
         x[is.na(x)] <- naValue;
      }
   }

   ## apply useRank
   if (useRank) {
      x1 <- apply(x, 2, rank, na.last=FALSE);
      if (any(is.na(x))) {
         x1[is.na(x)] <- NA;
      }
      x <- x1;
      rm(x1);
   }

   ## Center the data in one step
   x_ctr <- centerFunc(x,
      na.rm=na.rm,
      controlSamples=controlSamples,
      centerGroups=centerGroups,
      useMedian=useMedian,
      returnGroups=TRUE,
      returnGroupedValues=TRUE,
      ...);
   rownames(x_ctr) <- rownames(x);

   ## Pull out grouped data
   if ("x_group" %in% names(attributes(x_ctr))) {
      x_grp <- attr(x_ctr, "x_group");
   } else {
      if (length(centerGroups) == 0) {
         grp_col <- vigrep("^group_me", colnames(x_ctr));
         x_grp <- x_ctr[,grp_col,drop=FALSE];
      } else if (all(centerGroups %in% colnames(x_ctr))) {
         grp_col <- provigrep(paste0("^", unique(centerGroups), "_me"), colnames(x_ctr));
         x_grp <- x_ctr[,grp_col,drop=FALSE];
      } else {
         stop("The group values could not be detected from centerFunc().");
      }
   }
   rownames(x_grp) <- rownames(x);

   ## Determine the appropriate x-axis value to use
   if (groupedX && ncol(x_grp) > 1 && length(unique(centerGroups)) > 1) {
      if (verbose) {
         jamba::printDebug("jammacalc(): ",
            "Using group values for MA-plot x-axis values");
      }
   } else {
      groupedX <- FALSE;
      if (useMedian) {
         if (verbose) {
            jamba::printDebug("jammacalc(): ",
               "using rowMedian() of all values for x.");
         }
         x_use <- rowMedians(x,
            na.rm=na.rm);
      } else {
         if (verbose) {
            jamba::printDebug("jammacalc(): ",
               "using rowMean() of all values for x.");
         }
         x_use <- rowMeans(x,
            na.rm=na.rm);
      }
      names(x_use) <- rownames(x);
   }

   ## Calculate per-column MAD values, the grouped MAD values
   ## First create a matrix of centered values
   x_mad_m <- as.matrix(x_ctr);
   ## blank entries where the x-axis value is not above mad_row_min
   if (groupedX) {
      x_mad_m[x_grp[,centerGroups] < mad_row_min] <- NA;
   } else {
      x_mad_m[x_use < mad_row_min] <- NA;
   }
   ## Compute MAD values for all columns
   x_mads <- jamba::nameVector(
      colMedians(abs(x_mad_m),
         na.rm=TRUE),
      colnames(x_ctr));
   if (grouped_mad && length(unique(centerGroups)) > 1) {
      if (verbose) {
         jamba::printDebug("jammacalc(): ",
            "Calculating grouped MAD values.");
      }
      x_grp_mads <- tapply(x_mads, centerGroups, median, na.rm=TRUE);
      x_mad_factors <- x_mads / x_grp_mads[as.character(centerGroups)];
      names(x_mad_factors) <- names(x_mads);
   } else {
      if (verbose) {
         jamba::printDebug("jammacalc(): ",
            "Calculating global MAD values.");
      }
      x_grp_mads <- median(x_mads, na.rm=TRUE);
      names(x_grp_mads) <- "global";
      x_mad_factors <- x_mads / x_grp_mads;
      names(x_mad_factors) <- names(x_mads);
   }

   ## Calculate the list of two-column matrices
   jammadata <- lapply(whichSamples, function(i){
      if (groupedX) {
         cbind(x=x_grp[,as.character(centerGroups[i])],
            y=x_ctr[,i]);
      } else {
         cbind(x=x_use,
            y=x_ctr[,i]);
      }
   });
   names(jammadata) <- colnames(x)[whichSamples];

   if ("tidy" %in% returnType) {
      jammadata <- rbindList(lapply(names(jammadata), function(i){
         j <- jammadata[[i]];
         data.frame(row=rownames(j),
            sample=i,
            colnum=match(i, colnames(x)),
            x=j[,"x"],
            y=j[,"y"]
         )
      }));
      ## Suitable for
      ## ggplot2::ggplot(jammmadata_df, ggplot2::aes(x=x, y=y)) +
      ## ggplot2::geom_point() + ggplot2::facet_wrap(~sample) + colorjam::theme_jam()
   }

   ## Add MAD data to jammadata
   attr(jammadata, "MADs") <- x_mads;
   attr(jammadata, "groupMADs") <- x_grp_mads;
   attr(jammadata, "MADfactors") <- x_mad_factors;

   return(jammadata);
}
