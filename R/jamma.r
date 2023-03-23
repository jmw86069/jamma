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
#' The jamma package creates MA-plots for omics data, and provides
#' important options to handle specific experiment designs and
#' strategies for data quality control.
#'
#' ## Overview of features:
#'
#' * MA-plots can be calculated using the mean or median signal.
#' * Data can be centered using a subset of reference samples.
#' * Data can be centered within groups of samples, useful to
#' assess within-group variability, or within-batch variability.
#' * Ranked MA-plots can be generated to show rank-difference,
#' useful to assess consistency of the rank ordered signal across
#' samples.
#' * Putative technical outliers can be defined using a MAD factor
#' threshold derived from the data itself, to highlight individual
#' samples with much higher variability than expected from biological
#' sources, which often highlight technical failures in upstream protocol.
#'
#' ## Data centering
#'
#' For example, it can be useful to generate MA-plots within
#' biological sample replicates, or even among technical replicates.
#' By this approach, MA-plots can effectively highlight technical
#' outliers, where variability in one sample is measurably
#' higher than that from other comparable samples. A MAD outlier
#' approach is available to identify samples whose median
#' variance is more than X times higher than that across other
#' samples.
#'
#' It is useful to center within sample types, for example brain
#' samples can be centered independently of kidney or liver
#' samples. This approach is especially useful when statistical
#' comparisons are not intended to be applied across brain
#' and kidney for example.
#'
#' In general, it is recommended to use
#' `centerGroups` to center data within meaningful experimental
#' subsets where there are no intended statistical
#' comparisons across these subsets.
#' We find it useful to generate MA-plots across all samples
#' even when there are distinct experimental subsets, because it provides
#' context to the signal profiles obtained overall.
#' For example it may be informative to recognize that signal from one
#' experimental subset is lower and/or more variable than signal
#' from another subset. It could be of biological or technical
#' importance.
#'
#' ## Data Normalization
#'
#' Lastly, the MA-plot approach is often effective at visualizing
#' the need for data normalization, which is equivalent to methods
#' such as log-ratio normalization. The underlying assumption is that
#' the median or mean log ratio (y-axis difference shown on MA-plots)
#' is zero.
#'
#' A normalization method `jammanorm()` provides this normalization.
#' Note that it also abides by the `centerGroups` and `controlSamples`
#' arguments. Additional argument `controlGenes` optionally defines
#' a specific subset of genes as normalizers, equivalent to using
#' housekeeper genes for normalization. Note that housekeeper normalization
#' in this case is defined by housekeeper genes having log ratio of zero,
#' and does not directly use the geometric mean expression of housekeepers,
#' although the result is very often nearly identical.
#'
#' ## Volcano plots
#'
#' Volcano plots are similar to MA-plots, with some useful distinctions:
#' * Volcano plots display **group** log fold change results versus P-value,
#' based upon a statistical test.
#' * MA-plots display **per-sample** log differences
#' from control, versus the mean signal. Often the P-value is related to
#' the mean signal, therefore these plots have some resemblance.
#' * It is possible to show **group** MA-plots, notably `DESeq2::plotMA()`,
#' although its purpose is to display grouped summary to indicate
#' the effect of signal on the fold change threshold for statistical
#' significance. It is not intended to assess consistent signal across
#' individual samples.
#'
#' @family jamma package
#'
#' ## Core functions:
#'
#' * `jammaplot()`
#' * `centerGeneData()`
#' * `jammanorm()`
#' * `jammacalc()`
#' * `volcano_plot()`
#'
#' ## Additional plot functions:
#' * `points2polygonHull()`
#' * `outer_legend()`
#'
#' @docType package
#' @name jamma
NULL


#' Produce MA-plot of omics data.
#'
#' Produce MA-plot of omics data, where `jammaplot()` uses base R graphics,
#' `ggjammaplot()` uses ggplot2 graphics.
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
#' The argument `noise_floor` provides a numeric lower threshold,
#' where individual values **at or below** this threshold are
#' set to a defined value, defined by argument `noise_floor_value`.
#' The default was updated in version `0.0.21.900` to
#' `noise_floor=0` and `noise_floor_value=NA`.
#' Values of zero `0` are set to `NA` and therefore are not included
#' in the MA-plot calculations. Only points above zero are included
#' as points in each MA-plot panel.
#'
#' Another useful alternative is to define `noise_floor_value=noise_floor`
#' which sets any measurement **at or below** the `noise_floor` to
#' this value. This option has the effect of reducing random noise from
#' points that are already below the noise threshold and therefore
#' are unreliable for this purpose.
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
#'       classes: `ExpressionSet`, `SummarizedExperiment`,
#'       `MultiExperimentSet`}
#'    \item{Make it efficient to convey group information, for example
#'       define `titleBoxColor` with group colors, allow `centerByGroup=TRUE`
#'       which would re-use known sample group information.}
#'    \item{Adjust the suffix to indicate when \code{centerGroups} are being
#'       used. For example indicate `'sampleID vs groupA'` instead of
#'       `'sampleID vs median'`.}
#' }
#'
#' @param x `numeric` object usually a `matrix` that contains
#'    values with measurement rows, and sample/observation columns.
#'    For example, with gene or protein expression data, the genes
#'    or proteins (or the assays of genes or proteins) are
#'    represented in rows, and obtained samples are represented
#'    in columns. Alternatively `x` can be `SummarizedExperiment`
#'    object, used alongside argument `assay_name`.
#' @param assay_name `character` used when `x` is a
#'    `SummarizedExperiment` object, to determine which assay
#'    `matrix` to use for the MA plots. When `assay_name=NULL`
#'    the first assay entry is used, for example `assays(x)[[1]]`.
#' @param maintitle `character` string with the title displayed above
#'    all individual MA-plot panels. It will appear in the top outer
#'    margin.
#' @param titleBoxColor,subtitleBoxColor `character` vector of
#'    R colors used as background color for each panel title text,
#'    or subtitle text respectively. The subtitle appears in the
#'    bottom-left corner, and usually indicates the center groups
#'    as defined by `centerGroups`.
#' @param centerGroups `character` vector of groups passed to
#'    `jamma::centerGeneData()` which determines how data is centered.
#'    Each group is centered independently, to enable visual
#'    comparisons within each relevant centering group.
#'    It is useful to center within batches or within
#'    subsets of samples that are not intended to be compared to
#'    one another.
#'    Another useful alternative is to center by each sample
#'    group in order to view the variability among group
#'    replicates, which should be much lower than variability across
#'    sample groups. See `centerGeneData()` for more specific examples.
#' @param controlSamples `character` vector of `colnames(x)` passed
#'    to `centerGeneData()` which defines the control samples during
#'    the data centering step.
#'    By default, and the most common practice, MA-plots are
#'    calculated across all samples, which effectively uses all
#'    `colnames(x)` as `controlSamples`.
#'    However, it is quite useful sometimes to provide a subset of
#'    samples especially if there are known quality samples, to which
#'    new samples of unknown quality are being compared.
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
#' @param maintitleCex `numeric` cex character expansion used to
#'    resize the `maintitle`.
#' @param subtitle `NULL` or `character` vector to be drawn at
#'    the bottom left corner of each plot panel, the location
#'    is defined by `subtitlePreset`.
#' @param subtitlePreset `character` string passed to `jamba::coordPresets()`.
#'    The default `subtitlePreset="bottomleft"` defines the bottom-left
#'    corner of each panel.
#' @param subtitleAdjPreset `character` string passed to `jamba::coordPresets()`.
#'    The default `subtitleAdjPreset="topright"` places labels to the
#'    top-right of the subtitle position, which by default is the bottom-left
#'    corner of each panel.
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
#' @param titlePreset `character` string passed to `jamba::coordPresets()`.
#'    The default `titlePreset="top"` defines the top edge of each panel.
#' @param titleAdjPreset `character` string passed to `jamba::coordPresets()`.
#'    The default `titleAdjPreset="top"` places labels above the `titlePreset`
#'    location, by default above the top edge of each panel.
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
#' @param useMedian `logical` indicates whether to center data
#'    using the `median` value, where `useMedian=FALSE` by default.
#'    For consistency, this argument is preferred to `useMean` which
#'    is deprecated and will be removed in future. The median is
#'    preferred in cases where outliers should not influence the
#'    centering. The mean is preferred in cases where the data
#'    should visualize data consistent with downstream parametric
#'    statistical analysis. When a particular sample
#'    is a technical outlier, one option is to define
#'    `controlSamples` to exclude the outlier sample(s), so
#'    the data centering will be applied using the non-outlier
#'    samples as reference.
#' @param useMean `logical` (deprecated), use `useMean`. This argument
#'    indicates whether to center data using the `mean` value.
#'    When `useMean=NULL` the argument `useMedian` is preferred.
#'    For backward compatibility, when `useMean` is not `NULL`,
#'    then `useMedian` is defined by `useMedian <- !useMean`.
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
#' @param noise_floor,noise_floor_value `numeric` to define a numeric
#'    floor, or `NULL` for no numeric floor. Values **at or below**
#'    `noise_floor` are set to `noise_floor_value`, intended for two
#'    potential uses:
#'    1. Filter out value below a threshold, so they do not affect centering.
#'       * This option is valuable to remove zeros when a zero `0` is considered
#'       "no measurement observed", typically for count data such as RNA-seq,
#'       NanoString, and especially single-cell protocols or other protocols
#'       that produce a large number of missing values.
#'       * One can typically tell whether input data includes zero `0`
#'       values by the presence of characteristic 45-degree angle lines
#'       originating from `x=0` angled toward the right. The points along
#'       this line are rows with more measurements of zero than non-zero,
#'       there this sample has a non-zero value.
#'
#'    2. Set values at a noise floor to the noise floor, to retain the
#'    measurement but minimize the effect during centering to the lowest
#'    realiable measurement for the platform technology.
#'       * This value may be set to a platform noise floor
#'       for something like microarray data where the intensity may be
#'       unreliable below a threshold; or
#'       * for quantitative PCR measurements where cycle threshold (Ct)
#'       values may become unreliable, for example above CT=40 or CT=35.
#'       Data is often transformed to abundance with `2 ^ (40 - CT)` then
#'       log2-transformed for analysis. In this case, to apply a `noise_floor`
#'       effective for CT=35, one would use `noise_floor=5`.
#' @param filterFloor,filterFloorReplacement (deprecated) in favor of
#'    `noise_floor`, and `noise_floor_replacement` respectively.
#' @param transFactor `numeric` adjustment to the visual density of
#'    smooth scatter points. For base R graphics, this argument is
#'    passed to `jamba::plotSmoothScatter()`. The argument value is based upon
#'    `graphics::smoothScatter()` argument `transformation`, which uses
#'    default `function(x)x^0.25`. The `transFactor` is equivalent to the
#'    exponential in the form: `function(x)x^transFactor`. Lower values
#'    make the point density more visually intense, higher values make the
#'    point density less visually intense.
#' @param nrpoints `integer` or `NULL` indicating the number of points
#'    to display on the extremity of the smooth scatter density,
#'    passed to `jamba::plotSmoothScatter()`.
#' @param smoothScatterFunc `function` used to produce a smooth scatter plot
#'    in base R graphics. The default `jamba::plotSmoothScatter()` controls
#'    the level of detail in the density calculation, and in the graphical
#'    resolution of that density in each plot panel. The custom function
#'    should accept argument `transformation` as described in `transFactor`,
#'    even if the argument is not used.
#' @param ablineH,ablineV `numeric` vector indicating position of
#'    horizontal and vertical lines in each MA-plot panel.
#' @param doTxtplot `logical` (not yet implemented in `jamma`),
#'    indicating to produce colored ANSI text plot output, for example
#'    to a text terminal.
#' @param blankPlotPos `NULL` or `integer` vector indicating
#'    plot panel positions to be drawn blank, and therefore skipped.
#'    Plot panels are drawn in the exact order of `colnames(x)` received.
#'    Blank panel positions are intended to help customize the visual
#'    alignment of MA-plot panels. The mechanism is similar to
#'    `ggplot2::facet_wrap()` except that blank positions can be manually
#'    defined by what makes sense for the experiment design.
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
#' @param groupedMAD `logical` indicating how the MAD calculation
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
#' @param fillBackground `logical` currently used for base R graphics
#'    output, and passed to `jamba::plotSmoothScatter()`,
#'    indicating whether to fill the plot panel using the
#'    first color in the color ramp for each MA-plot panel, or when
#'    a plot panel is an outlier, it uses `outlierColor`.
#'    This argument is mainly useful to highlight outlier panels,
#'    although it is also useful when the color ramp has non-white
#'    base color, for example `viridis::viridis()`.
#' @param ma_method character string indicating whether to perform
#'    MA-plot calculations using the old method `"old"`; or `"jammacalc"`
#'    which uses the function `jammacalc()`.
#' @param panel_hook_function `function` or `NULL`, with custom function
#'    to be called as a "hook" after each MA-plot panel has been drawn.
#'    This `panel_hook_function` is recycled to the number of samples,
#'    defined internally as `nsamples`, although the same function can
#'    be called on all panels. If `panel_hook_function` is supplied as
#'    a `list`, it is recycled to the number of samples `nsamples`,
#'    and any element in the list which has `length==0` or is `NA` will
#'    not be called as a function.
#'    This function should accept at least two arguments:
#'    * `i` - an `integer` indicating the sample to be plotted in order,
#'    as defined by `colnames(x)`, and `whichSamples` in the event samples
#'    are subsetted or re-ordered with `whichSamples`.
#'    * `...` additional arguments passed by `...` in this function.
#'    * Any arguments of `jammaplot()` are available inside the panel
#'    hook function as a by-product of calling this function within
#'    the environment of the active `jammaplot()`, therefore any
#'    argument values will be available for use inside that function.
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
#' @return `list` of `numeric` `matrix` objects, one for each MA-plot,
#'   with colnames `"x"` and `"y"`. This `list` is sufficient input
#'   to `jammaplot()` to re-create the full set of MA-plots.
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
#'
#'    jammaplot(edata,
#'       whichSamples=c(1, 2));
#'
#'    jammaplot(edata,
#'       sample_labels=paste("Sample", colnames(edata)));
#'
#'    jammaplot(edata,
#'       controlIndicator="titlestar");
#'
#'    jammaplot(edata,
#'       controlIndicator="none");
#'
#'    jammaplot(edata,
#'       panel_hook_function=function(i,...){box("figure")});
#'
#'    jammaplot(edata,
#'       useRank=TRUE,
#'       maintitle="Rank MA-plots");
#' }
#'
#' @export
jammaplot <- function
(x,
 assay_name=NULL,
 maintitle=NULL,
 titleBoxColor="#DDBB9977",
 subtitleBoxColor=titleBoxColor,
 centerGroups=NULL,
 controlSamples=colnames(x),
 controlFloor=NA,
 naControlAction=c("row", "floor", "min", "na"),
 naControlFloor=0,
 controlIndicator=c(
    "labelstar",
    "titlestar",
    "none"),
 sample_labels=NULL,
 useMedian=FALSE,
 useMean=NULL,
 ylim=c(-4,4),
 xlim=NULL,
 highlightPoints=NULL,
 outlierMAD=5,
 outlierRowMin=5,
 displayMAD=FALSE,
 groupedMAD=TRUE,
 colramp=c("white", "lightblue", "blue", "navy", "orange", "orangered2"),
 colrampOutlier=NULL,
 outlierColor="lemonchiffon",
 whichSamples=NULL,
 maintitleCex=1.8,
 subtitle=NULL,
 subtitlePreset="bottomleft",
 subtitleAdjPreset="topright",
 titleCexFactor=1,
 titleCex=NULL,
 doTitleBox=TRUE,
 titleColor="black",
 titleFont=2,
 titlePreset="top",
 titleAdjPreset="top",
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
 margins=c(2.5, 2.5, 2.0, 0.2),
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
 noise_floor=0,
 noise_floor_value=NA,
 filterFloor=NULL,
 filterFloorReplacement=NULL,
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
 panel_hook_function=NULL,
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
   if (jamba::check_pkg_installed("matrixStats")) {
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
   controlIndicator <- match.arg(controlIndicator);
   naControlAction <- match.arg(naControlAction);

   if (doTxtplot) {
      blankPlotPos <- NULL;
      smoothScatterFunc <- function(...){
         plotTextSmoothScatter(height=20,
            width=80,
            doLegend=FALSE,
            ...);
      }
   }
   if (length(useMean) > 0 && is.logical(useMean)) {
      useMedian <- !useMean;
      if (verbose) {
         jamba::printDebug("jammaplot(): ",
            "useMedian defined by !useMean, useMedian=",
            useMedian);
      }
   }
   if (length(useMedian) == 0) {
      useMedian <- FALSE;
   }

   transformation <- function(x){
      x^transFactor;
   }
   if ("list" %in% class(x)) {
      nsamples <- length(x);
      x_names <- names(x);
   } else if ("SummarizedExperiment" %in% class(x)) {
      x <- get_se_assaydata(x,
         assay_name=assay_name,
         verbose=verbose);
      nsamples <- ncol(x);
      x_names <- colnames(x);
      if (length(x) == 0) {
         stop("assays(x)[[assay_name]] did not produce a usable data matrix.");
      }
   } else {
      nsamples <- ncol(x);
      x_names <- colnames(x);
   }
   if (nsamples == 0) {
      stop("No samples are available to plot.");
   }

   # optional hook function for each panel
   if (length(panel_hook_function) > 0) {
      if (is.function(panel_hook_function)) {
         panel_hook_function <- list(panel_hook_function);
      }
      panel_hook_function <- rep(panel_hook_function,
         length.out=nsamples);
   }

   if (length(x_names) == 0) {
      if (length(sample_labels) == nsamples) {
         if ("list" %in% class(x)) {
            names(x) <- sample_labels;
         } else {
            colnames(x) <- sample_labels;
         }
         x_names <- sample_labels;
      } else {
         x_names <- jamba::makeNames(rep("V", ncol(x)));
         if ("list" %in% class(x)) {
            names(x) <- x_names;
         } else {
            colnames(x) <- x_names;
         }
      }
   }
   if (length(sample_labels) == 0) {
      sample_labels <- x_names;
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
            c(0, 0, 1.5 * (maintitle_nlines + 3), 0)));
      }
      # if highlightPoints, add margin at the bottom
      if (length(highlightPoints) > 0 && doHighlightLegend) {
         par("oma"=pmax(par("oma"),
            c(3, 0, 0, 0)));
      }
      # if useRank=TRUE add minor margin to left
      if (TRUE %in% useRank) {
         margins <- margins + c(0, 2, 0, 0);
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

   if (doPar) {
      par("mar"=margins);
   }
   if (length(whichSamples) == 0) {
      whichSamples <- seq_len(nsamples);
   } else {
      # 0.0.25.900 updated to accept character,factor and not change numeric
      if (class(whichSamples) %in% c("character", "factor")) {
         whichSamples <- match(as.character(whichSamples),
            x_names);
         if (any(is.na(whichSamples))) {
            if (all(is.na(whichSamples))) {
               stop("whichSamples supplied as character/factor does not match any names in x");
            }
            whichSamples <- whichSamples[!is.na(whichSamples)];
         }
      }
      if (is.numeric(whichSamples)) {
         if (any(whichSamples > nsamples)) {
            whichSamples <- whichSamples[whichSamples <= nsamples];
         }
         # whichSamples <- x_names[whichSamples]
      }
   }
   names(whichSamples) <- x_names[whichSamples];
   gaveMVA <- FALSE;
   if ("data.frame" %in% class(x)) {
      x <- as.matrix(x);
   }

   # Handle noise_floor,noise_floor_value, in place of deprecated arguments
   # filterFloor,filterFloorReplacement
   deprecate_warn <- character(0);
   if (length(filterFloor) > 0) {
      deprecate_warn <- c(deprecate_warn, "filterFloor");
      noise_floor <- filterFloor;
   }
   if (length(filterFloorReplacement) > 0) {
      deprecate_warn <- c(deprecate_warn, "filterFloorReplacement");
      noise_floor_value <- filterFloorReplacement;
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
         jamba::printDebug("jammaplot(): ",
            sep="",
            c("deprecated method ", "ma_method='old'",
               " will be removed in future, please use ",
               "ma_method='jammacalc'"),
            fgText="red")
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
         if (length(noise_floor) > 0) {
            x[!is.na(x) & x <= noise_floor] <- noise_floor_value;
         }

         # apply rank here
         if (TRUE %in% useRank) {
            if (is.list(x)) {
               stop("Cannot combine useRank=TRUE with list input, it requires numeric matrix input.");
            }
            ## Convert x to rank per column
            if (verbose) {
               jamba::printDebug("jammaplot(): ",
                  "Converting input matrix to per-column rank values.");
            }
            x <- matrix_to_column_rank(x, keepNA=TRUE)
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
               useMedian=useMedian,
               controlSamples=controlSamples,
               centerGroups=centerGroups,
               returnValues=FALSE,
               returnGroupedValues=TRUE);
            centerGroups <- jamba::nameVector(centerGroups, x_names);
            ## TODO: Apply rowGroupsMeans() for custom y-values per group
            if (useMedian) {
               y <- rowMedians(yM, na.rm=TRUE);
            } else {
               y <- rowMeans(yM, na.rm=TRUE);
            }
         } else {
            if (length(customFunc) > 0) {
               y <- customFunc(x[,controlSamples,drop=FALSE], na.rm=TRUE);
            } else if (useMedian) {
               y <- rowMedians(x[,controlSamples,drop=FALSE], na.rm=TRUE);
            } else {
               y <- rowMeans(x[,controlSamples,drop=FALSE], na.rm=TRUE);
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
         }
         mvaDatas <- lapply(jamba::nameVector(whichSamples), function(i){NA});
      }
   }

   if (ma_method %in% "jammacalc") {
      ## Newer method using jammacalc()
      if (verbose) {
         jamba::printDebug("jammaplot(): ",
            "Calling jammacalc().",
            "useMedian:", useMedian);
      }
      mvaDatas <- jammacalc(x=x,
         na.rm=TRUE,
         useMedian=useMedian,
         controlSamples=controlSamples,
         centerGroups=centerGroups,
         controlFloor=controlFloor,
         naControlAction=naControlAction,
         naControlFloor=naControlFloor,
         groupedX=groupedX,
         grouped_mad=groupedMAD,
         whichSamples=whichSamples,
         noise_floor=noise_floor,
         noise_floor_value=noise_floor_value,
         naValue=filterNAreplacement,
         centerFunc=centerGeneData,
         returnType="ma_list",
         mad_row_min=outlierRowMin,
         useRank=useRank,
         verbose=verbose,
         ...);
      # for useRank=TRUE adjust default ylim
      if (TRUE %in% useRank && all(abs(ylim) == 4)) {
         ylim <- c(-1,1) * round(nrow(x) * .40);
      }
   } else {
      if (length(centerGroups) > 0) {
         objectCtr <- centerGeneData(x,
            centerGroups=centerGroups,
            useMedian=useMedian,
            returnGroupedValues=FALSE,
            returnValues=TRUE,
            controlSamples=controlSamples,
            verbose=verbose,
            ...);
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
      jamba::printDebug("jammaplot(): ",
         "mvaMADs:");
      print(format(digits=2, mvaMADs));
      jamba::printDebug("jammaplot(): ",
         "mvaMADfactors:");
      print(format(digits=2, mvaMADfactors));
      jamba::printDebug("jammaplot(): ",
         "mvaMADoutliers:",
         mvaMADoutliers);
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
    highlightCex=NULL,
    highlightPch=NULL,
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
               col=unname(unlist(highlightColor)),
               pch=unname(unlist(highlightPch)),
               pt.cex=unname(unlist(highlightCex)),
               pt.bg=NULL);
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

   ######################################
   ## Iterate each plot panel
   ##
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
            highlightPch=highlightPch,
            highlightCex=highlightCex,
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
                  highlightPch=highlightPch,
                  highlightCex=highlightCex,
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
         # titleText <- x_names[i];
         titleText <- sample_labels[i];

         ## Calculate the MAD, i.e. the median absolute deviation from zero
         mvaMAD <- mvaMADs[whichSamples[i]];
         #mvaMAD <- median(abs(mvaData[,"y"]), na.rm=TRUE);

         ## Determine color ramp per panel
         if (length(outlierMAD) > 0 && x_names[i] %in% mvaMADoutliers) {
            colrampUse <- colrampOutlier[[i]];
         } else {
            colrampUse <- colramp[[i]];
         }

         if (length(mvaData) == 0 ||
               nrow(mvaData) == 0 ||
               all(is.na(mvaData[,"y"]))) {
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
               yaxt="n",
               transformation=transformation,
               #col=smoothPtCol,
               useRaster=useRaster,
               nrpoints=nrpoints,
               fillBackground=fillBackground,
               applyRangeCeiling=applyRangeCeiling,
               ...);
            y_axis_at <- pretty(ylim);
            axis(2,
               las=2,
               at=y_axis_at,
               labels=format(y_axis_at,
                  trim=TRUE,
                  big.mark=","),
               ...)
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
         ## Optionally indicate control with asterisk in the plot title
         if (x_names[i] %in% controlSamples &&
               "titlestar" %in% controlIndicator) {
            titleText <- paste0(titleText,
               "*",
               groupSuffix[i]);
         } else {
            titleText <- paste0(titleText,
               groupSuffix[i]);
         }
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
               adjPreset=titleAdjPreset,
               panelWidth="minimum",
               txt=titleText,
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
            title(main=titleText,
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
               adjPreset=subtitleAdjPreset,
               txt=subtitle[i],
               boxColor=subtitleBoxColor[i],
               boxBorderColor=jamba::makeColorDarker(subtitleBoxColor[i]),
               labelCol=jamba::setTextContrastColor(subtitleBoxColor[i]),
               labelCex=titleCex[i]*0.9,
               drawBox=TRUE,
               font=titleFont[i],
               ...);
         }
         ## Optionally print control indicator inside the plot panel
         if (x_names[i] %in% controlSamples &&
               "labelstar" %in% controlIndicator) {
            preset_star <- "topleft";
            adjPreset_star <- "bottomright";
            if (doTitleBox) {
               if (grepl("bottom", titleAdjPreset)) {
                  adjPreset_star <- gsub("bottom",
                     "top",
                     adjPreset_star)
               }
            }
            jamba::drawLabels(preset=preset_star,
               adjPreset=adjPreset_star,
               panelWidth="default",
               txt=" * ",
               drawBox=FALSE,
               labelCol="black",
               labelCex=titleCex[i]*1.5,
               font=titleFont[i],
               ...);
         }

         ## Optionally print the MAD factor in the bottom right corner
         panelCol <- head(colrampUse(51), 1);
         labelBaseCol <- jamba::setTextContrastColor(panelCol);
         labelBaseColL <- jamba::col2hcl(labelBaseCol)["L",];
         labelCol <- ifelse(
            mvaMADfactors[whichSamples[i]] >= outlierMAD,
            ifelse(labelBaseColL > 50,
               "gold",
               "red3"),
            labelBaseCol);
         if (length(outlierMAD) > 0 && displayMAD == 1) {
            jamba::drawLabels(preset="bottomright",
               labelCex=titleCex[i]*0.9,
               font=titleFont[i],
               txt=paste0("MAD x",
                  format(digits=2,
                     mvaMADfactors[whichSamples[i]])
               ),
               drawBox=FALSE,
               labelCol=labelCol
               );
         } else if (length(outlierMAD) > 0 && displayMAD == 2) {
            jamba::drawLabels(preset="bottomright",
               labelCex=titleCex[i]*0.9,
               font=titleFont[i],
               txt=paste0("MAD:",
                  format(digits=2,
                     mvaMADs[whichSamples[i]])),
               drawBox=FALSE,
               labelCol=labelCol
            );
         }
         #mvaData;
         ## Optional per-panel hook function
         if (length(panel_hook_function) > 0 &&
               length(panel_hook_function[[i]]) > 0 &&
               is.function(panel_hook_function[[i]])) {
            panel_hook_function[[i]](i,
               ...)
         }
      }
      ## End of the per-panel MVA plot loop
      check_panel_page(iPanelNumber=prod(par("mfrow")) + 1,
         maintitle,
         doHighlightLegend=doHighlightLegend,
         highlightColor=highlightColor,
         highlightPch=highlightPch,
         highlightCex=highlightCex,
         verbose=verbose,
         ...);
   }

   invisible(mvaDatas);
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
#' @family jam utility functions
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


