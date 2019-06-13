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
#' fundamentally different signal distributions.
#'
#' Another very useful technique is to center by sample group,
#' for example `centerGroups=sample_group`.
#' This technique produces MA-plots that depict the
#' "difference from group", and is very useful for visualizing
#' the variance of replicates in sample groups. Outlier samples
#' are sometimes very apparent at this step. See `displayMAD` and
#' `outlierMAD` below for detection of potential outlier samples.
#'
#' Note: filterNeg=TRUE will floor log2 values at 0, which has the
#' effect of adding the "45 degree lines" originating from c(0,0).
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
#' The argument `displayMAD` will display the per-sample MAD factor
#' relative to its `centerGroups` value, if provided. The MAD threshold
#' is calculated using values whose mean/median is at least
#' `outlierRowMin`, based upon the median absolute difference (MAD)
#' from the mean/median. The median MAD is defined as the reference MAD,
#' and all other MAD factors are multiples of that reference. For example,
#' 1xMAD is equal to the median MAD value, and 2xMAD is two times higher
#' than the median MAD value.
#'
#' Potential outlier samples should not be evaluated when supplying
#' `controlSamples` unless those samples represent more than
#' experimantal control samples.
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
#' @param x numeric matrix typically containing log-normal measurements,
#'    with measurement rows, and sample columns.
#' @param colramp one of several inputs recognized by
#'    \code{\link{getColorRamp}}. It typically recognizes either the name of
#'    a color ramp from RColorBrewer, the name of functions from the
#'    \code{\link[viridis]{viridis}} package, or single R colors, or
#'    a vector of R colors.
#' @param colrampOutlier one of several inputs recognized by
#'    \code{\link{getColorRamp}}, to be used only in panels determined to be
#'    outliers based upon MAD outlier logic. Typically the color ramp is
#'    identical to \code{colramp} except the background color indicates
#'    an error, for example yellow background instead of white.
#' @param outlierColor character R color, used when `colrampOutlier is NULL`
#'    to substitute as the first color from `colramp`. This method keeps the
#'    color ramp consistent, but changes the background color to highlight
#'    the plot panel as an outlier.
#' @param applyRangeCeiling see \code{\link{jamba::plotSmoothScatter}} for details,
#'    logical whether to apply a noise floor and ceiling to display points
#'    outside the plot region at the boundaries of the plot. This parameter
#'    is useful to set FALSE when cropping a plot, in order to remove outlier
#'    points from the density calculation. For example when there is a
#'    very large number of zero values, it can be helpful to define
#'    \code{xlim=c(0.1,20), applyRangeCeiling=FALSE}.
#' @param whichSamples NULL or integer vector, representing an index of samples
#'    to include in the MA-plots. Note that all samples are used to create
#'    the MA-plot data, but only these samples will be plotted, which is useful
#'    when wanting to view only a subset of samples in detail, but where the
#'    MA-plot data is identical to the full dataset.
#' @param maintitle the main title for the overall plot, printed at the
#'    top of each page in the top outer margin.
#' @param maintitleCex numeric cex character expansion used to resize the
#'    `maintitle` when supplied.
#' @param subtitle NULL or character text to be drawn at the bottom center
#'    of each plot panel.
#' @param subtitlePreset character value describing where to position the
#'    subtitle, using terms valid in `jamba::coordPresets()`. The default
#'    "bottomleft" places the subtitle at the bottom-left corner of each
#'    panel.
#' @param titleBoxColor,subtitleBoxColor vector of colors applied to title
#'    text, or subtitle text, respectively. When `doTitleBox=TRUE`
#'    one or no value is supplied, it defines colors using
#'    `jamba::setTextContrastColor` to use a contrasting color.
#' @param titleColor vector of colors applied to title text. When doTitleBox=TRUE
#'    one or no value is supplied, it defines colors using
#'    \code{\link{jamba::setTextContrastColor}} to use a contrasting color.
#' @param doTitleBox logical whether to draw plot titles using a colored box.
#' @param titleFont integer font compatible with \code{par("font")}. Values
#'    are recycled across panels.
#' @param titlePreset character value describing where to position the
#'    subtitle, using terms valid in `jamba::coordPresets()`. The default is
#'    "top" which centers the label at the top of each panel.
#' @param xlab,ylab character x- and y-axis labels.
#' @param xlabline,ylabline numeric number of text lines as used by
#'    \code{title} to position labels relative to the plot area.
#' @param groupSuffix character text appended to each plot panel title,
#'    useful for indicating how a sample was processed, beyond just using
#'    colnames. By default, `groupSuffix` is blank, except when
#'    supplying `centerGroups`, in which case the `centerGroups` value is
#'    used as a suffix.
#' @param highlightPoints NULL, character vector, or list of character vectors,
#'    containing points to highlight based upon \code{rownames(object)}. If
#'    a list is supplied, then colors from \code{highlightColors} are recycled
#'    to \code{length(highlightPoints)} so that each list element can receive
#'    its own color. This coloring may be useful for coloring several subsets
#'    of rows, for example housekeeper genes, ribosomal genes, or platform
#'    positive- and negative-control genes.
#' @param highlightColor,highlightCex values recycled to the length of
#'    highlightPoints. If highlightPoints is a list, then a list or vector
#'    can be supplied.
#' @param doHighlightPolygon logical whether to draw a polygon encompassing
#'    highlighted points. If highlightPoints is a list, a polygon is drawn
#'    around each list element subset of points, using the appropriate
#'    highlightColors for each list element. The polygon is defined by
#'    `grDevices::chull()`.
#' @param highlightPolygonAlpha numeric value between 0 and 1 indicating the
#'    alpha transparency to use for the polygon color, based upon
#'    highlightColor.
#' @param smoothPtCol color used when nrpoints>0, and
#'    \code{\link{jamba::plotSmoothScatter}} draws this many points in the extremities.
#' @param margins vector of margins compatible with \code{par("mar")}. Defaults
#'    are applied, but provided here for convenient override.
#' @param useRaster logical whether to define \code{useRaster=TRUE} in
#'    \code{\link{jamba::plotSmoothScatter}}, which ultimately calls
#'    \code{\link{imageDefault}}. When TRUE, it creates a much smaller
#'    plot object, because it returns a raster image instead of a set of
#'    rectangles when rendering the image.
#' @param ncol,nrow integer number of columns and rows used by
#'    \code{par("mfrow")} when \code{doPar=TRUE} is defined. These values
#'    are helpful in defining a fixed layout for multiple MA-plot panels.
#' @param doPar logical whether \code{par("mfrow")} should be defined, in
#'    order to fit multiple panels in one plot window. Set to FALSE if the plot
#'    layout is already defined, or if plotting one panel per page, for
#'    example.
#' @param las integer value 1 or 2, indicating whether axis labels should be
#'    parallel or perpendicular to the axes, respectively. This default value
#'    is defined and provided here for convenience to override it.
#' @param ylim,xlim NULL or vector length 2 indicating the y- and x-axis
#'    ranges, respectively. These values are useful to define upfront when
#'    it is best to focus on fixed, consistent ranges across all samples. The
#'    default \code{ylim=c(-4,4)} represents 16-fold range in normal space, and
#'    is typically a reasonable starting point for most purposes. If values are
#'    all between -1.5 and 1.5, then it is still good to keep that range in
#'    context of -4 and 4, to indicate that the observed noise is lower than
#'    typically observed. It is also helpful to use xlim to trim off zeros
#'    for some data where there might be many undetected entries, combined
#'    with \code{applyRangeCeiling=FALSE} which also removed outlier points
#'    from the plot density calculation. The result can notably boost the
#'    detail displayed in the plot range.
#' @param controlSamples vector of colnames(object) which should be used
#'    as control samples in the data centering step. Typical MA-plots are
#'    calculated relative to the overall mean, however it can be insightful
#'    to calculate values relative to a known control group.
#' @param centerGoups vector of groups, of \code{length(colnames(object))}
#'    indicating optional subgroups to use when performing data centering.
#'    For example it may be appropriate to center samples of one cell type
#'    only relative to other samples of the same cell type. See
#'    \code{\link{centerGeneData}} for more specifics.
#' @param useMean logical whether to perform centering using mean or
#'    median. The default median is intended to represent the data without
#'    being skewed by outliers, however it can be informative to use the mean
#'    values when there are few replicates, or when there are no particular
#'    sample outliers.
#' @param customFunc an optional function that is used instead of mean or
#'    median. It should take a matrix input, and return a numeric vector
#'    output summarizing each row in \code{object}. It is intended to
#'    provide custom row statistics, for example geometric mean, or other
#'    row summary function.
#' @param filterNA,filterNeg logical indicating whether to remove NA or
#'    negative values before creating the MA-plot. It can be useful to
#'    filter negative values, which are often defined as noise in upstream
#'    normalization steps, so they do not impart visible features into the
#'    MA-plot since they are based upon noise. Negative values are set to
#'    zero when filtered, which gives a characteristic 45 degree angle
#'    symmetry above and below zero.
#' @param filterFloor,filterFloorReplacement numeric or NULL, when numeric
#'    values in \code{object} are below \code{filterFloor}, they are replaced
#'    with \code{filterFloorReplacement}. This mechanism is an alternative to
#'    \code{filterNeg}, in that it allows defining a filter above or below
#'    zero. For some platform data , it might be useful to define a
#'    \code{filterFloor} roughly equivalent to its noise threshold.
#' @param transFactor numeric adjustment used by
#'    \code{\link{jamba::plotSmoothScatter}}. The default is defined, but available
#'    here as a convenient override. The value is used in a power function in
#'    the form \code{transformation=function(x)x^transFactor}. Lower values
#'    make the scatterplot point density more intense.
#' @param nrpoints integer or NULL, the number of points to display on the
#'    extremity, as used by \code{\link{jamba::plotSmoothScatter}}.
#' @param smoothScatterFunc function used to plot smooth scatter plot, by
#'    default \code{\link{jamba::plotSmoothScatter}}.
#' @param ablineH,ablineV numeric vector indicating where to draw indicator
#'    lines, at these horizontal or vertical positions, respectively.
#' @param doTxtplot logical whether to plot results in colored text output,
#'    currently under development and is still experimental.
#' @param blankPlotPos NULL or integer vector indicating which plot panels
#'    should be used as blank filler positions. These blank panels can be
#'    used to help organize multiple plot panels to indicate sample groups.
#'    A mechanism similar to ggplot2 facets, except with the ability to
#'    insert known blank positions for organization. The facet_grid mechanic
#'    may be used in future, however is not always ideal depending upon the
#'    experiment design, and so for now a simple manual method is provided.
#' @param displayMAD logical whether to display MAD (median absolute deviation)
#'    values for each panel. MAD values are defined for each centerGroup if
#'    supplied, or globally if not. The MAD value is CALCULATED by taking the
#'    median deviation from zero as displayed on the y-axis of each MA-plot
#'    panel, then taking the median of those values. Then a MAD factor is
#'    calculated for each MA-plot panel by taking its median devation divided
#'    by the MAD value for its centerGroup. Most samples should therefore be
#'    close to 1. A value of 2 is interpreted as "two times the MAD factor for
#'    its center group." An outlierMAD of 5x is a reasonable cutoff, which
#'    roughly says that sample has 5 times higher median deviation than other
#'    samples in its center group. If all samples are noisy, it would define
#'    a high MAD value, and therefore all samples would be expected to have
#'    similar MAD factors, near 1.
#' @param groupedMAD logical indicating whether the MAD calculation should
#'    operate only within `centerGroups` (`groupedMAD=TRUE`) or whether the
#'    MAD calculation is shared across all samples (`groupedMAD=FALSE`).
#' @param outlierMAD the MAD factor threshold to define a sample an outlier.
#' @param outlierRowMin the minimum value as displayed on the x-axis, for
#'    a row in \code{object} to be used in the MAD outlier calculations. This
#'    threshold is intended to define a noise threshold, so that outlier
#'    samples are determined by using rows above the noise threshold. Thus
#'    a sample outlier would need to have greater than five times median
#'    absolute deviation compared to its center group, based only upon
#'    measurements above noise.
#' @param fillBackground logical sent to \code{\link{jamba::plotSmoothScatter}}
#'    indicating whether to fill the plot panel with the background color,
#'    which is typically the first value in the \code{colramp}. Set to TRUE
#'    when axis ranges may be defined which are outside the range of values
#'    for any particular column in \code{object}.
#' @param ma_method character string indicating whether to perform
#'    MA-plot calculations using the old method `"old"`; or `"jammacalc"`
#'    which uses the newer function `jammacalc()`.
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
 colramp=c("white", "lightblue", "blue", "navy", "orange", "orangered2"),
 colrampOutlier=NULL,
 outlierColor="palegoldenrod",
 whichSamples=NULL,
 maintitle=NULL,
 maintitleCex=1.8,
 subtitle=NULL,
 subtitlePreset="bottomleft",
 titleCexFactor=1,
 titleCex=NULL,
 doTitleBox=TRUE,
 titleBoxColor="#DDBB9977",
 subtitleBoxColor=titleBoxColor,
 titleColor="black",
 titleFont=2,
 titlePreset="top",
 xlab="",
 xlabline=2,
 ylab="",
 ylabline=1.5,
 groupSuffix=NULL,
 highlightPoints=NULL,
 highlightPch=21,
 highlightCex=1,
 highlightColor="#00AAAA66",
 doHighlightPolygon=FALSE,
 highlightPolygonAlpha=0.3,
 smoothPtCol="#00000055",
 margins=c(3.5, 3, 0.7, 0.5),
 useRaster=TRUE,
 ncol=NULL,
 nrow=NULL,
 doPar=TRUE,
 las=2,
 ylim=c(-4,4),
 xlim=NULL,
 controlSamples=colnames(x),
 centerGroups=NULL,
 groupedX=TRUE,
 useMean=FALSE,
 customFunc=NULL,
 filterNA=TRUE,
 filterNAreplacement=NA,
 filterNeg=TRUE,
 filterFloor=0,
 filterFloorReplacement=filterFloor,
 transFactor=0.18,
 nrpoints=50,
 smoothScatterFunc=jamba::plotSmoothScatter,
 applyRangeCeiling=TRUE,
 doTxtplot=FALSE,
 ablineV=0, ablineH=c(-2,0,2),
 blankPlotPos=NULL,
 outlierMAD=5,
 outlierRowMin=5,
 groupedMAD=TRUE,
 displayMAD=FALSE,
 fillBackground=TRUE,
 ma_method=c("old", "jammacalc"),
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
   if (!suppressPackageStartupMessages(require(jamba))) {
      stop("jammaplot() requires the jamba package.");
   }
   if (suppressPackageStartupMessages(require(matrixStats))) {
      doMS <- TRUE;
   } else {
      doMS <- FALSE;
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

   ## If colramp is "character" we use getColorRamp() to work it out
   if ("character" %in% class(colramp)) {
      colramp <- colorRampPalette(getColorRamp(colramp));
   }
   firstColor <- head(colramp(10), 1);
   if (length(colrampOutlier) > 0) {
      if (jamba::igrepHas("character", class(colrampOutlier))) {
         colrampOutlier <- colorRampPalette(jamba::getColorRamp(colrampOutlier));
      }
   }
   if (length(colrampOutlier) == 0 ||
         !is.function(colrampOutlier)) {
      firstColorL <- jamba::col2hcl(firstColor)["L",];
      if (length(outlierColor) < 1) {
         if (jamba::col2hcl(firstColor)["L",] < 70) {
            outlierColor <- "#551100";
         } else {
            outlierColor <- "palegoldenrod";
         }
      }
      colrampOutlier <- colorRampPalette(
         c(outlierColor,
            tail(colramp(101), 100)
      ));
   }

   ## Define groupSuffix using centerGroups, unless otherwise specified
   if (length(groupSuffix) > 0) {
      groupSuffix <- rep(groupSuffix,
         length.out=ncol(x));
   } else if (length(centerGroups) > 0) {
      groupSuffix <- paste0(" vs ", centerGroups);
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
   }
   if (length(maintitle) > 0) {
      maintitle_nlines <- length(unlist(strsplit(maintitle, "\n")));
      par("oma"=pmax(par("oma"),
         c(0, 0, 1.5+1.5*maintitle_nlines, 0)));
   }
   if (length(titleCex) == 0) {
      titleCex <- 1;
      #if (length(maintitle) == 0) {
      #   titleCex <- titleCexFactor;
      #} else {
      #   titleCex <- titleCexFactor + 1/(log2(10 + nchar(maintitle)));
      #}
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
            centerGroups <- nameVector(centerGroups, x_names);
            ## TODO: Apply rowGroupsMeans() for custom y-values per group
            if (useMean) {
               y <- rowMeans(yM, na.rm=TRUE);
            } else {
               if (doMS) {
                  y <- rowMedians(yM, na.rm=TRUE);
               } else {
                  y <- apply(yM, 1, function(i){
                     median(i, na.rm=TRUE);
                  });
               }
            }
         } else {
            if (length(customFunc) > 0) {
               y <- customFunc(x[,controlSamples,drop=FALSE], na.rm=TRUE);
            } else if (useMean) {
               y <- rowMeans(x[,controlSamples,drop=FALSE], na.rm=TRUE);
            } else {
               if (doMS) {
                  y <- rowMedians(x[,controlSamples,drop=FALSE], na.rm=TRUE);
               } else {
                  y <- apply(x[,controlSamples,drop=FALSE], 1, function(i){
                     median(i, na.rm=TRUE);
                  });
               }
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
         mvaDatas <- lapply(nameVector(whichSamples), function(i){NA});
      }
   }

   if (ma_method %in% "jammacalc") {
      ## Newer method using jammacalc()
      if (verbose) {
         printDebug("jammaplot(): ",
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
         verbose=verbose,
         ...);
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
      mvaMADs <- sapply(nameVector(whichSamples, x_names[whichSamples]), function(i){
         if (verbose) {
            jamba::printDebug("   i:",
               i);
         }
         #mvaData <- mvaDatas[[i]];
         mvaData <- mvaDatas[[x_names[i]]];
         iWhich <- (!is.na(mvaData[,"x"]) & abs(mvaData[,"x"]) >= outlierRowMin);
         rmNA(
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
   check_panel_page <- function(iPanelNumber, maintitle) {
      if (length(maintitle) == 0) {
         return();
      }
      max_panels <- prod(par("mfrow"));
      if (iPanelNumber > 1 && ((iPanelNumber-1) %% max_panels) == 0) {
         # print maintitle
         if (length(maintitle) > 0) {
            maintitle_txt <- paste(unlist(strsplit(maintitle, "\n")),
               collapse="\n");
            title(outer=TRUE,
               main=maintitle_txt,
               cex.main=maintitleCex);
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
         check_panel_page(iPanelNumber, maintitle);
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
               check_panel_page(iPanelNumber, maintitle);
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

         colrampUse <- colramp;
         if (length(outlierMAD) > 0 && x_names[i] %in% mvaMADoutliers) {
            colrampUse <- colrampOutlier;
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
               colramp=colrampUse,
               xlab="",
               ylab="",
               las=las,
               xlim=xlim,
               ylim=ylim,
               transformation=transformation,
               col=smoothPtCol,
               useRaster=useRaster,
               nrpoints=nrpoints,
               fillBackground=fillBackground,
               applyRangeCeiling=applyRangeCeiling,
               ...);
         }
         ## Add axis labels
         title(xlab=xlab, line=xlabline);
         title(ylab=ylab, line=ylabline);

         ## Add some axis lines across the plot for easy visual reference
         if (!is.null(ablineH)) {
            abline(h=ablineH, col="#44444488", lty="dashed", lwd=1, ...);
         }
         if (!is.null(ablineV) && any(xlim[1] < ablineV & xlim[2] > ablineV)) {
            ablineVuse <- ablineV[(xlim[1] < ablineV & xlim[2] > ablineV)];
            abline(col="grey30", v=ablineVuse, lty="dashed", ...);
         }

         ## Optionally highlight some subset of points
         if (!is.null(highlightPoints)) {
            if (!class(highlightPoints) %in% c("list")) {
               highlightPoints <- list(highlightPoints=highlightPoints);
            }
            if (!class(highlightColor) %in% c("list")) {
               highlightColor <- as.list(highlightColor);
            }
            highlightColor <- rep(highlightColor,
               length.out=length(highlightPoints));
            if (!class(highlightPch) %in% c("list")) {
               highlightPch <- as.list(highlightPch);
            }
            highlightPch <- rep(highlightPch,
               length.out=length(highlightPoints));
            if (!class(highlightCex) %in% c("list")) {
               highlightCex <- as.list(highlightCex);
            }
            highlightCex <- rep(highlightCex,
               length.out=length(highlightPoints));

            hp1 <- lapply(seq_along(highlightPoints), function(highI){
               highP <- highlightPoints[[highI]];
               hiData <- mvaData[which(rownames(mvaData) %in% highP),,drop=FALSE];

               ## Make sure to restrict the y-values to fit within the axis limits,
               ## consistent with smoothScatterFunc().
               yValues <- noiseFloor(hiData[,"y"], minimum=min(ylim),
                  ceiling=max(ylim), ...);

               ## Optionally draw a polygon hull around highlighted points.
               if (doHighlightPolygon && length(yValues) > 2) {
                  hiM <- cbind(x=hiData[,"x"],
                     y=yValues);
                  hiHull <- hiM[grDevices::chull(hiM),,drop=FALSE];
                  #hiHull <- points2polygonHull(data.frame(x=hiData[,"x"],
                  #   y=yValues),
                  #   returnClass="matrix");
                  if (is.null(highlightPolygonAlpha)) {
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
               #boxCexAdjust=c(1.1,1.3),
               font=titleFont[i],
               ...);

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
                  #boxCexAdjust=c(1.1,1.3),
                  font=titleFont[i],
                  ...);
            }
            par("xpd"=parXpd);
         } else {
            titleLine <- margins[3] - 1.5;
            title(main="",#maintitle,
               sub=subtitle,
               cex.sub=titleCex[i],
               cex.main=titleCex[i],
               line=titleLine,
               col.main=titleBoxTextColor,
               font.main=titleFont[i], ...);
            title(main=paste(titleText, groupSuffix[i]),
               cex.main=titleCex[i],
               line=titleLine-1,
               col.main=titleBoxTextColor,
               font.main=titleFont[i],
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
      ## Now check to see if we should pad blank panels at the end of the sequence
      if (length(blankPlotPos) > 0) {
         iPanelNumber <- iPanelNumber + 1;
         if (verbose) {
            jamba::printDebug("iPanelNumber: ",
               iPanelNumber);
         }
         while(iPanelNumber %in% blankPlotPos) {
            jamba::nullPlot(doBoxes=FALSE);
            if (verbose) {
               jamba::printDebug("Inserted a blank panel at position: ",
                  iPanelNumber);
            }
            iPanelNumber <- iPanelNumber + 1;
         }
      }
   }
   if (length(maintitle) > 0) {
      maintitle_txt <- paste(unlist(strsplit(maintitle, "\n")),
         collapse="\n");
      title(outer=TRUE,
         main=maintitle_txt,
         cex.main=maintitleCex);
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
#'    \code{\link[matrixStats]{rowMedians}} for calculations, otherwise
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
   useMS <- suppressPackageStartupMessages(require(matrixStats));

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
      rmNULL(
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
            if (useMS) {
               indataSds <- rowSds(centeredData[,controls,drop=FALSE], na.rm=TRUE);
            } else {
               indataSds <- apply(indata[,controls, drop=FALSE], 1, function(x){
                  sd(x, na.rm=TRUE);
               });
            }
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
      if (useMS) {
         indataMedians <- rowMedians(indata[,controls,drop=FALSE], na.rm=TRUE);
      } else {
         indataMedians <- apply(indata[,controls,drop=FALSE], 1, function(x){
            median(x, na.rm=TRUE);
         });
      }
      if (returnGroupedValues && !returnValues) {
         centeredData <- matrix(indataMedians,
            ncol=1,
            dimnames=list(rownames(indata), groupPrefix));
      } else {
         centeredData <- sweep(indata, 1, indataMedians);
         if (scale %in% "row") {
            if (useMS) {
               indataMads <- rowMads(centeredData[,controls, drop=FALSE], na.rm=TRUE);
            } else {
               indataMads <- apply(indata[,controls, drop=FALSE], 1, function(x){
                  mad(x, na.rm=TRUE);
               });
            }
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
         printDebug("centerGeneData_new(): ",
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
         rmNA(centerGroups,
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
   x_group <- rowGroupMeans(x[,controls_v, drop=FALSE],
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
      printDebug("jammacalc(): ",
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
   if (length(noise_floor) == 1 && !is.infinite(noise_floor)) {
      if (all(noise_floor == noise_floor_value)) {
         if (rmNA(naValue=0, any(x < noise_floor))) {
            if (verbose) {
               printDebug("jammacalc(): ",
                  "flooring values below ",
                  noise_floor,
                  " to ",
                  noise_floor_value);
            }
            x[x < noise_floor] <- noise_floor_value;
         }
      } else if (rmNA(naValue=0, any(x <= noise_floor))) {
         if (verbose) {
            printDebug("jammacalc(): ",
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
            printDebug("jammacalc(): ",
               "replacing NA with ",
               naValue);
         }
      }
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
         printDebug("jammacalc(): ",
            "Using group values for MA-plot x-axis values");
      }
   } else {
      groupedX <- FALSE;
      if (useMedian) {
         if (verbose) {
            printDebug("jammacalc(): ",
               "using rowMedian() of all values for x.");
         }
         x_use <- rowMedians(x,
            na.rm=na.rm);
      } else {
         if (verbose) {
            printDebug("jammacalc(): ",
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
   x_mads <- nameVector(colMedians(abs(x_mad_m), na.rm=TRUE),
      colnames(x_ctr));
   if (grouped_mad && length(unique(centerGroups)) > 1) {
      if (verbose) {
         printDebug("jammacalc(): ",
            "Calculating grouped MAD values.");
      }
      x_grp_mads <- tapply(x_mads, centerGroups, median, na.rm=TRUE);
      x_mad_factors <- x_mads / x_grp_mads[as.character(centerGroups)];
      names(x_mad_factors) <- names(x_mads);
   } else {
      if (verbose) {
         printDebug("jammacalc(): ",
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
