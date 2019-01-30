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
#' \code{jammaplot} takes a numeric matrix, typically of gene expression data,
#' and produces an MA-plot (Bland-Altman plot).
#'
#' Purpose is to provide an MA plot, based upon the MVAplotMED function from
#' the now-defunct Agi4x44PreProcess package,  except in this case using
#' smoothScatter since 20k to 40k datapoints is too many to show without
#' them being occluded.
#'
#' You must manually specify the number columns and rows, and set doPar=TRUE
#' for this function to define multi-panel plotting, e.g.
#' \code{"nrow=5, ncol=4, doPar=TRUE"}
#' Otherwise, you can manually set these parameters with
#' \code{"par('mfrow'=c(5,4));"}
#'
#' This function also uses "useRaster=TRUE" by default which causes the
#' \code{\link{imageDefault}} function used by \code{\link{smoothScatter}}
#' to produce a temporary raster image which is then resized with
#' interpolation to produce a properly-blended representation of the full
#' image. This process also substantially reduces the image size when saving
#' as a vector output format like PDF or SVG.
#'
#' To highlight certain points, e.g. marked outliers, or statistical hits,
#' supply rownames into highlightPoints as a vector, or a list of vectors.
#' If using a list of vectors, also supply a list of colors in
#' highlightColor, so each subset is colored uniquely.
#'
#' ctrlSamples is deprecated in favor of controlSamples.
#'
#' controlSamples indicates which samples should be used for centering
#' data, by default all samples. However, centering by control samples
#' visually indicates changes from those control samples.
#'
#' centerGroups is used to apply the centering to subsets of columns
#' of data, for example if two very different cell types should be separately
#' centered for closer scrutiny of the within-cell variability.
#' One can take this further, and center each sample group, which has
#' the helpful effect of displaying variability of replicates per
#' sample group. This method is quite useful for reviewing whether any
#' replicates are notably different than counterparts in the same
#' sample group.
#'
#' Note: filterNeg=TRUE will floor log2 values at 0, which has the
#' effect of adding the "45 degree lines" originating from c(0,0).
#'
#' blankPlotPos is intended for specific plot layouts, to keep the plots
#' organized by some fixed criteria, e.g. certain sample groups per
#' row or column, but where there may not be consistent replicates per
#' group.
#' For ncol and nrow, blankPlotPos indicates which panels should be
#' left empty, where each panel is numbered by row.  For example,
#' if ncol=4, then the panels are numbered 1-4, 5-8, 9-12, etc.
#'
#' outlierMAD is only viable when centerEachGroup=TRUE; in which case
#' outlierMAD is a threshold for the MAD factor of outliers.  Each MA
#' plot panel is given a MAD score, using probes whose x-values are
#' greater than outlierRowMin;, the median of the absolute value is
#' used as the MAD.
#' These MAD values are used to calculate the global MAD, and individual
#' MADs are divided by global MAD to determine a MAD factor.
#' The outlierMAD is the threshold above which a MAD factor is considered
#'  an outlier.
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
#' @param object numeric matrix typically containing log-normal measurements,
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
#' @param maintitle the main title for each panel, either printed inside each
#'    title box when doTitleBox=TRUE, or atop each panel when doTitleBox=FALSE.
#' @param maintitleSep character separator used to combine multiple values
#'    for the maintitle, used to keep all information either on one line or
#'    split on multiple lines.
#' @param subtitle NULL or character text to be drawn at the bottom center
#'    of each plot panel.
#' @param titleBoxColor vector of colors applied to title text. When doTitleBox=TRUE
#'    one or no value is supplied, it defines colors using
#'    \code{\link{jamba::setTextContrastColor}} to use a contrasting color.
#' @param titleColor vector of colors applied to title text. When doTitleBox=TRUE
#'    one or no value is supplied, it defines colors using
#'    \code{\link{jamba::setTextContrastColor}} to use a contrasting color.
#' @param doTitleBox logical whether to draw plot titles using a colored box.
#' @param titleFont integer font compatible with \code{par("font")}. Values
#'    are recycled across panels.
#' @param xlab,ylab character x- and y-axis labels.
#' @param xlabline,ylabline numeric number of text lines as used by
#'    \code{title} to position labels relative to the plot area.
#' @param groupSuffix character text appended to each plot panel title,
#'    useful for indicating how a sample was processed, beyond just using
#'    colnames. For example, \code{groupSuffix=" vs. group median"} could
#'    indicate that a sample is centered relative to its group median, rather
#'    than the global median.
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
#'    \code{\link{points2polygonHull}}.
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
#' \dontrun{
#' jammaplot(x)
#' }
#'
#' @export
jammaplot <- function
(object,
 colramp=c("white", "lightblue", "blue", "navy", "orange", "orangered2"),
 colrampOutlier=NULL,
 outlierColor="palegoldenrod",
 whichSamples=NULL,
 maintitle=NULL,
 subtitle=NULL,
 maintitleSep="\n",
 titleCexFactor=1,
 titleCex=NULL,
 doTitleBox=TRUE,
 titleBoxColor="#DDBB9977",
 titleColor="black",
 titleFont=2,
 xlab="A", xlabline=2,
 ylab="M", ylabline=1.5,
 groupSuffix=" vs. Median",
 highlightPoints=NULL,
 highlightPch=21,
 highlightCex=1,
 highlightColor="#00AAAA66",
 doHighlightPolygon=FALSE,
 highlightPolygonAlpha=0.3,
 smoothPtCol="#00000055",
 margins=c(3.5, 3, 0.7, 0.5),
 useRaster=TRUE,
 ncol=NULL, nrow=NULL, doPar=TRUE,
 las=2,
 ylim=c(-4,4),
 xlim=NULL,
 controlSamples=colnames(object),
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
   if (doTxtplot) {
      blankPlotPos <- NULL;
      smoothScatterFunc <- function(...){
         plotTextSmoothScatter(height=20, width=80, doLegend=FALSE, ...);
      }
   }
   transformation <- function(x){
      x^transFactor;
   }
   if (is.null(colnames(object))) {
      colnames(object) <- makeNames(rep("V", ncol(object)));
   }

   ## If colramp is "character" we use getColorRamp() to work it out
   if ("character" %in% class(colramp)) {
      colramp <- colorRampPalette(getColorRamp(colramp));
   }
   firstColor <- head(colramp(10), 1);
   if (length(colrampOutlier) > 0) {
      if (jamba::igrepHas("character", class(colrampOutlier))) {
         colrampOutlier <- colorRampPalette(getColorRamp(colrampOutlier));
      }
   }
   if (length(colrampOutlier) == 0 ||
         !is.function(colrampOutlier)) {
      firstColorL <- col2hcl(firstColor)["L",];
      if (length(outlierColor) < 1) {
         if (col2hcl(firstColor)["L",] < 70) {
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

   groupSuffix <- rep(groupSuffix, length.out=ncol(object));

   nsamples <- ifelse(is.null(ncol(object)), length(object), ncol(object));
   if (doPar) {
      oPar <- par(no.readonly=TRUE);
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
   if (is.null(titleCex)) {
      if (is.null(maintitle)) {
         titleCex <- titleCexFactor;
      } else {
         titleCex <- titleCexFactor + 1/(log2(10 + nchar(maintitle)));
      }
   }
   titleCex <- rep(titleCex, length.out=nsamples);
   if (verbose) {
      jamba::printDebug("jammaplot(): ",
         "titleCex:", head(titleCex, 10), ", length=", length(titleCex));
   }
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

   titleFont <- rep(titleFont, length.out=nsamples);
   par("mar"=margins);
   nARR <- nsamples;
   if (is.null(whichSamples)) {
      whichSamples <- 1:nARR;
   } else {
      if (is.numeric(whichSamples)) {
         whichSamples <- colnames(object)[whichSamples]
      }
      whichSamples <- match(whichSamples, colnames(object));
   }
   names(whichSamples) <- colnames(object)[whichSamples];
   gaveMVA <- FALSE;
   if (class(object) %in% "data.frame") {
      object <- as.matrix(object);
   }
   if (class(object) %in% "list" &&
      all(c("x", "y") %in% colnames(object[[1]]))) {
      ## We are given MVA data from previous run, or from another source.
      gaveMVA <- TRUE;
      if (is.null(xlim)) {
         xlim <- range(sapply(object, function(iObj){
            range(c(range(iObj[,"x"])));
         }));
         #xlim <- c(0, xlim1);
      }
      mvaDatas <- object;
   } else {
      if (verbose) {
         jamba::printDebug("jammaplot(): ",
            "Processing matrix input type, dim(object)", dim(object));
      }
      if (filterNA) {
         object[is.na(object)] <- filterNAreplacement;
      }
      if (filterNeg) {
         object[object < 0] <- 0;
      }
      if (!is.null(filterFloor)) {
         object[!is.na(object) & object < filterFloor] <- filterFloorReplacement;
      }
      if (is.null(controlSamples)) {
         controlSamples <- colnames(object);
      } else {
         controlSamples <- intersect(controlSamples, colnames(object));
         if (length(controlSamples) == 0) {
            controlSamples <- colnames(object);
         }
      }
      ## Calculate the summary value per row of input data object,
      ## used as the x-axis value for each panel
      if (groupedX && length(unique(centerGroups)) > 1) {
         yM <- centerGeneData(object,
            controlSamples=controlSamples,
            centerGroups=centerGroups,
            returnValues=FALSE,
            returnGroupedValues=TRUE);
         centerGroups <- nameVector(centerGroups, colnames(object));
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
         if (!is.null(customFunc)) {
            y <- customFunc(object[,controlSamples,drop=FALSE], na.rm=TRUE);
         } else if (useMean) {
            y <- rowMeans(object[,controlSamples,drop=FALSE], na.rm=TRUE);
         } else {
            if (doMS) {
               y <- rowMedians(object[,controlSamples,drop=FALSE], na.rm=TRUE);
            } else {
               y <- apply(object[,controlSamples,drop=FALSE], 1, function(i){
                  median(i, na.rm=TRUE);
               });
            }
         }
      }
      if (verbose) {
         jamba::printDebug("jammaplot(): ",
            "Processing matrix input type, dim(object)", dim(object));
         jamba::printDebug("jammaplot(): ",
            "Processing matrix input type, dim(y):", dim(y));
      }
      if (is.null(xlim)) {
         ## Calculate x range using the range of data overall
         xlim <- range(c(range(object, na.rm=TRUE), range(y, na.rm=TRUE)));
         #xlimMax <- mean(c(max(object, na.rm=TRUE), max(y, na.rm=TRUE)));
         #xlim <- c(filterFloor, xlimMax);
      }
      mvaDatas <- lapply(nameVector(whichSamples), function(i){NA});
   }

   if (!is.null(centerGroups)) {
      objectCtr <- centerGeneData(indata=object,
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
         jamba::printDebug("head(y): ", round(head(y), digits=2),
            fgText=c("orange", "lightblue"));
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
            iCol <- colnames(object)[i];
            iGroup <- centerGroups[iCol];
            yUse <- yM[,iGroup];
            #x <- objectCtr[,i] + y;
            x <- objectCtr[,i] + yUse;
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
            x <- object[,i];
            yUse <- y;
         }
         mvaValues <- c(yUse, (x - yUse));
         mvaData <- matrix(data=mvaValues,
            nrow=length(x),
            byrow=FALSE,
            dimnames=list(genes=rownames(object),
               MAcoords=c("x", "y")));
         mvaDatas[[i]] <- mvaData;
         names(mvaDatas)[i] <- colnames(object)[i];
      }
   }

   ## Optionally calculate the MAD per panel
   ## assumptions: that zero should be the median, all deviations are absolute values since
   ## we are comparing with zero, and we are not subtracting an actual median value
   if (length(outlierMAD) > 0) {
      if (verbose) {
         jamba::printDebug("Calculating MVA data MAD factors.");
      }
      mvaMADs <- sapply(nameVector(whichSamples, colnames(object)[whichSamples]),
         function(i){
         if (verbose) {
            jamba::printDebug("   i:", i, c("orange", "lightblue"));
         }
         mvaData <- mvaDatas[[i]];
         median(abs(mvaData[abs(mvaData[,"x"]) >= outlierRowMin,"y",drop=FALSE]),
            na.rm=TRUE);
      });

      ## Calculate MAD factor, a ratio to the median
      if (groupedMAD && length(centerGroups) > 0) {
         ## if grouped, use the median of each center group
         mvaMADfactors <- unlist(unname(
            tapply(mvaMADs, centerGroups, function(i){
               i / median(i);
            })))[names(mvaMADs)];
      } else {
         ## if not grouped, use the overall median across all samples
         mvaMAD <- median(mvaMADs);
         mvaMADfactors <- mvaMADs/mvaMAD;
      }

      #mvaMADthreshold <- outlierMAD * mvaMAD;
      names(mvaMADfactors) <- whichSamples;
      names(mvaMADs) <- whichSamples;
      mvaMADoutliers <- colnames(object)[whichSamples][which(mvaMADfactors >= outlierMAD)];
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
      attr(mvaDatas, "outlierMAD") <- outlierMAD;
      attr(mvaDatas, "outlierRowMin") <- outlierRowMin;
   } else {
      mvaMADoutliers <- NULL;
   }

   ## iPanelNumber keeps track of the numbered panels as they are plotted,
   ## so we can insert blank panels at the specified positions.
   if (doPlot) {
      iPanelNumber <- 0;
      for(i in whichSamples) {
         mvaData <- mvaDatas[[i]];
         iPanelNumber <- iPanelNumber + 1;
         if (verbose) {
            jamba::printDebug("iPanelNumber: ", iPanelNumber,
               fgText=c("orange", "lightgreen"));
         }
         if (length(blankPlotPos) > 0) {
            while(iPanelNumber %in% blankPlotPos) {
               nullPlot(doBoxes=FALSE);
               if (verbose) {
                  jamba::printDebug("   Inserted a blank panel at position: ",
                     iPanelNumber, fgText=c("orange", "lightblue"));
               }
               iPanelNumber <- iPanelNumber + 1;
               if (verbose) {
                  jamba::printDebug("new iPanelNumber: ", iPanelNumber,
                     fgText=c("orange", "lightgreen"));
               }
            }
         }
         par("mar"=margins);
         if (!is.null(maintitle)) {
            what <- paste(c(maintitle, colnames(object)[i], groupSuffix[i]),
               collapse=maintitleSep);
            groupName <- paste(c(colnames(object)[i], groupSuffix[i]),
               collapse="");
         } else {
            what <- paste(c(colnames(object)[i], groupSuffix[i]),
               collapse=maintitleSep);
            groupName <- paste(c(colnames(object)[i], groupSuffix[i]),
               collapse="");
         }
         titleText <- names(mvaDatas)[i];

         ## Calculate the MAD, i.e. the median absolute deviation from zero
         mvaMAD <- median(abs(mvaData[,"y"]), na.rm=TRUE);

         colrampUse <- colramp;
         if (!is.null(outlierMAD) && colnames(object)[i] %in% mvaMADoutliers) {
            colrampUse <- colrampOutlier;
         }
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
                  hiHull <- points2polygonHull(data.frame(x=hiData[,"x"],
                     y=yValues),
                     returnClass="matrix");
                  if (is.null(highlightPolygonAlpha)) {
                     highlightPolyAlpha <- jamba::col2alpha(highlightColor[[highI]]) / 2;
                  }
                  highlightPolyColor <- jamba::alpha2col(highlightColor[[highI]],
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
            if (1 == 2) {
               legend("topright",
                  inset=0.02,
                  cex=titleCex[i],
                  text.font=titleFont[i],
                  legend=paste(titleText, groupSuffix[i]),
                  adj=c(0,0.5),
                  bg=titleBoxColor[i],
                  title.col=titleBoxTextColor,
                  text.col=titleBoxTextColor,
                  title=maintitle);
            } else {
               ## Use drawLabels() for more control over the text box
               drawLabels(preset="top",
                  txt=paste(titleText, groupSuffix[i]),
                  boxColor=titleBoxColor[i],
                  boxBorderColor=jamba::makeColorDarker(titleBoxColor[i]),
                  labelCol=titleColor[i],
                  labelCex=titleCex[i],
                  drawBox=doTitleBox,
                  boxCexAdjust=c(1.1,1.3),
                  font=titleFont[i]);
            }
            ## Subtitle using tricks to center the label
            if (length(subtitle) > 0) {
               if (1 == 2) {
                  b1 <- legend("bottom", inset=0.02, plot=FALSE,
                     cex=titleCex[i],
                     text.font=titleFont[i],
                     legend=subtitle,
                     title=NULL);
                  rect(xleft=b1$rect$left,
                     ybottom=b1$rect$top-b1$rect$h,
                     xright=b1$rect$left+b1$rect$w,
                     ytop=b1$rect$top,
                     col=titleBoxColor[i],
                     border="black");
                  text(x=b1$rect$left+b1$rect$w/2,
                     y=b1$rect$top-b1$rect$h/2,
                     labels=subtitle,
                     cex=titleCex[i],
                     font=titleFont[i],
                     col=titleBoxTextColor);
               } else {
                  ## Use drawLabels() for more control over the text box
                  drawLabels(preset="bottom",
                     txt=subtitle[i],
                     boxColor=titleBoxColor[i],
                     boxBorderColor=jamba::makeColorDarker(titleBoxColor[i]),
                     labelCol=titleColor[i],
                     labelCex=titleCex[i],
                     drawBox=doTitleBox,
                     boxCexAdjust=c(1.1,1.3),
                     font=titleFont[i]);
               }
            }
            par("xpd"=parXpd);
         } else {
            titleLine <- margins[3] - 1.5;
            title(main=maintitle,
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
               font.main=titleFont[i], ...);
         }

         ## Optionally print the MAD factor in the bottom right corner
         if (length(outlierMAD) > 0 && displayMAD == 1) {
            legend("bottomright", inset=0.02, cex=titleCex[i],
               text.font=titleFont[i],
               legend=paste0("MAD:",
                  format(digits=2, mvaMADfactors[as.character(i)])),
               bg="transparent", box.lty=0,
               text.col=ifelse(mvaMADfactors[as.character(i)] > outlierMAD,
                  "red3", "grey40"),
               title.col=ifelse(mvaMADfactors[as.character(i)] > outlierMAD,
                  "red3", "grey40"));
         } else if (length(outlierMAD) > 0 && displayMAD == 2) {
            legend("bottomright", inset=0.02, cex=titleCex[i],
               text.font=titleFont[i],
               legend=paste0("MAD:", format(digits=2, mvaMADs[as.character(i)])),
               bg="transparent", box.lty=0,
               text.col=ifelse(mvaMADfactors[as.character(i)] > outlierMAD,
                  "red3", "grey40"),
               title.col=ifelse(mvaMADfactors[as.character(i)] > outlierMAD,
                  "red3", "grey40"));
         }
         mvaData;
      }

      ## End of the per-panel MVA plot loop
      ## Now check to see if we should pad blank panels at the end of the sequence
      if (length(blankPlotPos) > 0) {
         iPanelNumber <- iPanelNumber + 1;
         if (verbose) {
            jamba::printDebug("iPanelNumber: ", iPanelNumber,
               fgText=c("orange", "lightgreen"));
         }
         while(iPanelNumber %in% blankPlotPos) {
            nullPlot(doBoxes=FALSE);
            if (verbose) {
               jamba::printDebug("Inserted a blank panel at position: ", iPanelNumber,
                  fgText=c("orange", "lightblue"));
            }
            iPanelNumber <- iPanelNumber + 1;
         }
      }
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
#' @param indata numeric matrix typically containing measurements (genes)
#'    as rows, and samples as columns.
#' @param floor optional numeric floor, below which values are set to the
#'    floor value, useful when one wants to avoid centering to values which
#'    may be below a noise threshold which might otherwise result in
#'    artificially inflated log fold changes.
#' @param controlSamples optional character vector of colnames(indata) to be
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
#'    it sets needsLog=TRUE and uses log2(indata) for centering.
#' @param mean logical indicating whether to use row means, or row medians.
#'    If the matrixStats package is available, it uses
#'    \code{\link[matrixStats]{rowMedians}} for calculations, otherwise
#'    falling back to apply(indata, 1, median) which is notably slower for
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
#' @export
centerGeneData <- function
(indata,
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
   ## centerGroups can be a list of vectors, containing colnames(indata),
   ## or it can be a vector of group names, in order of colnames(indata).
   ##
   ## If centerGroups is a named vector, and all colnames(indata) are contained in named(centerGroups)
   ## then centerGroups[colnames(indata)] will be used, to ensure they are both
   ## in consistent order, and to allow using a subset of indata accordingly.
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
   if (!returnValues && !returnGroupedValues) {
      stop("One must be TRUE: returnValues or returnGroupedValues.");
   }
   useMS <- FALSE;
   if (suppressPackageStartupMessages(require(matrixStats, quietly=TRUE))) {
      useMS <- TRUE;
   }
   hasCharColumns <- FALSE;
   if (!class(indata) %in% c("matrix", "numeric")) {
      colClass <- sapply(1:ncol(indata), function(i){class(indata[,i])});
      colClassChar <- which(!colClass %in% c("numeric", "integer"));
      if (length(colClassChar) > 0) {
         hasCharColumns <- TRUE;
         indataChar <- indata[,colClassChar,drop=FALSE];
         indata <- as.matrix(indata[,-colClassChar]);
      }
   }
   #printDebug("controlSamples:", controlSamples, c("orange", "lightgreen"));
   if (missing(controlSamples) || length(controlSamples) == 0
       || is.na(controlSamples) || is.null(controlSamples)) {
      controls <- rep(TRUE, ncol(indata));
   } else {
      if (class(controlSamples) %in% c("integer", "numeric") && any(controlSamples <= ncol(indata))) {
         controls <- rep(FALSE, ncol(indata));
         controls[controlSamples[controlSamples <= ncol(indata)]] <- TRUE;
      } else if (any(controlSamples %in% colnames(indata))) {
         controls <- (colnames(indata) %in% controlSamples);
      } else {
         controls <- rep(TRUE, ncol(indata));
      }
   }
   controlCols <- colnames(indata)[controls];

   ## We will try to auto-detect whether to log2 transform the data
   if (length(needsLog) == 0) {
      needsLog <- max(indata, na.rm=TRUE) > 100;
   }
   if (needsLog) {
      if (verbose) {
         jamba::printDebug("centerGeneData(): ",
            "applying log2 transform:",
            "log(1 + x)");
      }
      indata <- log2(1+indata);
   }

   ## Create summary data.frame
   centerDF <- as.data.frame(rmNULL(list(sample=colnames(indata),
      centerGroups=centerGroups,
      centerControl=colnames(indata) %in% controlCols)));
   #centerDF <- mmixedOrderDF(centerDF, byCols=c(2,-3));
   if (verbose) {
      jamba::printDebug("centerGeneData(): ",
         "dim(indata)", dim(indata));
      jamba::printDebug("centerGeneData(): ",
         "dim(centerDF)", dim(centerDF));
      print(head(centerDF));
   }
   rownames(centerDF) <- centerDF[,"sample"];

   ## optionally show summary for visual confirmation
   if (showGroups && length(centerGroups) > 0) {
      print(centerDF);
      ## We avoid including ch() and colsHead() for now.
      #ch(centerDF, maxRows=nrow(centerDF), orderBy=c(-2,3));
   }

   ## Optionally center groups separately
   if (length(centerGroups) > 0) {
      if (!jamba::igrepHas("list", class(centerGroups))) {
         if (length(names(centerGroups)) > 0) {
            if (all(colnames(indata) %in% names(centerGroups))) {
               centerGroups <- centerGroups[colnames(indata)];
            } else {
               stop("colnames(indata) must be contained in names(centerGroups) when centerGroups has names.");
            }
         }
         centerGroups <- split(colnames(indata), centerGroups);
      }
      ## Now iterate through the list
      #centeredSubsets <- lapply(centerGroups, function(iGroup){
      centeredSubsets <- lapply(nameVectorN(centerGroups), function(iGroupN){
         iGroup <- centerGroups[[iGroupN]];
         iGroupCols <- colnames(indata[,iGroup,drop=FALSE]);
         iControls <- intersect(iGroupCols, controlCols);
         iM <- centerGeneData(indata[,iGroup,drop=FALSE],
            controlSamples=iControls,
            needsLog=needsLog,
            floor=floor,
            mean=mean,
            returnGroupedValues=returnGroupedValues,
            returnValues=returnValues,
            groupPrefix=iGroupN,
            centerGroups=NULL,
            scale=scale);
         #if (!returnValues) {
         #   iM <- iM[,!colnames(iM) %in% iGroup,drop=FALSE];
         #}
         iM;
      });
      centeredData <- do.call(cbind, centeredSubsets);
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
      indataMeans <- rowMeans(indata[,controls,drop=FALSE]);
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
      jamba::printDebug("Keeping non-numeric columns as-is.");
      centeredData <- cbind(indataChar, centeredData)
   }
   attr(centeredData, "centerDF") <- centerDF;
   return(centeredData);
}


#' define a polygon hull around points
#'
#' Defines a polygon hull to encompass a set of data points
#'
#' A polygon hull can be useful as a visual indicator of the location
#' of a set of points. This function creates a polygon around all points,
#' and not just the most dense region of points, using
#' `rgeos::gConvexHull()`. Since that function returns
#' an object of class SpatialPolygons, this function is mainly a wrapper
#' which returns the set of coordinates, either as a matrix or as the
#' SpatialPolygons object type.
#'
#' @param x two-column numeric matrix of data points, or an object of
#'    class SpatialPoints.
#' @param returnClass character class name to return, either matrix or
#'    SpatialPolygons.
#' @param verbose logical whether to print verbose output'
#'
#' @return
#' A numeric matrix of polygon coordinates, unless the parameter
#' `returnClass="SpatialPolygons"` is set, in which case it returns
#' a `sp::SpatialPolygons` object.
#'
#' @seealso `sp::SpatialPolygons`, `rgeos::gConvexHull()`
#'
#' @family jam spatial functions
#'
#' @export
points2polygonHull <- function
(x,
 returnClass=c("matrix","SpatialPolygons"),
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
   if (!suppressPackageStartupMessages(require(sp))) {
      stop("points2polygonHull() requires the sp package SpatialPoints class.");
   }
   if (!suppressPackageStartupMessages(require(rgeos))) {
      stop("points2polygonHull() requires the rgeos package gConvexHull().");
   }
   if (jamba::igrepHas("matrix|data.*frame|tbl", class(x))) {
      if (ncol(x) != 2) {
         stop(paste0("points2polygonHull() requires x as a",
            " matrix,data.frame,tbl",
            " to contain 2 columns."));
      }
      x <- SpatialPoints(x);
   } else if (jamba::igrepHas("SpatialPoints", class(x))) {
      #SpatialPoints(x)
   }

   ## Determine convex hull coordinates
   gH <- gConvexHull(x);
   if (returnClass %in% "matrix") {
      if (jamba::igrepHas("polygons", class(gH))) {
         gHcoords <- gH@polygons[[1]]@Polygons[[1]]@coords;
      } else if (jamba::igrepHas("lines", class(gH))) {
         gHcoords <- gH@lines[[1]]@Lines[[1]]@coords;
      } else if (jamba::igrepHas("points", class(gH))) {
         gHcoords <- gH@coords;
      }
   } else {
      gHcoords <- gH;
   }
   gHcoords;
}

#' Draw text labels on a base R plot
#'
#' Draw text labels on a base R plot
#'
#' This function takes a vector of coordinates and text labels,
#' and draws the labels with colored rectangles around each label
#' on the plot. Each label can have unique font, cex, and color,
#' and are drawn using vectorized operations.
#'
#' TODO: In future allow rotated text labels. Not that useful within
#' a plot panel, but sometimes useful when draw outside a plot, for
#' example axis labels.
#'
#' @return invisible data.frame containing label coordinates used
#' to draw labels. This data.frame can be manipulated and provided
#' as input to `drawLabels()` for subsequent customized label
#' positioning.
#'
#' @family jam plot functions
#'
#' @param newCoords optional `data.frame` typically as a result of
#'    a previous call to `drawLabels()`. In general, it should contain
#'    colnames equivalent to the function parameters of `drawLabels()`.
#' @param x,y vector of x- and y- coordinates.
#' @param txt character vector of labels, length equal to `x` and `y`.
#' @param lx,ly optional vector of segment endpoint coordinates, used
#'    to draw a line from x,y coordinates to the segment lx,ly coordinate.
#' @param segmentLwd,segmentCol vector of segment line widths, and colors,
#'    respectively. Each vector will be recycled to `length(txt)` as
#'    needed.
#' @param drawSegments logical whether to draw segments, where applicable.
#' @param boxBorderColor vector of colors used for the box border around
#'    each label.
#' @param boxColor vector of colors used for the box background behind
#'    each label.
#' @param boxLwd vector of box line widths, sent to `graphics::rect()`,
#'    this vector will be recycled to `length(txt)`.
#' @param drawBox logical whether to draw boxes behind each text label.
#' @param drawLabels logical whether to draw each text label.
#' @param font vector of font values as described in `graphics::par()`,
#'    where 1=normal, 2=bold, 3=italics, 4=bold-italics.
#' @param labelCex vector of cex values used for text labels. This vector
#'    will be recycled to `length(txt)` as needed.
#' @param boxCexAdjust numeric vector length=2, used to expand the x-width
#'    and y-height of the box around around text labels.
#' @param labelCol vector of label colors, by default it calls
#'    `jamba::setTextContrastColor()` to generate a color to contrast
#'    the background box color.
#' @param doPlot logical whether to perform any plot operations. Set
#'    `FALSE` to calculate coordinates and return a `data.frame` of
#'    label coordinates, which can then be manipulated before calling
#'    `drawLabels()` again.
#' @param xpd value compatible with `par("xpd")`, where NA allows labels
#'    anywhere in the device region, TRUE retricts labels within the figure
#'    region, and FALSE restricts labels within the plot region.
#' @param preset vector of text values used to position text labels relative
#'    to the x,y coordinate, where "topleft" will position the label so the
#'    entire label box is top-left of the point, therefore the point will be
#'    at the bottom-right corner of the label box. When `preset` is anything
#'    by `"none"` the `adjX` and `adjY` values are ignored.
#' @param adjX,adjY the text adjustment of labels relative to the x,y
#'    coordinate. The values are recycled to `length(txt)`.
#' @param trimReturns logical whether to trim leading and trailing return
#'    (newline) characters from labels.
#' @param verbose logical whether to print verbose output.
#' @param ... additional arguments are passed to `graphics::segments()` when
#'    segments are drawn, to `graphics::rect()` when label boxes are drawn,
#'    and to `graphics::text()` when text labels are drawn.
#'
#' @examples
#' jamba::nullPlot();
#' drawLabels(x=par("usr")[1],
#'    y=par("usr")[4],
#'    txt="Top-left\nof plot",
#'    preset="bottomright",
#'    boxColor="blue4");
#'
#' drawLabels(x=par("usr")[2],
#'    y=par("usr")[3],
#'    txt="Bottom-right\nof plot",
#'    preset="topleft",
#'    boxColor="green4");
#'
#' drawLabels(x=mean(par("usr")[1:2]),
#'    y=mean(par("usr")[3:4]),
#'    txt="Center\nof plot",
#'    preset="center",
#'    boxColor="purple");
#'
#' points(x=c(par("usr")[1], par("usr")[2],
#'       mean(par("usr")[1:2])),
#'    y=c(par("usr")[4], par("usr")[3],
#'       mean(par("usr")[3:4])),
#'    pch=20,
#'    col="red",
#'    xpd=NA);
#'
#' @export
drawLabels <- function
(newCoords=NULL,
 x=NULL,
 y=NULL,
 txt=NULL,
 lx=NULL,
 ly=NULL,
 segmentLwd=1,
 segmentCol="#00000088",
 drawSegments=TRUE,
 boxBorderColor="#000000AA",
 boxColor="#DDAA77",
 boxLwd=1,
 drawBox=TRUE,
 drawLabels=TRUE,
 font=1,
 labelCex=0.8,
 boxCexAdjust=c(1,1.2),
 labelCol=jamba::alpha2col(alpha=0.8, jamba::setTextContrastColor(boxColor)),
 doPlot=TRUE,
 xpd=NA,
 preset=c("default",
   "top",
   "topright",
   "right",
   "bottomright",
   "topleft",
   "left",
   "bottomleft",
   "bottom",
   "center"),
 adjPreset=c("default",
   "top",
   "topright",
   "right",
   "bottomright",
   "topleft",
   "left",
   "bottomleft",
   "bottom",
   "center"),
 adjX=0.5,
 adjY=0.5,
 trimReturns=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to wrapper only the last portion of addNonOverlappingLabels()
   ## which draws the labels, boxes, and segments after positions are determined
   ## by addNonOverlappingLabels().
   if (length(boxCexAdjust) == 0) {
      boxCexAdjust <- 1;
   }
   boxCexAdjust <- rep(boxCexAdjust, length.out=2);
   preset <- match.arg(preset);
   adjPreset <- match.arg(adjPreset);


   if (length(newCoords) == 0) {
      ## Create a basic data.frame
      if (jamba::igrepHas("top|bottom|left|right|center", c(preset, adjPreset))) {
         if (verbose) {
            jamba::printDebug("drawLabels(): ",
               "Processing non-default preset and adjPreset values.");
         }
         presetL <- coordPresets(
            preset=preset,
            x=x,
            y=y,
            adjPreset=adjPreset,
            adjX=adjX,
            adjY=adjY);
         x <- presetL$x;
         y <- presetL$y;
         adjX <- presetL$adjX;
         adjY <- presetL$adjY;
         preset <- presetL$preset;
         adjPreset <- presetL$adjPreset;
      }
      newCoords <- data.frame(
         x=x,
         y=y,
         txt=txt,
         w=strwidth(txt,
            font=font,
            cex=labelCex),
         h=strheight(txt,
            font=font,
            cex=labelCex),
         hNudge=0,
         labelCex=labelCex,
         adjX=adjX,
         adjY=adjY,
         preset=preset,
         adjPreset=adjPreset);
      if (length(lx) > 0 && length(ly) > 0) {
         newCoords$lx <- rep(lx, length.out=nrow(newCoords));
         newCoords$ly <- rep(ly, length.out=nrow(newCoords));
      } else {
         drawSegments <- FALSE;
      }

      ## Make additional adjustments based upon leading or trailing "\n"
      if (trimReturns) {
         ## First determine the height of the buffer between two lines of text
         labelCexU <- unique(newCoords$labelCex);
         labelCexUh <- sapply(labelCexU, function(iCex){
            bufferH <- (strheight("|\n|", cex=iCex) -
                  2*(strheight("|", cex=iCex)));
         });
         #labelCexUh[match(newCoords$labelCex, labelCexU)]
         if (jamba::igrepHas("\n$", newCoords$txt)) {
            if (verbose) {
               jamba::printDebug("drawLabels(): ",
                  "adjusting adjY for trailing return characters.");
               print(newCoords);
            }
            while(jamba::igrepHas("\n$", newCoords$txt)) {
               iAdj <- jamba::igrep("\n$", newCoords$txt);
               adjY[iAdj] <- adjY[iAdj] - newCoords$adjY[iAdj];
               newCoords$txt[iAdj] <- sub("\n$", "", newCoords$txt[iAdj]);
            }
         }
         if (jamba::igrepHas("^\n", newCoords$txt)) {
            if (verbose) {
               jamba::printDebug("drawLabels(): ",
                  "adjusting adjY for leading return characters.");
            }
            while(jamba::igrepHas("^\n", newCoords$txt)) {
               iAdj <- jamba::igrep("^\n", newCoords$txt);
               adjY[iAdj] <- adjY[iAdj] + newCoords$adjY[iAdj];
               newCoords$txt[iAdj] <- sub("^\n", "", newCoords$txt[iAdj]);
            }
         }
         newCoords$adjY <- adjY;
         newCoords$h <- strheight(newCoords$txt, cex=labelCex);
      } else {
         newCoords$txt <- gsub("\n", "\n|", newCoords$txt);
      }

      ## Now widen the label box using boxCexAdjust
      if (verbose) {
         jamba::printDebug("drawLabels(): ",
            "adjusting w,h using boxCexAdjust");
         print(newCoords);
      }
      ## New strategy intended to keep the bottom-left edge fixed
      newCoords$h <- newCoords$h * boxCexAdjust[2];
      newCoords$w <- newCoords$w * boxCexAdjust[1];
      newCoords$x <- newCoords$x - adjX * newCoords$w;
      newCoords$y <- newCoords$y - adjY * newCoords$h;
   } else {
      preset <- newCoords$preset;
      adjPreset <- newCoords$adjPreset;
      if (jamba::igrepHas("top|bottom|left|right|center", c(preset, adjPreset))) {
         if (verbose) {
            jamba::printDebug("drawLabels(): ",
               "Processing non-default preset and adjPreset values.");
         }
         presetL <- coordPresets(
            preset=preset,
            x=newCoords$x,
            y=newCoords$y,
            adjPreset=adjPreset,
            adjX=newCoords$adjX,
            adjY=newCoords$adjY);
         x <- presetL$x;
         y <- presetL$y;
         adjX <- presetL$adjX;
         adjY <- presetL$adjY;
         preset <- presetL$preset;
         adjPreset <- presetL$adjPreset;
      }
      newCoords$x <- x;
      newCoords$y <- y;
      newCoords$adjX <- adjX;
      newCoords$adjY <- adjY;
      newCoords$preset <- preset;
      newCoords$adjPreset <- adjPreset;
   }

   ## Add height and width if not supplied
   if (!"h" %in% names(newCoords)) {
      newCoords$h <- strheight(txt, cex=labelCex*1.1);
   }
   if (!"w" %in% names(newCoords)) {
      newCoords$w <- strwidth(txt, cex=labelCex*1.1);
   }

   x <- newCoords$ptX;
   y <- newCoords$ptY;
   if (!"labelCex" %in% names(newCoords)) {
      if (verbose) {
         jamba::printDebug("drawLabels(): ",
            "Defining all fixedLabels=FALSE");
      }
      newCoords[,"fixedLabels"] <- FALSE;
   }
   if (!"fixedLabels" %in% names(newCoords)) {
      newCoords[,"fixedLabels"] <- TRUE;
   }
   whichLabels <- which(!newCoords$fixedLabels > 1);
   if (verbose) {
      jamba::printDebug("drawLabels(): ",
         "head(whichLabels):",
         head(whichLabels));
   }

   segmentCol <- rep(segmentCol, length.out=nrow(newCoords));
   segmentLwd <- rep(segmentLwd, length.out=nrow(newCoords));

   if (!"labelCex" %in% names(newCoords)) {
      newCoords$labelCex <- rep(labelCex,
         length.out=nrow(newCoords));
   }
   if (!"labelCol" %in% names(newCoords)) {
      newCoords$labelCol <- rep(labelCol,
         length.out=nrow(newCoords));
   }

   if (!"boxColor" %in% names(newCoords)) {
      newCoords$boxColor <- rep(boxColor,
         length.out=nrow(newCoords));
   }
   boxColor <- newCoords$boxColor;
   if (!"boxBorderColor" %in% names(newCoords)) {
      newCoords$boxBorderColor <- rep(boxBorderColor,
         length.out=nrow(newCoords));
   }
   boxBorderColor <- newCoords$boxBorderColor;

   ## Draw segments first so the box and labels will be printed on top of them
   if (doPlot &&
         drawSegments &&
         any(whichLabels) &&
         all(c("lx","ly") %in% names(newCoords))) {
      ## Determine closest corner or flat side for each segment
      x0l <- newCoords$x[whichLabels];
      x0m <- newCoords$x[whichLabels] + newCoords$w[whichLabels]/2;
      x0r <- newCoords$x[whichLabels] + newCoords$w[whichLabels];
      y0b <- newCoords$y[whichLabels];
      y0m <- newCoords$y[whichLabels] + newCoords$h[whichLabels]/2;
      y0t <- newCoords$y[whichLabels] + newCoords$h[whichLabels];
      x1 <- newCoords$lx[whichLabels];
      y1 <- newCoords$ly[whichLabels];
      sx0 <- ifelse(x0r < x1, x0r,
         ifelse(x0l > x1, x0l,
            x1));
      sy0 <- ifelse(y0t < y1, y0t,
         ifelse(y0b > y1, y0b,
            y1));
      if (verbose) {
         headN <- 10;
         jamba::printDebug("drawLabels(): ",
            "Drawing line segments.");
         jamba::printDebug("sx0: ",
            head(sx0, headN));
         jamba::printDebug("x1: ",
            head(x1, headN));
         jamba::printDebug("sy0: ",
            head(sy0, headN));
         jamba::printDebug("y1: ",
            head(y1, headN));
         jamba::printDebug("segmentCol[whichLabels]: ",
            head(c(segmentCol[whichLabels]), headN), c("orange", "lightblue"));
         jamba::printDebug("segmentLwd[whichLabels]: ",
            head(c(segmentLwd[whichLabels]), headN), c("orange", "lightblue"));
      }
      segments(x0=sx0,
         y0=sy0,
         x1=newCoords$lx[whichLabels],
         y1=newCoords$ly[whichLabels],
         col=segmentCol[whichLabels],
         lwd=segmentLwd[whichLabels],
         xpd=xpd,
         ...);
   } else {
      if (verbose) {
         jamba::printDebug("drawLabels(): ",
            "Drawing no segments.");
      }
   }

   ## Due to a small bug in add.labels in the snippets package, the box and text
   ## labels are not properly adjusted when the adj(1,1) parameter is supplied,
   ## causing the boxes to be slightly offset in the y-axis direction relative
   ## to the text.
   ## The code below calculates the text to be the exact center of each box,
   ## then forces the text() method below to use adj=c(0.5,0.5) which centers
   ## the text exactly at this coordinate. Looks much nicer than other
   ## options.
   boxX1 <- newCoords$x;
   boxX2 <- boxX1 + newCoords$w;
   boxY1 <- newCoords$y - newCoords$hNudge;
   boxY2 <- boxY1 + newCoords$h;
   textX <- (boxX1 + boxX2) / 2;
   textY <- (boxY1 + boxY2) / 2;

   ## Draw boxes before labels so the labels will be printed on top
   if (drawBox && doPlot && any(whichLabels)) {
      if (verbose) {
         jamba::printDebug("drawLabels(): ",
            "Drawing rectangles.");
         jamba::printDebug("drawLabels(): ",
            "First rectangle:",
            c(boxX1[whichLabels][1], boxX2[whichLabels][1],
               boxY1[whichLabels][1], boxY2[whichLabels][1]));
      }
      ## Only draw boxes where there are characters to be printed
      rect(xleft=boxX1[whichLabels],
         ybottom=boxY1[whichLabels],
         xright=boxX2[whichLabels],
         ytop=boxY2[whichLabels],
         col=boxColor[whichLabels],
         border=boxBorderColor[whichLabels],
         lwd=boxLwd,
         xpd=xpd,
         ...);
   }
   if (drawLabels && doPlot && any(whichLabels)) {
      if (verbose) {
         jamba::printDebug("drawLabels(): ",
            "Printing labels head(txt): ",
            paste(head(txt), collapse=", "));
      }
      text(x=textX[whichLabels],
         y=textY[whichLabels],
         font=font,
         labels=newCoords$txt[whichLabels],
         col=newCoords$labelCol[whichLabels],
         cex=newCoords$labelCex[whichLabels],
         outer=TRUE,
         adj=c(0.5,0.5),
         xpd=xpd,
         ...);
   }
   invisible(newCoords);
}

#' Process coordinate adjustment presets
#'
#' Process coordinate adjustment presets
#'
#' This function is intended to be a convenient way to define
#' coordinates using preset terms like "topleft", "bottom", "center".
#'
#' Similarly, it is intended to help define corresponding text
#' adjustments, using `adj` compatible with `graphics::text()`,
#' using preset terms like "bottomright", "center".
#'
#' When `preset` is `"default"`, the original `x,y` coordinates
#' are used. Otherwise the `x,y` coordinates are defined using the
#' plot region coordinates, where `"left"` uses `par("usr")[1]`,
#' and `"top"` uses `par("usr")[4]`.
#'
#' When `adjPreset` is `"default"` it will use the `preset` to
#' define a reciprocal text placement. For example when `preset="topright"`
#' the text placement will be equivalent to `adjPreset="bottomleft"`.
#' The `adjPreset` terms `"top"`, `"bottom"`, `"right"`, `"left"`,
#' and `"center"` refer to the text label placement relative to
#' `x,y` coordinate.
#'
#' If both `preset="default"` and `adjPreset="default"` the original
#' `adjX,adjY` values are returned.
#'
#' The function is vectorized, and uses the longest input argument,
#' so one can supply a vector of `preset` and it will return coordinates
#' and adjustments of length equal to the input `preset` vector.
#' The `preset` value takes priority over the supplied `x,y` coordinates.
#'
#' @return data.frame with colnames `x,y,adjX,adjY,preset,adjPreset`
#' after adjustment, where the number of rows is determined by the
#' longest input argument.
#'
#' @param preset character vector of coordinate positions, or "default"
#'    to use the `x,y` coordinates.
#' @param x,y numeric vectors indicating the default coordinates `x,y`.
#' @param adjPreset character vector of text label positions, or
#'    "default" to use `preset`, or when `preset="default"` the
#'    `adjX,adjY` values are used.
#' @param adjX,adjY numeric vectors indicating default text adjustment
#'    values, as described for `adj` in `graphics::text()`.
#' @param adjOffsetX,adjOffsetY numeric vector used to apply an offset
#'    value to the `adjX,adjY` values, where positive values would
#'    place a label farther away from center. Note these units are
#'    relative to the text label size, when used with `graphics::text()`,
#'    larger labels will be adjusted more than smaller labels.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' # make sure to prepare the plot region first
#' jamba::nullPlot(plotAreaTitle="");
#'
#' # determine coordinates
#' presetV <- c("top","bottom","left","right", "topleft");
#' cp1 <- coordPresets(preset=presetV);
#' cp1;
#'
#' points(cp1$x, cp1$y, pch=20, cex=2, col="red");
#'
#' # unfortunately graphics::text() does not have vectorized adj
#' for (i in seq_along(presetV)) {
#'    text(cp1$x[i], cp1$y[i],
#'       labels=presetV[i],
#'       adj=c(cp1$adjX[i], cp1$adjY[i]));
#' }
#'
#' # drawLabels() is vectorized for example
#' jamba::nullPlot(plotAreaTitle="");
#' presetV2 <- c("topleft", "topright", "bottomleft", "bottomright");
#' cp2 <- coordPresets(preset="center", adjPreset=presetV2, adjOffsetX=0.1, adjOffsetY=0.4);
#' points(cp2$x, cp2$y, pch=20, cex=2, col="red");
#' drawLabels(x=cp2$x, y=cp2$y, adjX=cp2$adjX, adjY=cp2$adjY, txt=presetV2, boxCexAdjust=c(1.15,1.6), labelCex=1.3, lx=rep(1.5, 4), ly=rep(1.5, 4))
#'
#' @export
coordPresets <- function
(preset="default",
 x=0,
 y=0,
 adjPreset="default",
 adjX=0.5,
 adjY=0.5,
 adjOffsetX=0,
 adjOffsetY=0,
 verbose=FALSE,
 ...)
{
   ## Takes vector of preset values, and x,y
   ## and re-defines x,y for anything with
   ## top, left, right, bottom, center
   ## Takes vector of adjPreset values, and adjX, adjY
   ## and re-defines adjX,adjY.
   ## When adjPreset="default" it uses opposite orientation
   ## of "preset", otherwise position is defined by adjPreset
   ## for top, bottom, left, right, center
   n <- max(lengths(list(
      x,
      y,
      adjX,
      adjY,
      preset,
      adjPreset
   )));
   x <- rep(x, length.out=n);
   y <- rep(y, length.out=n);
   preset <- rep(preset, length.out=n);
   adjPreset <- rep(adjPreset, length.out=n);
   adjX <- rep(adjX, length.out=n);
   adjY <- rep(adjY, length.out=n);
   if (verbose) {
      jamba::printDebug("coordPresets(): ",
         "n:",
         n);
      jamba::printDebug("coordPresets(): ",
         "preset:",
         preset);
      jamba::printDebug("coordPresets(): ",
         "adjPreset:",
         adjPreset);
   }

   ## Process the preset logic
   x <- ifelse(grepl("right", preset), par("usr")[2],
      ifelse(grepl("left", preset), par("usr")[1],
         ifelse(grepl("center|top|bottom", preset), mean(par("usr")[1:2]),
            x)));
   if (verbose) {
      jamba::printDebug("coordPresets(): ",
         "x:",
         x);
   }
   adjX <- ifelse(grepl("default", adjPreset),
      ifelse(grepl("left", preset), 0-adjOffsetX,
         ifelse(grepl("right", preset), 1+adjOffsetX,
            ifelse(grepl("center|top|bottom", preset), 0.5,
               adjX))),
      ifelse(grepl("right", adjPreset), 0-adjOffsetX,
         ifelse(grepl("left", adjPreset), 1+adjOffsetX,
            ifelse(grepl("center|top|bottom", adjPreset), 0.5,
               adjX))));
   y <- ifelse(grepl("top", preset), par("usr")[4],
      ifelse(grepl("bottom", preset), par("usr")[3],
         ifelse(grepl("center|left|right", preset), mean(par("usr")[3:4]),
            y)));
   adjY <- ifelse(grepl("default", adjPreset),
      ifelse(grepl("top", preset), 1+adjOffsetY,
         ifelse(grepl("bottom", preset), 0-adjOffsetY,
            ifelse(grepl("center|left|right", preset), 0.5,
               adjY))),
      ifelse(grepl("top", adjPreset), 0-adjOffsetY,
         ifelse(grepl("bottom", adjPreset), 1+adjOffsetY,
            ifelse(grepl("center|left|right", adjPreset), 0.5,
               adjY))));
   return(data.frame(x=x,
      y=y,
      adjX=adjX,
      adjY=adjY,
      preset=preset,
      adjPreset=adjPreset));
}
