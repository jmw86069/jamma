
#' @describeIn jammaplot
#'
#' @family jam plot functions
#'
#' @examples
#' if (jamba::check_pkg_installed("SummarizedExperiment") &&
#'    jamba::check_pkg_installed("farrisdata")) {
#'    suppressPackageStartupMessages(require(SummarizedExperiment));
#'
#'    GeneSE <- farrisdata::farrisGeneSE;
#'
#'    titleBoxColor <- jamba::nameVector(
#'       farrisdata::colorSub[as.character(
#'          SummarizedExperiment::colData(GeneSE)$groupName)],
#'       colnames(GeneSE));
#'    options("warn"=FALSE);
#'
#'    gg <- ggjammaplot(GeneSE,
#'       ncol=6,
#'       base_size=12,
#'       maintitle="Farris raw RNAseq data",
#'       titleBoxColor=jamba::rgb2col(col2rgb(titleBoxColor)),
#'       assay_name="raw_counts")
#'
#'    gg <- ggjammaplot(GeneSE,
#'       ncol=6,
#'       assay_name="counts",
#'       useRank=TRUE,
#'       ylim=c(-11000, 11000),
#'       maintitle="MA-plots by rank and rank difference",
#'       titleBoxColor=titleBoxColor)
#'
#'    gg <- ggjammaplot(GeneSE,
#'       ncol=6,
#'       assay_name="counts",
#'       titleBoxColor=titleBoxColor,
#'       base_size=10,
#'       maintitle="MA-plots showing MAD factor",
#'       displayMAD=TRUE)
#'
#'    gg <- ggjammaplot(GeneSE,
#'       ncol=6,
#'       assay_name="counts",
#'       titleBoxColor=titleBoxColor,
#'       maintitle="MA-plot omitting one panel, then using blankPlotPos",
#'       whichSamples=colnames(GeneSE)[c(1:21, 23:24)],
#'       blankPlotPos=22,
#'       displayMAD=TRUE)
#'
#'    if (FALSE) {
#'    ggdf <- ggjammaplot(GeneSE,
#'       assay_name="counts",
#'       whichSamples=c(1:3, 7:9),
#'       return_type="data",
#'       titleBoxColor=titleBoxColor)
#'    highlightPoints1 <- names(jamba::tcount(subset(ggdf, mean > 15 & difference < -1)$item, 2))
#'    highlightPoints2 <- subset(ggdf, name %in% "CA1CB492" &
#'       difference < -4.5)$item;
#'    highlightPoints <- list(
#'       divergent=highlightPoints1,
#'       low_CA1CB492=highlightPoints2);
#'
#'    ggdf_h <- ggjammaplot(GeneSE,
#'       assay_name="counts",
#'       highlightPoints=highlightPoints,
#'       whichSamples=c(1:3, 7:9),
#'       return_type="data",
#'       titleBoxColor=titleBoxColor)
#'
#'    # you can use output from `jammaplot()` as input to `ggjammaplot()`:
#'    jp2 <- jammaplot(GeneSE,
#'       outlierMAD=2,
#'       doPlot=FALSE,
#'       assay_name="raw_counts",
#'       filterFloor=1e-10,
#'       filterFloorReplacement=NA,
#'       centerGroups=colData(GeneSE)$Compartment,
#'       subtitleBoxColor=farrisdata::colorSub[as.character(colData(GeneSE)$Compartment)],
#'       useRank=FALSE);
#'
#'    gg1 <- ggjammaplot(jp2,
#'       ncol=6,
#'       titleBoxColor=titleBoxColor);
#'    print(gg1);
#'    }
#' }
#'
#' @inheritParams jammaplot
#' @param x one of the following inputs:
#'    * `numeric` matrix
#'    * `SummarizedExperiment` object, where the
#'    assay data is defined using `assays(x)[[assay_name]]`. Accordingly,
#'    `assay_name` can be either an integer index or `character` string
#'    matching the `names(assays(x))`.
#'    * `list` output from `jammacalc()` or `jammaplot()`, where each
#'    element in the list is a two-column `matrix` with colnames `c("x", "y")`.
#' @param detail_factor `numeric` used to adjust the level
#'    of detail, as a multiplier for `nbin_factor` and `bw_factor`.
#' @param nbin_factor `numeric` value used to adjust the number of bins
#'    used to display the MA-plots, where values higher than `1` increase
#'    the resolution and level of detail, and values below `1` decrease
#'    the resolution. Note the number of bins are already adjusted based
#'    upon the square root of the number of plot panels, and `nbin_factor`
#'    applied to that value.
#' @param bw_factor `numeric` used to adjust the resolution of the
#'    2-dimensional bandwidth calculation, where higher values create
#'    more detailed density, and lower values create a smoother density
#'    across the range of data.
#'    In some cases, the `ggplot` panel aspect ratio diverges from 1:1,
#'    in which case `bw_factor` can be used to expand the bandwidth by
#'    y-axis, or y-axis, respectively. For example, if the density
#'    appears short-wide, try `bw_factor=c(1.2, 1)`, if the density
#'    appears tall-skinny, try `bw_factor=c(1, 1.2)`.
#' @param assay_name relevant only when `x` is `SummarizedExperiment`,
#'    one of these input types:
#'    * `character` string that matches `names(assays(x))`
#'    * `integer` index for `assays(x)`, where any value higher than
#'    `length(assays(x))` is adjusted to `length(assays(x))`, which makes
#'    it convenient to select the last element in the list of `assays(x)`
#'    by using `assay_name = Inf`.
#' @param useMedian `logical` indicating whether calculations should use
#'    `median`, or when `useMedian=FALSE` the `mean` is used. The median
#'    has the benefit of reducing effect of outliers, however the mean
#'    has the advantage that it represents data consistent with most
#'    parametric statistical analyses.
#' @param controlSamples `character` vector of `colnames(x)` to use as
#'    the control when calculating centered data. By default, all samples
#'    are used, so the classic MA-plot is the value of each sample,
#'    subtracting the median or mean value calculated across all samples.
#'    It is sometimes useful to define a subset of known samples for this
#'    calculation, which can be beneficial in avoiding outliers, or for
#'    consistency by selecting high quality control samples.
#' @param centerGroups `character` vector with length equal to `ncol(x)`,
#'    which defines subgroups of `colnames(x)` to be treated independently
#'    during the MA-plot calculation.
#' @param groupedX `logical` indicating whether the x-axis value, which
#'    represents the median or mean value, should be calculated independently
#'    for each group when `centerGroups` is used with multiple groups.
#'    Typically `groupedX=TRUE` is recommended, however it can be beneficial
#'    to share an overall x-axis value in specific circumstances.
#' @param grouped_mad `logical` indicating whether the MAD factor calculation
#'    of variability among samples should be performed independently
#'    for each group when `centerGroups` is used with multiple groups.
#'    Typically `grouped_max=TRUE` is recommended, however it can be beneficial
#'    to share an overall MAD factor threshold across all samples
#'    in specific circumstances.
#' @param outlierMAD `numeric` indicating the MAD factor threshold above
#'    which a particular sample is considered an outlier.
#' @param mad_row_min `numeric` value indicating the minimum x-axis
#'    value, calculated using either median or mean as defined by
#'    argument `useMedian`, at or above which a measurement is used in the
#'    MAD factor calculation. This threshold is useful to restrict the
#'    MAD variability calculation to measurements (rows in `x`) with
#'    signal that meets a minimum noise threshold.
#' @param displayMAD `logical` indicating whether to display the MAD factor
#'    in the bottom right corner of each MA-plot panel.
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
#' @param naValue `character` string used to convert values of `NA` to
#'    something else. This argument is useful when a numeric matrix may
#'    contain `NA` values but would prefer them to be, for example, `0`.
#' @param centerFunc `function` used to supply a custom data centering
#'    function. In practice this argument should rarely be changed.
#' @param whichSamples `integer` index of samples in `colnames(x)` to be
#'    plotted, however all samples in `colnames(x)` will be used for the
#'    MA-plot calculations and data centering. This argument is intended
#'    to help zoom in to inspect a specific subset of samples, without
#'    having to plot all samples in `x`.
#' @param useRank `logical` indicating whether to plot rank on the x-axis,
#'    rank-difference on the y-axis for each sample. This transformation
#'    is rather useful, especially when downstream analysis tools may
#'    also refer to the rank value of particular measurements.
#' @param titleBoxColor `character` vector of R colors, where
#'    `titleBoxColor` is equal to `ncol(x)`, or where
#'    `names(titleBoxColor)` matches `colnames(x)`. When supplied, each
#'    plot panel strip background will be colored accordingly.
#'
#' @export
ggjammaplot <- function
(x,
 detail_factor=1,
 nbin_factor=1,
 bw_factor=1,
 assay_name=1,
 useMedian=FALSE,
 controlSamples=NULL,
 centerGroups=NULL,
 controlFloor=NA,
 naControlAction=c("row", "floor", "min", "na"),
 naControlFloor=0,
 colramp=c("transparent",
    "lightblue",
    "blue",
    "navy",
    "orange",
    "orangered2"),
 groupedX=TRUE,
 grouped_mad=TRUE,
 outlierMAD=5,
 mad_row_min=4,
 displayMAD=FALSE,
 noise_floor=0,
 noise_floor_value=NA,
 naValue=NA,
 centerFunc=centerGeneData,
 whichSamples=NULL,
 useRank=FALSE,
 apply_transform_limit=40,
 titleBoxColor="lightgoldenrod1",
 titleCex=1,
 outlierColor="lemonchiffon",
 fillBackground=TRUE,
 maintitle=NULL,
 subtitle=NULL,
 summary="Mean",
 difference="Difference",
 transFactor=0.25,
 doPlot=TRUE,
 highlightPoints=NULL,
 highlightPch=21,
 highlightCex=1.5,
 highlightColor=NULL,
 doHighlightLegend=TRUE,
 ablineH=c(-2, 0, 2),
 base_size=12,
 panel.grid.major.colour="grey90",
 panel.grid.minor.colour="grey95",
 axis.text.x.angle=90,
 return_type=c("ggplot",
    "data"),
 xlim=NULL,
 ylim=c(-6, 6),
 ncol=NULL,
 nrow=NULL,
 blankPlotPos=NULL,
 verbose=FALSE,
 ...)
{
   # deal with missing colnames
   if (!is.list(x) && length(dim(x)) >= 2) {
      if (length(colnames(x)) == 0) {
         colnames(x) <- paste0("column",
            as.character(seq_len(ncol(x))));
      }
      # deal with missing rownames
      if (length(rownames(x)) == 0) {
         rownames(x) <- paste0(#"row",
            as.character(seq_len(nrow(x))));
      }
   }

   # expand titleBoxColor as needed
   if (length(titleBoxColor) == 0) {
      titleBoxColor <- "lightgoldenrod1";
   }
   if (length(names(titleBoxColor)) == 0) {
      if (is.list(x)) {
         titleBoxColor <- jamba::nameVector(
            rep(titleBoxColor,
               length.out=length(x)),
            names(x));
      } else {
         titleBoxColor <- jamba::nameVector(
            rep(titleBoxColor,
               length.out=ncol(x)),
            colnames(x));
      }
   }
   # jamba::printDebugI(titleBoxColor);# debug
   naControlAction <- match.arg(naControlAction);

   # newer method of calling jammacalc()
   if (is.list(x) &&
         "matrix" %in% class(x[[1]]) &&
         all(c("x", "y") %in% colnames(x[[1]]))) {
      # input is previous jammaplot()
      jp2 <- x;
      rm(x);
   } else {
      if (inherits(x, "SummarizedExperiment")) {
         x <- get_se_assaydata(x,
            assay_name=assay_name,
            verbose=verbose);
      }
      ## Optionally apply log2(1 + x) transform?
      if (!TRUE %in% useRank &&
            any(x > apply_transform_limit)) {
         if (verbose) {
            jamba::printDebug("ggjammaplot(): ",
               "Applied transform due to values above ",
               apply_transform_limit, ": ",
               "log2(1 + x)");
         }
         x <- jamba::log2signed(x,
            offset=1);
      }
      # call jammacalc()
      jp2 <- jammacalc(x,
         na.rm=TRUE,
         useMedian=useMedian,
         controlSamples=controlSamples,
         centerGroups=centerGroups,
         controlFloor=controlFloor,
         naControlAction=naControlAction,
         naControlFloor=naControlFloor,
         groupedX=groupedX,
         grouped_mad=grouped_mad,
         whichSamples=whichSamples,
         noise_floor=noise_floor,
         noise_floor_value=noise_floor_value,
         naValue=naValue,
         centerFunc=centerFunc,
         returnType="ma_list",
         mad_row_min=mad_row_min,
         useRank=useRank,
         verbose=verbose,
         ...);
      whichSamples <- names(jp2);
   }

   # jamba::printDebug("whichSamples:");print(whichSamples);# debug
   titleBoxColor <- titleBoxColor[whichSamples];
   # jamba::printDebug("titleBoxColor:");print(titleBoxColor);# debug

   if (length(whichSamples) == 0) {
      whichSamples <- names(jp2);
   }
   if (any(c("integer", "numeric") %in% class(whichSamples))) {
      whichSamples <- whichSamples[whichSamples %in% seq_along(jp2)];
      whichSamples <- names(jp2)[whichSamples];
   }

   return_type <- match.arg(return_type);
   # quick pivot
   # jamba::printDebug("head(jp2[[1]], 10):");print(head(jp2[[1]], 10));# debug
   jp2tall <- jamba::rbindList(
      lapply(unname(whichSamples), function(iname){
         idf <- jp2[[iname]];
         data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            name=rep(iname, length.out=nrow(idf)),
            item=rownames(idf),
            idf)
      })
   )
   rownames(jp2tall) <- NULL;
   head(jp2tall)
   jp2tall <- jamba::renameColumn(jp2tall,
      from=c("x", "y"),
      to=c("mean", "difference"));

   # optional subset of samples
   if (length(whichSamples) > 0) {
      jp2tall <- subset(jp2tall,
         name %in% whichSamples);
   } else {
      whichSamples <- names(jp2);
   }

   # optional blank plot panels
   if (length(blankPlotPos) > 0) {
      sample_order <- character(0);
      total_panels <- max(c(
         length(whichSamples) + length(blankPlotPos),
         blankPlotPos));
      sample_order <- character(total_panels);
      blank_names <- jamba::makeNames(
         rep("blank", length(blankPlotPos)));
      sample_order[blankPlotPos] <- blank_names;
      sample_order[sample_order %in% ""] <- whichSamples;
      # order facet names in same order as jp2
      if (verbose) {
         jamba::printDebug("total_panels: ", total_panels);
         jamba::printDebug("sample_order: ", sample_order);
      }
      jp2tall$name <- factor(jp2tall$name,
         levels=sample_order);
   } else {
      # order facet names in same order as jp2
      jp2tall$name <- factor(jp2tall$name,
         levels=whichSamples);
   }

   # outliers
   mad_outliers <- NULL;
   if (any(!is.na(attr(jp2, "MADfactors")) &
         attr(jp2, "MADfactors") > outlierMAD)) {
      mad_outliers <- names(which(attr(jp2, "MADfactors") > outlierMAD));
   }
   jp2tall$outlier <- ifelse(jp2tall$name %in% mad_outliers,
      TRUE,
      FALSE);

   # subset for non-NA values
   jp2tall <- subset(jp2tall, !is.na(mean) & !is.na(difference))

   # titleBoxColor
   if (length(unique(titleBoxColor)) <= 1) {
      strip_bg <- unique(titleBoxColor);
   } else {
      strip_bg <- "transparent";
   }

   # color gradient
   if (length(colramp) > 0) {
      smooth_colors <- colramp;
   } else {
      smooth_colors <- c("transparent",
         "lightblue",
         "blue",
         "navy",
         "orange",
         "orangered2");
   }
   smooth_colors1 <- jamba::getColorRamp(smooth_colors, n=51);
   baseColor <- head(smooth_colors1, 1);
   if (length(outlierColor) == 1 && all(jamba::isColor(outlierColor))) {
      smooth_colors2 <- c(outlierColor, #"#FF000000", #"#EEE8AA7F",
         tail(jamba::getColorRamp(smooth_colors, n=51), -1));
   } else {
      smooth_colors2 <- jamba::getColorRamp(outlierColor, n=51);
      outlierColor <- head(smooth_colors2, 1);
   }

   # expanded x- and y-axis ranges
   if (length(xlim) == 0) {
      xlims <- lapply(jp2, function(i){
         range(i[,"x"], na.rm=TRUE);
      })
      xlim <- round(digits=2,
         range(unlist(xlims), na.rm=TRUE));
   }
   x_exp <- xlim + diff(range(xlim)) * c(-0.6, 0.6);

   # adjust ylim for useRank=TRUE
   if ((all(c(-6, 6) %in% ylim) || length(ylim) == 0) && useRank) {
      if (verbose > 1) {
         jamba::printDebug("ggjammaplot(): ",
            "over-riding ylim for useRank=TRUE, using ",
            "0.30 * nrow(x)");
      }
      ylim <- c(-1,1) * round(nrow(jp2[[1]]) * 0.30);
   }
   # if more than 40% of all points are outside ylim, then adjust
   # jamba::printDebug("head(jp2tall, 10):");print(head(jp2tall, 10));# debug
   if (length(ylim) == 2 &&
         (sum(abs(jp2tall[, "difference"]) > max(ylim)) / nrow(jp2tall)) > 0.4) {
      if (verbose > 1) {
         jamba::printDebug("ggjammaplot(): ",
            "over-riding ylim because >40% values were outside range.");
         jamba::printDebug("ggjammaplot(): ",
            jamba::formatInt(sum(abs(jp2tall[,"difference"]) > max(ylim))),
            " values out of ",
            jamba::formatInt(nrow(jp2tall)));
      }
      ylims <- jamba::rmNA(unlist(lapply(jp2, function(i){
         range(i[,"y"], na.rm=TRUE);
      })));
      ylim <- max(na.rm=TRUE, c(
         median(ylims[ylims >= 0]),
         median(abs(ylims)[ylims < 0]))) * c(-1, 1);
   }
   if (length(ylim) == 0 || 0 %in% ylim) {
      ylim <- c(-6, 6);
   }
   ylim <- round(digits=2,
      ylim);
   y_exp <- ylim + diff(range(ylim)) * c(-0.6, 0.6);
   if (verbose) {
      jamba::printDebug("ggjammaplot(): ",
         "xlim: ", sapply(xlim, function(i){format(digits=3, i)}),
         ", ylim: ", sapply(ylim, function(i){format(digits=3, i)}));
      jamba::printDebug("ggjammaplot(): ",
         "x_exp: ", format(digits=3, x_exp),
         ", y_exp: ", format(digits=3, y_exp));
      jamba::printDebug("ggjammaplot(): ",
         smooth_colors,
         fgText=list("darkorange1", smooth_colors));
      jamba::printDebug("ggjammaplot(): ",
         smooth_colors1,
         fgText=list("darkorange1", smooth_colors1));
   }

   # ablineH
   if (useRank && all(c(-2, 0, 2) %in% ablineH)) {
      ablineH <- sort(c(0,
         c(-1,1) * round(nrow(jp2[[1]]) * 0.05)));
   }

   # h values for MASS::kde2()
   bw_factor <- rep(bw_factor,
      length.out=2);
   hx <- diff(range(xlim, na.rm=TRUE)) / (20 * bw_factor[1] * detail_factor);
   hy <- diff(range(ylim, na.rm=TRUE)) / (20 * bw_factor[2] * detail_factor);
   # hx <- max(c(hx, hy));
   # hy <- hx;

   ## updated default values
   # nbin <- detail_factor * nbin_factor * 400 / sqrt(length(jp2)) * 2;
   nbin <- detail_factor * (nbin_factor / 3) * 400 / sqrt(length(jp2)) * 2;
   if (verbose > 1) {
      jamba::printDebug("ggjammaplot(): ",
         "kde2d() argument hx: ", format(hx, digits=2),
         ", hy: ", format(hy, digits=2),
         ", n: ", format(nbin[1], digits=2));
   }

   # ggplot MA-plot
   nrow_x <- nrow(x);
   p <- ggplot2::ggplot(
      jp2tall,
      ggplot2::aes(x=mean,
         y=difference)) +
      ggplot2::lims(
         x=xlim,
         y=ylim) +
      ggplot2::xlab(summary) +
      ggplot2::ylab(difference);
      # ggplot2::scale_y_continuous(name=summary)
   if (fillBackground) {
      p <- p +
         ggplot2::geom_rect(
            fill=baseColor,
            data=subset(jp2tall, !outlier & !duplicated(name)),
            alpha=1,
            xmin=x_exp[1],
            xmax=x_exp[2],
            ymin=y_exp[1],
            ymax=y_exp[2])
   }
   p <- p +
      ggplot2::stat_density_2d(
         geom="raster",
         # stat="density_2d",
         data=subset(jp2tall, !outlier),
         ggplot2::aes(fill=(ggplot2::after_stat(count) / nrow_x)^transFactor),
         n=nbin,
         show.legend=FALSE,
         h=c(hx, hy),
         contour=FALSE) +
      # ggplot2::facet_wrap(~name,
      #    drop=FALSE,
      #    nrow=nrow,
      #    ncol=ncol) +
      colorjam::theme_jam(
         strip.background.fill=strip_bg,
         strip.text.size=ggplot2::rel(0.6 * titleCex),
         panel.grid.major.colour=panel.grid.major.colour,
         panel.grid.minor.colour=panel.grid.minor.colour,
         axis.text.x.angle=axis.text.x.angle,
         base_size=base_size);

   # attempt to detect the gradient density range used thus far
   if (length(mad_outliers) > 0) {
      layer_data_2 <- ggplot2::layer_data(p, 2);
      gradient_range_1 <- (range(layer_data_2$count, na.rm=TRUE) / nrow_x) ^ transFactor;
      if (verbose > 1) {
         jamba::printDebug("range(gradient_range_1): ",
            gradient_range_1);
      }
      p <- p +
         ggplot2::scale_fill_gradientn(colours=smooth_colors1,
            limits=gradient_range_1);
   } else {
      p <- p +
         ggplot2::scale_fill_gradientn(colours=smooth_colors1);
   }

   # optionally mark outliers separately
   if (length(mad_outliers) > 0) {
      if (verbose) {
         jamba::printDebug("ggjammaplot(): ",
            "plotting mad_outliers: ",
            mad_outliers);
      }
      if (fillBackground) {
         p <- p +
            ggplot2::geom_rect(
               fill=outlierColor,#"lemonchiffon",
               data=subset(jp2tall, outlier & !duplicated(name)),
               alpha=1,
               xmin=x_exp[1],
               xmax=x_exp[2],
               ymin=y_exp[1],
               ymax=y_exp[2])
      }
      p <- p +
         ggnewscale::new_scale_fill() +
         ggplot2::stat_density_2d(
            geom="raster",
            # stat="density_2d",
            data=subset(jp2tall, outlier),
            ggplot2::aes(fill=(ggplot2::after_stat(count) / nrow_x)^transFactor),
            n=nbin,
            show.legend=FALSE,
            h=c(hx, hy),
            contour=FALSE) +
         ggplot2::scale_fill_gradientn(colours=smooth_colors2,
            limits=gradient_range_1);
   }

   # optionally display MAD factors
   if (displayMAD && length(attr(jp2, "MADfactors")) > 0) {
      mad_factors <- attr(jp2, "MADfactors");
      mad_factors <- mad_factors[names(mad_factors) %in% names(jp2)];
      mad_factors_f <- factor(names(mad_factors),
         levels=levels(jp2tall$name));
      mad_df <- data.frame(name=mad_factors_f,
         mad_factors=mad_factors,
         label=paste0("MAD x",
            format(digits=2,
               mad_factors)),
         x=Inf,
         y=-Inf
      );
      alt_baseColor <- jamba::setTextContrastColor(baseColor);
      alt_outlierColor <- jamba::setTextContrastColor(outlierColor);
      alt_outlierColorL <- jamba::col2hcl(alt_outlierColor)["L",];
      alt_outlierColor <- ifelse(alt_outlierColorL > 50,
         "gold",
         "red3");
      p <- p +
         ggplot2::geom_text(
            data=subset(mad_df, mad_factors < outlierMAD),
            ggplot2::aes(x=x,
               y=y,
               label=label),
            col=alt_baseColor,
            hjust=1.1,
            vjust=-1.0);
      if (length(mad_outliers) > 0) {
         p <- p +
            ggplot2::geom_text(
               data=subset(mad_df, mad_factors >= outlierMAD),
               ggplot2::aes(x=x,
                  y=y,
                  label=label),
               col=alt_outlierColor,
               hjust=1.1,
               vjust=-1.0);
      }
   }

   # optional ablines
   if (length(ablineH) > 0) {
      p <- p +
         ggplot2::geom_hline(yintercept=ablineH,
            linetype="dashed",
            color="#111111",
            alpha=0.5) +
         ggplot2::geom_hline(yintercept=ablineH,
            linetype="dashed",
            color="#FFFFFF",
            alpha=0.5);
   }

   # optional highlight points
   if (length(highlightPoints) > 0) {
      if (verbose) {
         jamba::printDebug("ggjammaplot(): ",
            "plotting ",
            jamba::formatInt(length(unlist(highlightPoints))),
            " highlightPoints.");
      }
      hpl <- handle_highlightPoints(highlightPoints,
         highlightColor=highlightColor,
         highlightPch=highlightPch,
         highlightCex=highlightCex,
         ...);
      highlightPoints <- hpl$highlightPoints;
      highlightColor <- hpl$highlightColor;
      highlightPch <- hpl$highlightPch;
      highlightCex <- hpl$highlightCex;

      highlightPointsTall <- data.frame(
         item=unlist(highlightPoints),
         highlight=rep(names(highlightPoints),
            lengths(highlightPoints)));

      highlightPointsTall2 <- merge(jp2tall,
         highlightPointsTall);
      highlightColorSub <- jamba::nameVector(
         highlightColor,
         names(highlightPoints));
      p <- p +
         ggplot2::geom_point(
            data=highlightPointsTall2,
            ggplot2::aes(x=mean,
               y=difference,
               color=highlight),
            show.legend=TRUE) +
         ggplot2::scale_color_manual(values=highlightColorSub);
   }

   # apply titleBoxColor to strip background and text colors
   facet_colors <- ggh4x::strip_themed(
      background_x=ggh4x::elem_list_rect(fill=titleBoxColor),
      text_x=lapply(titleBoxColor, function(i){
         ggplot2::element_text(
            colour=jamba::setTextContrastColor(
            i))
      }))
   # add to ggplot
   p <- p +
      ggh4x::facet_wrap2(~name,
         drop=FALSE,
         nrow=nrow,
         ncol=ncol,
         strip=facet_colors)

   # optional title
   if (length(maintitle) > 0 || length(subtitle) > 0) {
      p <- p +
         ggplot2::ggtitle(label=maintitle,
            subtitle=subtitle);
   }

   # optionally return ggplot object here
   if (doPlot) {
      if (verbose) {
         jamba::printDebug("ggjammaplot(): ",
            "drawing ggplot2")
      }
      print(p);
      doPlot <- FALSE;
   }
   if ("ggplot" %in% return_type || doPlot) {
      return(invisible(p));
   }

   attr(jp2tall, "MADoutliers") <- mad_outliers;
   return(invisible(jp2tall));
}


#' Get SummarizedExperiment assay matrix data
#'
#' @family jam utility functions
#'
#' @param x `SummarizedExperiment` object
#' @param assay_name `character` string that should match one entry
#'    in `names(assays(x))`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
get_se_assaydata <- function
(x,
 assay_name=NULL,
 verbose=FALSE,
 ...)
{
   if ("SummarizedExperiment" %in% class(x)) {
      if (!jamba::check_pkg_installed("SummarizedExperiment")) {
         stop("The 'SummarizedExperiment' is required for SummarizedExperiment input x.");
      }
      #assay_name <- intersect(assay_name,
      #   names(SummarizedExperiment::assays(x)));
      if (length(assay_name) == 0) {
         assay_name <- head(names(SummarizedExperiment::assays(x)), 1);
         if (verbose) {
            jamba::printDebug("get_se_assaydata(): ",
               c("Using first assay_name: '", assay_name, "'"),
               sep="");
         }
      }
      if (is.numeric(assay_name)) {
         if (assay_name > length(SummarizedExperiment::assays(x))) {
            assay_name <- length(SummarizedExperiment::assays(x));
         }
         assay_name <- names(SummarizedExperiment::assays(x))[assay_name];
      }
      x <- SummarizedExperiment::assays(x)[[assay_name]];
      x_names <- colnames(x);
      nsamples <- length(x_names);
      if (length(x) == 0) {
         stop("assays(x)[[assay_name]] did not produce a usable data matrix.");
      }
   } else {
      if (!(is.matrix(x) && is.numeric(x))) {
         stop("x must be a numeric matrix or SummarizedExperiment.");
      }
   }
   x;
}

#' Handle highlightPoints argument to jammaplot()
#'
#' @family jam utility functions
#'
handle_highlightPoints <- function
(highlightPoints=NULL,
   highlightColor=NULL,
   highlightPch=21,
   highlightCex=1.5,
   ...)
{
   #
   if (length(highlightPoints) == 0) {
      return(NULL)
   }

   # convert to list
   #
   # - if there are enough colors for each vector element, convert each
   # to their own list
   #
   # - otherwise convert vector to one list, color them all the same
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
   if (length(highlightColor) == 0) {
      highlightColor <- jamba::nameVector(
         colorjam::rainbowJam(
            n=length(highlightPoints)),
         names(highlightPoints));
   }
   #if (!is.list(highlightColor)) {
   #   highlightColor <- as.list(highlightColor);
   #}
   if (length(names(highlightColor)) > 0 && all(names(highlightColor) %in% names(highlightPoints))) {
      highlightColor <- highlightColor[names(highlightPoints)];
   }
   if (length(highlightColor) != length(highlightPoints)) {
      highlightColor <- rep(highlightColor,
         length.out=length(highlightPoints));
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
   if (!c("list") %in% class(highlightCex)) {
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
   return(list(
      highlightPoints=highlightPoints,
      highlightColor=highlightColor,
      highlightPch=highlightPch,
      highlightCex=highlightCex
   ))
}
