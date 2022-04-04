
#' MA-plots using ggplot2
#'
#' MA-plots using ggplot2
#'
#' This method is under active development and may change as
#' features are implemented.
#'
#' It is currently fully functional and is being documented.
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
#'       farrisdata::colorSub[as.character(colData(GeneSE)$groupName)],
#'       colnames(GeneSE));
#'    options("warn"=FALSE);
#'
#'    gg <- ggjammaplot(GeneSE,
#'       ncol=6,
#'       base_size=12,
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
#' @param x one of the following inputs:
#'    * `numeric` matrix
#'    * `SummarizedExperiment` object, where the
#'    assay data is defined using `assays(x)[[assay_name]]`. Accordingly,
#'    `assay_name` can be either an integer index or `character` string
#'    matching the `names(assays(x))`.
#'    * `list` output from `jammacalc()` or `jammaplot()`, where each
#'    element in the list is a two-column `matrix` with colnames `c("x", "y")`.
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
#' @param noise_floor `numeric` value indicating the minimum numeric value,
#'    below which the value is converted to `noise_floor_value`. Typically
#'    this threshold is useful to convert values `0` into `NA` so they will
#'    not be used in the MA-plots nor the associated calculations. It is
#'    also useful to establish a noise floor above zero, for example
#'    in QPCR if noise threshold for calculated expression is roughly 50,
#'    a useful noise floor might be `noise_floor=log2(50)` and
#'    `noise_floor_value=log2(50)`, which will set any expression value
#'    at or below `log2(50)` to that value. It prevents noise from becoming
#'    meaningful, when it may be just instrument noise with no relevance.
#'    Those values could also reasonably be converted to `NA` for this
#'    purpose.
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
   nbin_factor=1,
   bw_factor=1,
   assay_name=1,
   useMedian=FALSE,
   controlSamples=NULL,
   centerGroups=NULL,
   groupedX=TRUE,
   grouped_mad=TRUE,
   outlierMAD=5,
   mad_row_min=4,
   displayMAD=FALSE,
   noise_floor=1e-10,
   noise_floor_value=NA,
   naValue=NA,
   centerFunc=centerGeneData,
   whichSamples=NULL,
   useRank=FALSE,
   titleBoxColor="lightgoldenrod1",
   outlierColor="palegoldenrod",
   maintitle=NULL,
   subtitle=NULL,
   summary="mean",
   difference="difference",
   transFactor=0.2,
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
   # expand titleBoxColor as needed
   if (length(names(titleBoxColor)) == 0) {
      if (is.list(x)) {
         titleBoxColor <- jamba::nameVector(
            rep(titleBoxColor,
               length.out=length(x)),
            names(x));
      } else {
         titleBoxColor <- jamba::nameVector(
            rep(titleBoxColor,
               length.out=length(colnames(x))),
            colnames(x));
      }
   }

   # newer method of calling jammacalc()
   if (is.list(x) &&
         "matrix" %in% class(x[[1]]) &&
         all(c("x", "y") %in% colnames(x[[1]]))) {
      # input is previous jammaplot()
      jp2 <- x;
      rm(x);
   } else {
      if ("SummarizedExperiment" %in% class(x)) {
         x <- get_se_assaydata(x,
            assay_name=assay_name,
            verbose=verbose);
      }

      jp2 <- jammacalc(x,
         na.rm=TRUE,
         useMedian=useMedian,
         controlSamples=controlSamples,
         centerGroups=centerGroups,
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

   titleBoxColor <- titleBoxColor[whichSamples];

   if (length(whichSamples) == 0) {
      whichSamples <- names(jp2);
   }
   if (any(c("integer", "numeric") %in% class(whichSamples))) {
      whichSamples <- whichSamples[whichSamples %in% seq_along(jp2)];
      whichSamples <- names(jp2)[whichSamples];
   }

   return_type <- match.arg(return_type);
   # quick pivot
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
   mad_outliers <- attr(jp2, "MADoutliers");
   if (!any(mad_outliers %in% names(jp2))) {
      mad_outliers <- NULL;
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
   smooth_colors <- c("transparent",
      "lightblue",
      "blue",
      "navy",
      "orange",
      "orangered2");
   smooth_colors1 <- jamba::getColorRamp(smooth_colors, n=51);
   smooth_colors2 <- c("#FF000000", #"#EEE8AA7F",
      tail(jamba::getColorRamp(smooth_colors, n=51), -1));

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
   if (length(ylim) == 2 &&
         (sum(abs(jp2tall[,"difference"]) > max(ylim)) / nrow(jp2tall)) > 0.4) {
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
   hx <- diff(range(xlim, na.rm=TRUE)) / 30 * bw_factor[1];
   hy <- diff(range(ylim, na.rm=TRUE)) / 30 * bw_factor[2];

   nbin <- nbin_factor * 400 / sqrt(length(jp2));
   if (verbose > 1) {
      jamba::printDebug("ggjammaplot(): ",
         "kde2d() argument hx: ", format(hx, digits=2),
         ", hy: ", format(hy, digits=2),
         ", n: ", format(nbin[1], digits=2));
   }

   # ggplot MA-plot
   p <- ggplot2::ggplot(
      jp2tall,
      ggplot2::aes(x=mean,
         y=difference)) +
      ggplot2::lims(
         x=xlim,
         y=ylim) +
      ggplot2::xlab(summary) +
      ggplot2::xlab(difference) +
      # ggplot2::scale_y_continuous(name=summary)
      ggplot2::stat_density_2d(geom="raster",
         stat="density_2d",
         data=subset(jp2tall, !outlier),
         ggplot2::aes(fill=(..density..)^transFactor),
         n=nbin,
         show.legend=FALSE,
         h=c(hx, hy),
         contour=FALSE) +
      ggplot2::scale_fill_gradientn(colours=smooth_colors1) +
      ggplot2::facet_wrap(~name,
         drop=FALSE,
         nrow=nrow,
         ncol=ncol) +
      colorjam::theme_jam(strip.background.fill=strip_bg,
         panel.grid.major.colour=panel.grid.major.colour,
         panel.grid.minor.colour=panel.grid.minor.colour,
         base_size=base_size);

   # optionally mark outliers separately
   if (length(mad_outliers) > 0) {
      if (verbose) {
         jamba::printDebug("ggjammaplot(): ",
            "plotting mad_outliers: ",
            mad_outliers);
      }
      p <- p +
         ggplot2::geom_rect(
            fill="lemonchiffon",
            data=subset(jp2tall, outlier & !duplicated(name)),
            alpha=0.7,
            xmin=x_exp[1],
            xmax=x_exp[2],
            ymin=y_exp[1],
            ymax=y_exp[2]) +
         ggnewscale::new_scale_fill() +
         ggplot2::stat_density_2d(geom="raster",
            stat="density_2d",
            data=subset(jp2tall, outlier),
            ggplot2::aes(fill=(..density..)^transFactor),
            n=nbin,
            show.legend=FALSE,
            h=c(hx, hy),
            contour=FALSE) +
         ggplot2::scale_fill_gradientn(colours=smooth_colors2)
   }

   # optionally display MAD factors
   if (displayMAD && length(attr(jp2, "MADfactors")) > 0) {
      mad_factors <- attr(jp2, "MADfactors");
      mad_factors <- mad_factors[names(mad_factors) %in% names(jp2)];
      mad_factors_f <- factor(names(mad_factors),
         levels=levels(jp2tall$name));
      mad_df <- data.frame(name=mad_factors_f,
         label=paste0("MAD x",
            format(digits=2,
               mad_factors)),
         x=Inf,
         y=-Inf
      );
      p <- p +
         ggplot2::geom_text(data=mad_df,
            ggplot2::aes(x=x,
               y=y,
               label=label),
            hjust=1.1,
            vjust=-1.0);
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
      print(head(highlightPointsTall2));
      highlightColorSub <- jamba::nameVector(
         highlightColor,
         names(highlightPoints));
      print(highlightColorSub);
      jamba::printDebug(highlightColorSub);
      p <- p +
         ggplot2::geom_point(
            data=highlightPointsTall2,
            ggplot2::aes(x=mean,
               y=difference,
               color=highlight),
            show.legend=TRUE) +
         ggplot2::scale_color_manual(values=highlightColorSub);
   }

   # manually adjust titleBoxColor with strip text background color
   if (length(unique(titleBoxColor)) > 1 &&
         jamba::check_pkg_installed("ggtext")) {
      p <- p +
         ggplot2::theme(
            strip.background=ggplot2::element_blank(),
            strip.text=element_textbox_colorsub(
               colorSub=titleBoxColor,
               box.colour="black"
            )
         );
   }

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
   return(list(
      highlightPoints=highlightPoints,
      highlightColor=highlightColor,
      highlightPch=highlightPch,
      highlightCex=highlightCex
   ))
}
