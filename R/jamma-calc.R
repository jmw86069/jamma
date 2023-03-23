
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
#' @param x `numeric` matrix typically containing log-normal measurements,
#'    with measurement rows, and sample columns.
#' @param na.rm `logical` indicating whether to ignore NA values
#'    during `numeric` summary functions.
#' @param controlSamples `character` vector containing values in
#'    `colnames(x)` to define control samples used during centering.
#'    These values are passed to `centerGeneData()`.
#' @param centerGroups `character` vector with length equal to `ncol(x)`
#'    which defines the group for each column in `x`. Data will
#'    be centered within each group.
#' @param groupedX `logical` indicating how to calculate the x-axis
#'    value when `centerGroups` contains multiple groups. When
#'    groupedX=TRUE, the mean of each group median is used, which
#'    has the effect of representing each group equally. When
#'    groupedX=FALSE, the median across all columns is used, which
#'    can have the effect of preferring sample groups with a larger
#'    number of columns.
#' @param useMedian `logical` indicating whether to use the median
#'    values when calculating the x-axis and during data centering.
#'    The median naturally reduces the effect of outlier points on
#'    the resulting MA-plots., when compared to using the mean.
#'    When useMedian=FALSE, the mean value is used.
#' @param useMean (deprecated) `logical` indicating whether to use the
#'    mean instead of the median value. This argument is being removed
#'    in order to improve consistency with other Jam package functions.
#' @param whichSamples `character` vector containing `colnames(x)`, or
#'    integer vector referencing column numbers in `x`. This argument
#'    specifies which columns to return, but does not change the columns
#'    used to define the group centering values. For example, the
#'    group medians are calculated using all the data, but only the
#'    samples in `whichSamples` are centered to produce MA-plot data.
#' @param noise_floor `numeric` value indicating the minimum numeric value
#'    allowed in the input matrix `x`. When `NULL` or `-Inf` no noise
#'    floor is applied. It is common to set `noise_floor=0` to limit
#'    MA-plot data to use values zero and above.
#' @param noise_floor_value single `numeric` value used to replace `numeric`
#'    values at or below `noise_floor` when `noise_floor` is not NULL.
#'    By default,
#'    `noise_floor_value=noise_floor` which means values at or below
#'    the noise floor are set to the floor. Another useful option is
#'    `noise_floor_value=NA` which has the effect of removing the point
#'    from the MA-plot altogether. This option is recommended for sparse
#'    data matrices where the presence of values at or below zero are
#'    indicative of missing data (zero-inflated data) and does not
#'    automatically reflect an actual value of zero.
#' @param naValue single `numeric` value used to replace any `NA` values in
#'    the input matrix `x`. This argument can be useful to replace
#'    `NA` values with something like zero.
#' @param mad_row_min `numeric` value defining the minimum group
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
#' @param grouped_mad `logical` indicating whether the MAD value
#'    should be calculated per group when `centerGroups` is
#'    supplied, from which the MAD factor values are derived.
#'    When `TRUE` it has the effect of highlighting outliers
#'    within each group using the variability in that group.
#'    When `FALSE` the overall MAD is calculated, and a
#'    particularly high variability group may have all its
#'    group members labeled with a high MAD factor.
#' @param centerFunc `function` used for centering data, by default
#'    one of the functions `centerGeneData()` or `centerGeneData_v1()`.
#'    This argument will be removed in the near future and is mainly
#'    intended to allow testing the two centering functions.
#'    The following arguments are passed to this function:
#'    * x: the input `numeric` data matrix
#'    * na.rm: `logical` whether to ignore NA value. Always use `na.rm=TRUE`.
#'    * controlSamples: `character` optional subset of `colnames(x)` to
#'    use as reference controls during centering
#'    * centerGroups: `character` vector of groups for `colnames(x)`
#'    * controlFloor: `numeric` optional minimum allowed value for control
#'    summary prior to centering
#'    * naControlAction: `character` string for how to handle entirely NA
#'    control groups during centering
#'    * naControlFloor: `numeric` used when `naControlAction="floor"` and
#'    all control values are `NA`. One `numeric` value is inserted into
#'    the control group.
#'    * useMedian: `logical` whether to use median (TRUE) or mean (FALSE)
#'    * returnGroups: `logical` whether to return summary of group assignment
#'    in attribute `"center_df"`
#'    * returnGroupedValues: `logical` whether to return group summary values
#'    in attribute `"x_group"`
#'    * ...: other arguments are passed along via `...`.
#' @param returnType `character` string indicating the format of data
#'    to return: `"ma_list"` is a list of MA-plot two-column
#'    numeric matrices with colnames `c("x","y")`; "tidy"
#'    returns a tall `data.frame` suitable for use in ggplot2.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#'
#' @export
jammacalc <- function
(x,
 na.rm=TRUE,
 controlSamples=NULL,
 centerGroups=NULL,
 controlFloor=NA,
 naControlAction=c("row", "floor", "min", "na"),
 naControlFloor=0,
 groupedX=TRUE,
 useMedian=TRUE,
 useMean=NULL,
 whichSamples=NULL,
 noise_floor=-Inf,
 noise_floor_value=noise_floor,
 naValue=NA,
 mad_row_min=0,
 grouped_mad=TRUE,
 centerFunc=centerGeneData,
 useRank=FALSE,
 returnType=c("ma_list", "tidy"),
 verbose=FALSE,
 ...)
{
   ## Purpose is to separate jammaplot() from the underlying
   ## math that produces the data for MA-plots.
   if (length(x) == 0 || ncol(x) == 0 || nrow(x) == 0) {
      return(NULL);
   }
   returnType <- match.arg(returnType);
   naControlAction <- match.arg(naControlAction);

   if (length(useMean) > 0 && is.logical(useMean)) {
      useMedian <- !useMean;
      if (verbose) {
         jamba::printDebug("jammacalc(): ",
            "useMedian defined by !useMean, useMedian=",
            useMedian);
      }
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
      whichSamples <- match(whichSamples, colnames(x));
   }
   if (verbose) {
      jamba::printDebug("jammacalc(): ",
         "whichSamples:",
         whichSamples);
   }

   if (!jamba::check_pkg_installed("matrixStats")) {
      rowMedians <- function(x, na.rm=TRUE) {
         apply(x, 1, median, na.rm=na.rm);
      }
      colMedians <- function(x, na.rm=TRUE) {
         apply(x, 2, median, na.rm=na.rm);
      }
   } else {
      rowMedians <- matrixStats::rowMedians;
      colMedians <- matrixStats::colMedians;
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
   if (TRUE %in% useRank) {
      x <- matrix_to_column_rank(x,
         keepNA=TRUE,
         ...)
      if (FALSE) {
         x1 <- apply(x, 2, function(xi){
            rank(
               jamba::rmNA(xi,
                  naValue=min(xi, na.rm=TRUE) - 1),
               na.last=FALSE);
         })
         if (any(is.na(x))) {
            x1[is.na(x)] <- NA;
         }
         x <- x1;
         rm(x1);
      }
   }

   ## Center the data in one step
   x_ctr <- centerFunc(x,
      na.rm=na.rm,
      controlSamples=controlSamples,
      centerGroups=centerGroups,
      controlFloor=controlFloor,
      naControlAction=naControlAction,
      naControlFloor=naControlFloor,
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
      jammadata <- jamba::rbindList(lapply(names(jammadata), function(i){
         j <- jammadata[[i]];
         data.frame(
            check.names=FALSE,
            stringsAsFactors=FALSE,
            item=rownames(j),
            name=i,
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
