
#' Center gene data
#'
#' Performs row centering on a matrix of data in log space
#'
#' This function is deprecated in favor of `centerGeneData()`
#' which includes more flexibility in data centering. It is
#' maintained here for backward compatibility.
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
#' @return `numeric` `matrix` with the same row and column dimensions
#' as input data. If `returnGroupedValues=TRUE`, the additional columns contain
#' the row median or mean values, dependent upon mean=FALSE or mean=TRUE,
#' respectively.
#' An attribute \code{centerGroups} is included, which describes the specific
#' relationship between each colname, and associated control sample colnames,
#' and optional centerGroups grouping of colnames. When columns are grouped
#' and centered to specific control samples, is it important to keep this
#' information during downstream scrutiny of results.
#'
#' @family jamma deprecated functions
#'
#' @importFrom matrixStats rowMedians
#'
#' @examples
#' x <- matrix(1:100, ncol=10);
#' colnames(x) <- letters[1:10];
#' # basic centering
#' centerGeneData_v1(x);
#'
#' # grouped centering
#' centerGeneData_v1(x,
#'    centerGroups=rep(c("A","B"), c(5,5)));
#'
#' # centering versus specific control columns
#' centerGeneData_v1(x,
#'    controlSamples=letters[c(1:3)]);
#'
#' # grouped centering versus specific control columns
#' centerGeneData_v1(x,
#'    centerGroups=rep(c("A","B"), c(5,5)),
#'    controlSamples=letters[c(1:3, 6:8)]);
#'
#' # confirm the centerGroups and controlSamples
#' x_ctr <- centerGeneData_v1(x,
#'    centerGroups=rep(c("A","B"), c(5,5)),
#'    controlSamples=letters[c(1:3, 6:8)],
#'    showGroups=TRUE);
#'
#' attr(x_ctr, "centerDF");
#'
#' @export
centerGeneData_v1 <- function
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
   if (jamba::check_pkg_installed("matrixStats")) {
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
         jamba::printDebug("centerGeneData_v1(): ",
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
      jamba::printDebug("centerGeneData_v1(): ",
         "dim(indata):",
         dim(indata));
      jamba::printDebug("centerGeneData_v1(): ",
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
      centeredSubsets <- lapply(jamba::nameVectorN(centerGroups), function(iGroupN){
         iGroup <- centerGroups[[iGroupN]];
         iGroupCols <- colnames(indata[,iGroup,drop=FALSE]);
         iControls <- intersect(iGroupCols, controlCols);
         iM <- centerGeneData_v1(indata[,iGroup,drop=FALSE],
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
      if (length(jamba::tcount(colnames(centeredData), minCount=2)) > 0) {
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
         jamba::printDebug("centerGeneData_v1(): ",
            "Keeping non-numeric columns as-is.");
      }
      centeredData <- cbind(indataChar, centeredData)
   }
   attr(centeredData, "centerDF") <- centerDF;
   return(centeredData);
}


#' Center gene data
#'
#' Performs per-row centering on a numeric matrix
#'
#' This function centers data by subtracting the median or
#' mean for each row.
#'
#' Columns can be grouped using argument `centerGroups`.
#' Each group group of columns defined by `centerGroups`
#' is centered independently.
#'
#' Data can be centered relative to specific control columns
#' using argument `controlSamples`.
#' When `controlSamples` is not supplied, the default behavior
#' is to use all columns. This process is consistent with
#' typical MA-plots.
#'
#' It may be preferred to define `controlSamples` in cases where
#' there are known reference samples, against which other samples
#' should be compared.
#'
#' The `controlSamples` logic is applied independently to each
#' group defined in `centerGroups`.
#'
#' You can confirm the `centerGroups` and `controlSamples` are
#' correct in the result data, by accessing the attribute
#' `"center_df"`, see examples below.
#'
#' Note: This function assumes input data is suitable for
#' centering by subtraction.
#' This data requirement is true for:
#'
#' * most log-transformed gene expression data
#' * quantitative PCR (QPCR) cycle threshold (CT) values
#' * other numeric data that has been suitably transformed
#' to meet reasonable parametric assumption of normality,
#' * rank-transformed data which results in difference in rank
#' * generally speaking, any data where the difference between 5 and 7 (2)
#' is reasonably similar to the difference between 15 and 17 (2).
#'
#' The data requirement is not true for:
#'
#'  * most gene expression data in normal space (hint: if any value is above 100)
#'  * numeric data that is strongly skewed
#'  * generally speaking, any data where the difference between 5 and 7
#'  is not reasonably similar to the difference between 15 and 17. If
#'  the percent difference is more likely to be the interesting measure,
#'  data may be log-transformed for analysis.
#'
#'  For special cases, `rowStatsFunc` can be supplied to perform
#'  specific group summary calculations per row.
#'
#' @param x `numeric` matrix of input data. See assumptions,
#'    that data is assumed to be log2-transformed, or otherwise
#'    appropriately transformed.
#' @param centerGroups `character` vector of group names, or
#'    `NULL` if there are no groups.
#' @param na.rm `logical` indicating whether NA values should be
#'    ignored for summary statistics. This argument is passed
#'    to the corresponding row stats function. Frankly, this
#'    value should be `na.rm=TRUE` for all stat functions by default,
#'    for example `mean(..., na.rm=TRUE)` should be default.
#' @param controlSamples `character` vector of values in `colnames(x)`
#'    which defines the columns to use when calculating group
#'    summary values.
#' @param useMedian `logical` indicating whether to use group median
#'    values when calculating summary statistics `TRUE`, or
#'    group means `FALSE`. In either case, when `rowStatsFunc`
#'    is provided, it is used instead.
#' @param rmOutliers `logical` indicating whether to perform outlier
#'    detection and removal prior to row group stats. This
#'    argument is passed to `jamba::rowGroupMeans()`. Note that
#'    outliers are only removed during the row group summary step,
#'    and not in the centered data.
#' @param madFactor `numeric` value passed to `jamba::rowGroupMeans()`,
#'    indicating the MAD factor threshold to use when `rmOutliers=TRUE`.
#'    The MAD of each row group is computed, the overall group median
#'    MAD is used to define 1x MAD factor, and any MAD more than
#'    `madFactor` times the group median MAD is considered an outlier
#'    and is removed. The remaining data is used to compute row
#'    group values.
#' @param rowStatsFunc `optional` function used to calculate row group
#'    summary values. This function should take a numeric matrix as
#'    input, and return a one-column numeric matrix as output, or
#'    a numeric vector with length `nrow(x)`. The function should
#'    also accept `na.rm` as an argument.
#' @param returnGroupedValues `logical` indicating whether to include
#'    the numeric matrix of row group values used during centering,
#'    returned in the attributes with name `"x_group"`.
#' @param returnGroups `logical` indicating whether to return the
#'    centering summary data.frame in attributes with name "center_df".
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `jamba::rowGroupMeans()`.
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
#'    returnGroups=TRUE);
#' attr(x_ctr, "center_df");
#'
#' @export
centerGeneData <- function
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
 mean=NULL,
 verbose=FALSE,
 ...)
{
   ## This function is a refactor of centerGeneData() to consolidate
   ## some logic into rowGroupMeans()
   if (length(x) == 0 || ncol(x) == 0 || nrow(x) == 0) {
      return(x);
   }

   if (length(mean) > 0 && is.logical(mean)) {
      useMedian <- !mean;
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
         jamba::printDebug("centerGeneData(): ",
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
