
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
   if (inherits(x, "SummarizedExperiment")) {
      if (!requireNamespace("SummarizedExperiment", quietly=TRUE)) {
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

#' Get SummarizedExperiment column  data from colData
#'
#' @family jam utility functions
#'
#' @param x `SummarizedExperiment` object
#' @param use_values `character` values, default NULL. Optionally used
#'    to validate input values in one of two ways:
#'    1. Vector of values typically named by `colnames(x)` that should
#'    be applied to each column in `x` matching the order of the colnames.
#'    2. Vector of values that matches `colnames(colData(x))`, which
#'    should be converted to named vector of values, with names matching
#'    `colnames(x)`.
#'
#'    Either way, when `use_values` is supplied, the output will be
#'    a named vector, with names defined using `colnames(x)`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @returns `data.frame` or `DFrame` when `use_values` is not supplied,
#'    or `character` vector named by `colnames(x)` when `use_values`
#'    is supplied.
#'
#' @examplesIf {requireNamespace("SummarizedExperiment", quietly=TRUE)}
#' set.seed(123)
#' m <- matrix(rnorm(8100), ncol=9)
#' colnames(m) <- head(letters, 9)
#' rownames(m) <- as.character(1:900)
#' countse <- SummarizedExperiment::SummarizedExperiment(
#'    assays=list(counts=m))
#' SummarizedExperiment::colData(countse)$type <- rep(LETTERS[1:3], each=3)
#' SummarizedExperiment::colData(countse)$class <- rep(c("WT", "KO"), c(6, 3))
#'
#' # use colname in colData
#' get_se_colData(countse, "type")
#'
#' # use two colnames in colData
#' get_se_colData(countse, c("class", "type"))
#'
#' # use three colnames in colData, one is missing
#' get_se_colData(countse, c("class", "type", "detail"))
#'
#' # use named vector, which put it into proper order colnames(x)
#' use_values <- setNames(tail(LETTERS, 9), rev(head(letters, 9)))
#' use_values
#' get_se_colData(countse, use_values)
#'
#' # use unnamed vector, which does not re-order values
#' unname(use_values)
#' get_se_colData(countse, unname(use_values))
#' @export
get_se_colData <- function
(x,
 use_values=NULL,
 verbose=FALSE,
 ...)
{
   if (!inherits(x, "SummarizedExperiment")) {
      stop("SummarizedExperiment input is required by get_se_colData()");
   }
   if (!requireNamespace("SummarizedExperiment", quietly=TRUE)) {
      stop("The 'SummarizedExperiment' is required for SummarizedExperiment input x.");
   }
   se_colData <- SummarizedExperiment::colData(x)
   # handle use_values when supplied
   if (length(use_values) > 0) {
      if (all(colnames(x) %in% names(use_values))) {
         use_values <- use_values[colnames(x)];
      } else {
         if (any(use_values %in% colnames(se_colData))) {
            centerGroups_cols <- intersect(use_values,
               colnames(se_colData));
            use_values <- jamba::nameVector(
               jamba::pasteByRow(
                  se_colData[, centerGroups_cols, drop=FALSE]),
               colnames(x))
         } else if (length(use_values) == ncol(x)) {
            if (length(names(use_values)) == 0) {
               names(use_values) <- colnames(x);
            } else {
               stop("colnames(x) must be present in names(use_values)");
            }
         } else {
            stop(paste0("names(use_values) must match colnames(x), or",
               "use_values values must be in colnames(colData(x))"))
         }
      }
      return(use_values);
   }
   return(se_colData);
}
