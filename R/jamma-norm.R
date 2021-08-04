#' Normalize data using MA-plot logic
#'
#' Normalize data using MA-plot logic
#'
#' This function normalizes data using an approach analogous
#' to data viewed in MA-plots. Normalization is applied by shifting
#' data along the y-axis so the mean (or median) expression among
#' control genes is zero, shown on an MA-plot as y=0.
#'
#' It effectively reinforces the assumption that the mean log fold
#' change for control genes is expected to be zero. When `useMedian=TRUE`
#' it reinforces the assumption that the log fold change for the
#' majority of control genes should be zero.
#'
#' It is useful to run `jammaplot()` data after `jammanorm()`
#' to visualize the effect of this normalization. If data is not
#' centered at y=0, the parameters should be adjusted.
#'
#' ## About the normalization
#'
#' This normalization is a "linear normalization" in that it uniformly
#' shifts data up or down relative to other samples, without
#' modifying the relative distribution of signal. It is very similar
#' to housekeeper normalization, geometric mean normalization, and
#' global signal scaling, which are all also "linear normalization"
#' methods. An example of non-linear normalization is quantile or
#' VSN normalization.
#'
#' The control genes can be defined upfront with `controlGenes`,
#' which can be housekeeper genes, or a custom subset of genes.
#' The `controlGenes` are filtered to require the mean or median
#' expression at or above `minimum_mean` in order to be used
#' during normalization.
#'
#' The `minimum_mean` threshold is useful and important to match
#' the variability seen in the MA-plots. For example data below
#' a certain x-axis value may have very high variability, and
#' should usually not be used for normalization.
#'
#' When the MA-plot after normalization does not show signal centered
#' at y=0, the most common and effective adjustment is to apply
#' `minimum_mean` to require `controGenes` to have expression
#' at or above this threshold. The next most effective option is
#' `useMedian=TRUE` which will center the majority of genes at
#' y=0 instead of the overall mean at y=0.
#'
#' ## Comparison to geometric mean normalization
#'
#' The end result is very similar to other housekeeper normalization
#' methods which typically define normalization factors by
#' calculating the geometric mean of log-transformed housekeeper
#' gene abundances. Such approaches usually work in part because
#' higher expressed housekeeper genes usually have lower variability,
#' which keeps the influence on geometric mean reasonably consistent
#' across a broad range of expression. That said, genes with higher
#' expression have more influence on the geometric mean than genes
#' with much lower expression.
#'
#' However, the strategy with `jammanorm()` is to assert that
#' the mean difference from average expression for the `controlGenes`
#' should be equal to zero. The effect is applied evenly across
#' control genes by evaluating the mean or median difference from y=0
#' for each sample.
#'
#' ## Assumptions
#'
#' In order for this approach to be valid:
#'
#' * Data on the MA-plots should be horizontal for all samples,
#' particularly for `controlGenes`. When data is not horizontal
#' across samples, data should instead be normalized using another
#' approach, typically something like quantile normalization which
#' is intended to impose a consistent signal distribution across
#' all samples.
#' * The majority of `controlGenes` should not be changing, therefore
#' the `controlGenes` should not have substantial variation along
#' the y-axis within each sample.
#' * For data with a small number of genes, normalization should use
#' a set of pre-defined `controlGenes` that were verified to have
#' no change across samples.
#'
#' ## Noise threshold
#'
#' Note that some platform technologies generate a noise threshold,
#' below which data may be skewed up or down relative to other samples.
#' It is recommended to ignore this type of skew below the threshold
#' when determining whether data is horizontal on MA-plots.
#'
#' For example Nanostring data includes a series of positive and
#' negative control probes, and a suitable noise threshold is either
#' the midpoint between the lowest positive and highest negative probe,
#' or the lowest positive probe. When this noise threshold is applied,
#' data above the noise threshold is typically horizontal, although
#' data below the threshold may be skewed up or down depending upon
#' the effective input RNA concentration.
#'
#'
#' @return numeric matrix whose columns have been normalized using
#' the y-axis mean offset following MA-plot calculations. The output
#' contains attributes: `"nf"` numeric vector of normalization factors;
#' `"hk"` list of controlGenes used for each sample; `"hk_count"` the
#' number of controlGenes used for each sample.
#'
#' @inheritParams jammaplot
#' @param x `numeric` matrix with expression data suitable for use
#'    by `jammaplot()`. Gene expression data is typically transformed
#'    using `log2(1+x)` to represent reasonably normal distribution.
#' @param minimum_mean `numeric` value used to filter `controlGenes`,
#'    a control gene must have at least this level of expression to
#'    be included in the normalization, where the expression is
#'    determined by the mean or median value analogous to the x-axis
#'    value in MA-plots.
#' @param noise_floor,noise_floor_value `numeric` values passed
#'    to `jammacalc()`. The `noise_floor` is the value below which
#'    a floor is applied. The floor sets all values below this floor
#'    to `noise_floor_value`. For example, one could apply `noise_floor=0`
#'    and `noise_floor_value=NA` which would change any value below 0
#'    to `NA`.
#'
#' @family jam matrix functions
#'
#' @export
jammanorm <- function
(x,
 controlGenes=NULL,
 minimum_mean=0,
 controlSamples=NULL,
 centerGroups=NULL,
 useMedian=FALSE,
 useMean=NULL,
 noise_floor=-Inf,
 noise_floor_value=noise_floor,
 verbose=FALSE,
 ...)
{
   ## Purpose is to use MA-plot logic to provide a normalization
   ## method designed to result in data centered at y=0.
   ##
   ## It optionally takes control_genes, such as housekeeper genes,
   ## to use as controls for normalization.
   ##
   ## It optionally filters control_genes to ensure the mean abundance
   ## is at least minimum_mean.
   ##
   ## This function essentially takes the output from jammaplot()
   ## and applies the inverse of the y-axis mean per sample.
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

   jpr <- jammacalc(x,
      controlSamples=controlSamples,
      centerGroups=centerGroups,
      useMean=useMean,
      noise_floor=noise_floor,
      noise_floor_value=noise_floor_value,
      returnType="ma_list");

   if (length(minimum_mean) == 0) {
      minimum_mean <- 0;
   }
   if (length(controlGenes) == 0) {
      controlGenes <- rownames(x);
   }
   if (verbose) {
      jamba::printDebug("jammanorm(): ",
         "length(controlGenes):",
         length(controlGenes));
   }

   ## Calculate HK genes and normalization factors
   jpr_hknf <- lapply(jpr, function(i){
      j <- as.data.frame(i);
      hk <- intersect(controlGenes,
         rownames(subset(j, x >= minimum_mean)));
      if (length(hk) == 0) {
         nf <- NA;
      } else {
         y <- j[hk,"y"];
         if (useMean) {
            nf <- mean(y, na.rm=TRUE);
         } else {
            nf <- median(y, na.rm=TRUE);
         }
      }
      list(nf=nf, hk=hk);
   });
   jpr_nf <- unlist(lapply(jpr_hknf, function(i){
      i$nf;
   }));
   jpr_hk <- lapply(jpr_hknf, function(i){
      i$hk;
   });

   ## Now adjust the input data using the nf values
   x_hk <- t(t(x) - jpr_nf);

   ## attributes to store NF and HK
   attr(x_hk, "hk") <- jpr_hk;
   attr(x_hk, "hk_count") <- lengths(jpr_hk);
   attr(x_hk, "nf") <- (- jpr_nf);

   x_hk;
}
