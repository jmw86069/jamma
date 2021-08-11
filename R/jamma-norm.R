#' Normalize data using MA-plot logic
#'
#' Normalize data using MA-plot logic
#'
#' This function normalizes data using an approach analogous
#' to data viewed in MA-plots. Normalization is applied by shifting
#' data along the y-axis so the mean (or median) expression among
#' control genes is zero, indicated on an MA-plot as y=0.
#'
#' **Note:** This method should be performed only after reviewing the MA-plots,
#' to ensure the assumptions are met. Similarly, data can also be
#' viewed in MA-plots after normalization to confirm and review the
#' effect of normalization.
#'
#' It is useful to run `jammaplot()` data after `jammanorm()`
#' to visualize the effect of this normalization. If data is not
#' centered at y=0, the parameters should be adjusted.
#'
#' ## Assumptions
#'
#' This method effectively reinforces the assumption that the mean log fold
#' change for control genes is expected to be zero. When `useMedian=TRUE`
#' it reinforces the assumption that the log fold change for the
#' majority of control genes should be zero.
#'
#' Therefore, the assumptions may be summarized as follows:
#'
#' * The principle assumption is that the set of `controlGenes`,
#' whose mean expression is at or above the `minimum_mean` value,
#' are unchanged within the respective `centerGroups`. For typical
#' whole genome transcript microarray and RNA-seq experiments,
#' this assumption is typically valid when using `useMedian=TRUE`.
#' For experiments with specific reference genes, or housekeeper
#' genes, this assumption may only be true for those specific genes.
#' * The data signal is assumed to be a roughly linear representation
#' of the relative abundance of each measured entity, which is usually
#' true for log-transformed microarray intensity, or log-transformed
#' RNA-seq read counts. For QPCR or TaqMan, this assumption is valid
#' for the direct CT cycle threshold values, or after log-transformation
#' of the exponentiated CT values, for example 2^(40-CT).
#' All that said, a straightforward way to visualize this assumption
#' is with MA-plots, to confirm that signal is horizontal across
#' the full signal range - either for the majority of all genes,
#' or for the specific `controlGenes` used for normalization.
#' * The variability among control genes should not be more than twice
#' the median absolute deviation across other samples within the
#' relevant `centerGroups`. Effectively this assumption means the
#' control genes on the MA-plot should not show wide spread along the
#' y-axis.
#'
#' In cases where some samples show non-horizontal signal across the
#' MA-plot, the data is not conforming to a consistent and proportional
#' signal across the response range of the experiment. In effect, it
#' means signal is being compressed, or expanded along the response
#' as compared to other samples in the same `centerGroups`. In this
#' scenario, the best normalization method may be
#' `limma::normalizeQuantiles()`, `limma::normalizeCyclicLoess()`,
#' or `vsn::vsn()`. These methods adjust the distribution of signal
#' to enforce consistency across samples.
#'
#' In general, the signal distribution itself should not be adjusted
#' unless necessary, in order to retain as much information from the
#' underlying technology as possible. This method `jammanorm()` is
#' intended to apply linear normalization, which effectively shifts
#' the entire signal for a sample up or down relative to other samples.
#'
#' This scenario is effective for technologies such as QPCR, TaqMan,
#' Nanostring counts, and RNA-seq counts or RNA-seq pseudocounts.
#'
#' When the MA-plot demonstrates non-horizontal
#' signal, it is most often the result one or both of these influences:
#'
#' 1. batch effect, imposed either by different upstream sample
#' processing steps among the samples being tested, or
#' 2. platform technology that tends to produce relative signal strength
#' but not absolute quantitative signal, commonly seen with
#' microarray hybridization technologies such as Affymetrix, Illumina,
#' Agilent, SomaLogic, or Myriad RBM.
#'
#' Note that any upstream sample amplification technique may also impose
#' non-linear effects on the molecules being measured.
#'
#' One method to test for a batch effect is to define `centerGroups`
#' to include batch, so the data will be centered for each batch
#' independently. If this centering resolves the non-horizontal
#' signal, then batch is very likely to be a component to be modeled
#' in the experiment. See `limma::removeBatchEffect()`. The batch effect
#' adjustment by  `limma::removeBatchEffect()` and `sva::ComBat()`
#' almost exactly subtract the batch component from the signal.
#'
#' That said, it may or may not be ideal to apply batch adjustment
#' prior to running downstream statistical tests, as opposed to
#' including batch as a covariate term in the statistical
#' model used for testing, example when using `DESeq2`.
#' The main benefit of applying batch adjustment at this step is to
#' visualize data downstream consistent with the method used by
#' those statistical tests, or when running a clustering technique
#' that does not have the capability of applying appropriate
#' batch effect modeling.
#'
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
#' @return `numeric` `matrix` whose columns have been normalized,
#' with the following named `attributes`:
#'    * `"nf"` numeric vector of normalization factors;
#'    * `"hk"` a `character` vector of controlGenes for each sample;
#'    * `"hk_count"` the `integer` number of controlGenes for each sample.
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
 useMedian=TRUE,
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
         jamba::printDebug("jammanorm(): ",
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
      useMedian=useMedian,
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
         if (useMedian) {
            nf <- median(y, na.rm=TRUE);
         } else {
            nf <- mean(y, na.rm=TRUE);
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
