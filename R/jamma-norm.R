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
#' @returns `numeric` `matrix` whose columns have been normalized,
#'    with the following `attributes`:
#'    * `'nf'`: `numeric` normalization factors, which
#'    were applied through addition `+` to original signal.
#'    A value `NA` indicates that no `controlGenes` remained
#'    after filtering by `minimum_mean` and `maximum_mean`.
#'    * `'hk'`: `character` vector of `controlGenes` used for
#'    each sample, after applying `minimum_mean` and `maximum_mean`
#'    thresholds.
#'    * `'hk_count'`: `integer` number of `controlGenes` for each
#'    sample.
#'    * `'params'`: `list` with other argument values provided.
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
#'    
#'    Some notes:
#'    * For RNA-seq data, and other similar count-based 'Omics data,
#'    a common default is `minimum_mean=5`. The value 5 in log2 space
#'    corresponds to approximately 32 read counts (or pseudocounts)
#'    in normal space counts.
#'    For normalization, this default is effective at
#'    avoiding "shot noise" from low integer counts.
#'    * For NanoString data, we often use the lowest positive
#'    control probe mean signal, or the point at which the positive
#'    controls are still horizontal, in the event that the lowest
#'    positive control probe becomes non-linear across the dilution
#'    range.
#'    * For hybridization-based platforms, it may be useful to review
#'    MA-plots for the point at which variability becomes high, or
#'    the point at which signal is no longer horizontal when compared
#'    with other samples.
#'    In the event the MA-plots are generally not consistently horizontal,
#'    other normalization options are recommended, such as
#'    quantile or LOESS normalization.
#' @param maximum_mean `numeric`, default `Inf` (infinite), to define
#'    the highest mean value for a row to be used in normalization.
#'    The `maximum_mean` is useful to avoid using rows with extremely
#'    high signal, which may have aberrant values that may therefore
#'    skew the normalization.
#' 
#'    Together with `minimum_mean` one can define the exact minimum
#'    and maximum mean values for control genes to be used.
#' @param noise_floor `numeric` threshold below which values are
#'    replaced with `noise_floor_value`. The default is `-Inf`
#'    (negative infinite) which applies no floor.
#'    Values are passed internally to `jammacalc()`.
#'    * A common option is `noise_floor=0` which replaces any value
#'    below `0` with the `noise_floor_value` whose default is
#'    the `noise_floor`. It therefore replaces negative values
#'    with `0`.
#'    * Another common option is to define a platform noise floor,
#'    a value below which the platform measurement is no longer
#'    trustworthy.  
#'    Values below the platform floor are replaced with the floor,
#'    which has the effect of making any reported measurement
#'    at or below the floor equivalent to the floor itself.
#'    * For QPCR, in some circumstances, a cycle threshold
#'    (Ct) value above 35 maybe considered outside of trusted
#'    range. When using virtual abundance calculation  
#'    `abundance = 2^(40 - Ct)`  
#'    it would suggest `noise_floor=5`, equivalent to setting
#'    any Ct higher than 35 to 35. It retains the measurement,
#'    but does not inflate a potential fold change calculation.
#' @param noise_floor_value `numeric` default is `noise_floor`
#'    which replaces values below the `noise_floor` with the
#'    `noise_floor` itself.
#'    * A common alternative is to use `noise_floor=0`, and
#'    `noise_floor_value=NA`, which replaces all values below `0`
#'    with `NA`.
#'
#' @family jam matrix functions
#'
#' @export
jammanorm <- function(
   x,
   controlGenes=NULL,
   minimum_mean=0,
   maximum_mean=Inf,
   controlSamples=NULL,
   centerGroups=NULL,
   useMedian=FALSE,
   useMean=NULL,
   noise_floor=-Inf,
   noise_floor_value=noise_floor,
   verbose=FALSE,
   ...
) {
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
   # use only one value overall
   useMedian <- head(useMedian, 1);

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
   if (length(maximum_mean) == 0) {
      maximum_mean <- Inf;
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
      hk <- intersect(
         controlGenes,
         rownames(
            subset(
               j,
               !is.na(x) &
               x >= minimum_mean &
               x <= maximum_mean
            )
         )
      )
      if (length(hk) == 0) {
         # no rows remain to normalize data,
         # return NA
         nf <- NA;
      } else {
         y <- j[hk, "y"];
         if (isTRUE(useMedian)) {
            nf <- median(y, na.rm=TRUE);
         } else {
            nf <- mean(y, na.rm=TRUE);
         }
      }
      list(
         nf=nf,
         hk=hk
      );
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
   attr(x_hk, "nf") <- (- jpr_nf);
   attr(x_hk, "hk") <- jpr_hk;
   attr(x_hk, "hk_count") <- lengths(jpr_hk);
   params <- list(
      minimum_mean=minimum_mean,
      maximum_mean=maximum_mean,
      useMedian=useMedian,
      centerGroups=centerGroups,
      controlSamples=controlSamples,
      noise_floor=noise_floor,
      noise_floor_value=noise_floor_value
   )
   attr(x_hk, "params") <- params;

   x_hk;
}
