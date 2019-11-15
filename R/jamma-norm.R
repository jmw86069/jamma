#' Normalize data using MA-plot logic
#'
#' Normalize data using MA-plot logic
#'
#' This function normalizes data using the MA-plot y-axis
#' offset to define the normalization factor. It can optionally
#' use only a subset of `controlGenes`, such as known
#' housekeeper genes. It can optionally apply a `minimum_mean`
#' threshold, to ensure the mean abundance of each row
#' is at least this value.
#'
#' The approach is very similar to other housekeeper normalization
#' methods, which typically use the geometric mean (or the mean
#' of the log-transformed abundances.) However, the strategy here
#' is to assert that the mean difference from average for the
#' selected `controlGenes` is equal to zero.
#'
#' In order for this approach to be valid:
#'
#' * Data on the MA-plots should be horizontal for all samples.
#' * The majority of controlGenes should not be changing, which
#' is seen as the controlGenes representing the mean expression.
#' * For data with a small number of genes, known controlGenes
#' should be selected and verified to be unchanging across samples.
#'
#' @return numeric matrix whose columns have been normalized using
#' the y-axis mean offset following MA-plot calculations. The output
#' contains attributes: `"nf"` numeric vector of normalization factors;
#' `"hk"` list of controlGenes used for each sample; `"hk_count"` the
#' number of controlGenes used for each sample.
#'
#' @param x numeric matrix with expression data suitable for use
#'    by `jammaplot()`. Gene expression data is typically transformed
#'    using `log2(1+x)` to roughly normally distributed data.
#'
#'
#' @export
jammanorm <- function
(x,
 controlGenes=NULL,
 minimum_mean=0,
 controlSamples=NULL,
 centerGroups=NULL,
 useMean=TRUE,
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
