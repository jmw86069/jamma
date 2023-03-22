
#' Convert numeric matrix to column-rank data
#'
#' Convert numeric matrix to column-rank data
#'
#' This function converts a numeric matrix to rank-based data,
#' where the higher rank value refers to higher numeric value in
#' the input data. It therefore behaves similar to continuous
#' numeric data:
#'
#' * higher values indicate higher signal
#' * difference in values indicates difference in rank
#'
#' @family jam matrix functions
#'
#' @return `numeric` matrix with same dimensions as input data.
#'
#' @param x `numeric` matrix data
#' @param ... additional arguments are ignored.
#'
#' @export
matrix_to_column_rank <- function
(x,
 keepNA=TRUE,
 usePercentile=FALSE,
 ...)
{
   #
   if (length(x) == 0) {
      return(x)
   }
   if (TRUE %in% keepNA && any(is.na(x))) {
      x_is_na <- is.na(x);
   } else {
      x_is_na <- FALSE;
   }
   x1 <- apply(x, 2, function(xi){
      xr <- rank(xi,
         na.last=FALSE);
      if (TRUE %in% keepNA & any(is.na(xi))) {
         xr[is.na(xi)] <- NA;
      }
      if (TRUE %in% usePercentile) {
         xrmin <- min(xr, na.rm=TRUE);
         xrmax <- max(xr, na.rm=TRUE);
         if (xrmax == xrmin) {
            xrmin <- 0
         }
         xr <- 100 * (xr - xrmin) /
            (xrmax - xrmin);
      }
      xr;
   })
   x[] <- x1;
   if (TRUE %in% keepNA && any(x_is_na)) {
      x[x_is_na] <- NA;
   }
   return(x)
}
