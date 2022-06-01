
#' Convert log2 fold change to signed fold change
#'
#' Convert log2 fold change to signed fold change
#'
#' This function takes log2 fold change values as input, and returns
#' normal space fold change values that retain the positive and negative
#' sign, and the magnitude.
#'
#' For example:
#'
#' * `log2 fold change = 2` becomes `fold change = 4`.
#' * `log2 fold change = -2` becomes `fold change = -4`.
#'
#' This function therefore differs from similar functions that convert
#' log2 fold change into a ratio. Instead, `log2fold_to_fold()`
#' specifically retains the magnitude of negative changes.
#'
#' @param x `numeric` vector
#'
#' @family jam utility functions
#'
#' @return `numeric` vector representing signed fold change values.
#'
#' @examples
#' x <- c(-3, -2, -1, 0, 1, 2, 3);
#' fc <- log2fold_to_fold(x);
#' fc;
#'
#' fold_to_log2fold(fc);
#'
#' @export
log2fold_to_fold <- function
(x,
 ...)
{
   x_sign <- sign(x);
   x_sign <- ifelse(x_sign == 0, 1, x_sign);
   2^(abs(x)) * x_sign;
}

#' Convert normal signed fold change to log2 fold change
#'
#' Convert normal signed fold change to log2 fold change
#'
#' This function takes fold change values as input, and returns
#' log2 fold change values.
#'
#' This function recognizes two forms of input:
#'
#' * ratio, which includes values between 0 and 1, but no negative values;
#' * fold change, as from `log2fold_to_fold()` which includes no values
#' between 0 and 1, but may include negative values.
#'
#' For example, for ratio input:
#'
#' * `ratio = 4` becomes `log2 fold change = 2`.
#' * `ratio = 0.25` becomes `log2 fold change = -2`.
#'
#' For example, for fold change input:
#'
#' * `fold change = 4` becomes `log2 fold change = 2`.
#' * `fold change = -4` becomes `log2 fold change = -2`.
#'
#' @param x `numeric` vector
#'
#' @return `numeric` vector representing log2 fold change values.
#'
#' @family jam utility functions
#'
#' @examples
#' x <- c(-3, -2, -1, 0, 1, 2, 3);
#' fc <- log2fold_to_fold(x);
#' fc;
#'
#' fold_to_log2fold(fc);
#'
#' @export
fold_to_log2fold <- function
(x,
 ...)
{
   is_ratio <- (x > 0 & x < 1);
   ifelse(is_ratio,
      log2(x),
      log2(abs(x)) * sign(x));
}
