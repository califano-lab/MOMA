#' Empirical Cumulative Distribution Helper function
#' 
#' This is a wrapper function for the ecdf function in order to get a vector of 
#' the resulting values
#' @importFrom stats ecdf
#' @importFrom tidyr replace_na
#' @param vals original list of p-values to convert
#' @param na_value value to substitute for NAs if present. Default is NA
#' @return vector of adjusted p-values
#' @keywords internal
cdf.pval <- function(vals, na_value = NA) {
  # could tweak this...
  vals.adj <- tidyr::replace_na(vals, na_value)
  fn <- stats::ecdf(tidyr::replace_na(vals, 1))
  res <- fn(vals.adj)
  res
}