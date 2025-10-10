# =====================================================
# Imports and linking configuration
# =====================================================

#' @useDynLib bayesQRsurvey, .registration = TRUE
#' @import Rcpp
#'
#' @importFrom stats acf mad median predict qnorm quantile rnorm runif sd var setNames
#' @importFrom stats delete.response model.frame model.matrix model.response terms coef
#' @importFrom graphics plot lines points axis legend arrows grid mtext par segments
#' @importFrom methods setClass setMethod new is
#' @importFrom utils globalVariables
#' @importFrom rlang .data
NULL

# =====================================================
# Global variables for NSE (Non-Standard Evaluation)
# =====================================================

utils::globalVariables(c(
  "tau", "x", "yhat_lo", "yhat_up",
  "xid", "y1", "y2", "y3", "tau_f",
  "med", "lo", "hi"
))
