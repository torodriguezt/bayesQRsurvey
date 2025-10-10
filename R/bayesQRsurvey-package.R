#' bayesQRsurvey: Bayesian Weighted Quantile Regression for complex survey designs
#' with EM and MCMC Algorithm
#'
#' The bayesQRsurvey package provides Bayesian quantile regression methods for complex
#' survey designs with two main functions:
#'
#' \itemize{
#'   \item \code{bqr.svy()}: Bayesian methods for estimating quantile regression models
#'    using MCMC methods
#'   \item \code{mo.bqr.svy()}: Bayesian approach to multiple-output quantile regression
#'    using EM algorithm
#' }
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{bqr.svy}}}{Fits Bayesian quantile regression for multiple quantiles using MCMC methods (ALD, Score, Approximate)}
#'   \item{\code{\link{mo.bqr.svy}}}{Fits Bayesian quantile regression for multiple quantiles using EM algorithm}
#'   \item{\code{\link{plot}}}{Standard plot method for bqr.svy objects}
#'   \item{\code{\link{summary}}}{Unified summary method for all bayesQRsurvey model objects}
#'   \item{\code{\link{prior}}}{Unified interface for creating prior distributions}
#' }
#'
#' @section MCMC Methods:
#' The bqr.svy function can estimate three types of models, where the quantile regression
#' coefficients are defined at the super-population level, and their estimators are
#' built upon the survey weights.:
#' \itemize{
#'   \item \strong{ALD (Asymmetric Laplace Distribution)}: Uses asymmetric Laplace likelihood
#'   \item \strong{Score}: Uses score-based approach
#'   \item \strong{Approximate}: Uses approximate methods for faster computation
#' }
#'
#' @section EM Algorithm:
#' Implements a Bayesian approach to multiple-output quantile regression for complex
#' survey data analysis.
#'
#' @references
#' Yu, K. and Moyeed, R. A. (2001). Bayesian quantile regression.
#' \emph{Statistics & Probability Letters}, 54(4), 437-447.
#'
#' Kozumi, H. and Kobayashi, G. (2011). Gibbs sampling methods for Bayesian
#' quantile regression. \emph{Journal of Statistical Computation and Simulation},
#' 81(11), 1565-1578.
#'
#' @author Marcus L. Nascimento, Kelly Cristina Mota Goncalves,
#'         Johnatan Cardona Jimenez, Tomas Rodriguez Taborda
#'
#' @keywords internal
#' @useDynLib bayesQRsurvey, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom grDevices adjustcolor
#' @importFrom graphics axis grid legend lines par points segments
"_PACKAGE"
