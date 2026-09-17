#' Bayesian quantile regression models for complex survey data
#'
#' The bayesQRsurvey package provides Bayesian quantile regression methods for complex
#' survey designs through two main functions.
#'
#' \itemize{
#'   \item \code{bqr.svy()} estimates single-output quantile regression models
#'    using MCMC methods.
#'   \item \code{mo.bqr.svy()} estimates multiple-output quantile regression
#'    models using an EM algorithm.
#' }
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{bqr.svy}}}{Fits Bayesian quantile regression for multiple quantiles using MCMC methods (ALD, Score, Approximate)}
#'   \item{\code{\link{mo.bqr.svy}}}{Fits Bayesian quantile regression for multiple quantiles using EM algorithm}
#'   \item{\code{\link{prior}}}{Unified interface for creating prior distributions}
#' }
#'
#' @section Methods for fitted models:
#' \describe{
#'   \item{\code{\link[=summary.bqr.svy]{summary}}}{Posterior summaries for \code{bqr.svy} and \code{mo.bqr.svy} objects}
#'   \item{\code{\link[=plot.bqr.svy]{plot}}}{Plot methods for \code{bqr.svy} and \code{mo.bqr.svy} objects}
#'   \item{\code{\link[=predict.bqr.svy]{predict}}}{Predictions at new covariate values}
#'   \item{\code{\link{diagnostics}}}{MCMC and EM convergence diagnostics}
#'   \item{\code{\link[=bqr.svy.posterior]{posterior_interval}}}{Posterior credible intervals}
#'   \item{\code{\link[=bqr.svy.methods]{bqr.svy methods}}}{Generic accessors for \code{bqr.svy} fits: \code{coef}, \code{vcov}, \code{fitted}, \code{sigma}, \code{nobs}, \code{formula}, \code{terms}, \code{model.matrix}, \code{weights}, \code{as.matrix} and \code{update}}
#'   \item{\code{\link[=mo.bqr.svy.methods]{mo.bqr.svy methods}}}{Generic accessors for \code{mo.bqr.svy} fits: \code{coef}, \code{vcov}, \code{sigma}, \code{nobs}, \code{formula}, \code{terms}, \code{model.matrix}, \code{weights} and \code{update}}
#' }
#'
#' @section MCMC Methods:
#' The bqr.svy function can estimate three types of models, where the quantile regression
#' coefficients are defined at the super-population level, and their estimators are
#' built upon the survey weights.
#' \itemize{
#'   \item \strong{ALD (Asymmetric Laplace Distribution)}, which uses the
#'    asymmetric Laplace likelihood.
#'   \item \strong{Score}, which uses a score-based approach.
#'   \item \strong{Approximate}, which uses approximate methods for faster
#'    computation.
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
