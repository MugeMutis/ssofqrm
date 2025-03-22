#'  Kim and Muller (2004) Estimator for Spatial Quantile Regression
#'
#' Estimates spatial quantile regression parameters using the two-stage approach proposed by Kim et al.
#'
#' @param y Numeric vector. Response variable of length \code{n}.
#' @param x Numeric matrix of size \code{n x p}. Predictor variables.
#' @param w Numeric matrix of size \code{n x n}. Spatial weight matrix.
#' @param tau Numeric. Quantile level (between 0 and 1).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{b0}{Estimated intercept.}
#'   \item{rho}{Estimated spatial autoregressive parameter.}
#'   \item{b}{Estimated regression coefficients (excluding intercept and \code{rho}).}
#' }
#'
#' @importFrom quantreg rq
#' @importFrom stats fitted
kim_fun <- function(y, x, w, tau){

  nvar <- dim(x)[2]
  wy <- w %*% y
  wx <- w %*% x
  w2x <- w %*% wx

  fit1 <- rq(wy ~ x + wx + w2x, tau = tau)
  wyhat <- fitted(fit1)
  fit2 <- rq(y ~ x + wyhat, tau = tau)
  coefs <- fit2$coefficients
  b0 <- coefs[1]
  rho <- coefs[length(coefs)]
  b <- coefs[c(-1, -(nvar + 2))]

  return(list(b0 = b0, rho = rho, b = b))
}
