#' Chernozhukov Estimator for Spatial Quantile Regression
#'
#' Implements the Chernozhukov and Hansen (2006) estimation method for spatial quantile regression
#' with instrumental variable correction and a grid search over possible spatial autocorrelation values.
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
#'   \item{b}{Estimated regression coefficients (excluding intercept).}
#' }
#'
#' @importFrom quantreg rq
#' @importFrom stats fitted
chernozhukov_fun <- function(y, x, w, tau){

  nvar <- dim(x)[2]
  rho_candidates <- seq(-0.95, 0.95, 0.001)
  nrho <- length(rho_candidates)

  wy <- w %*% y
  wx <- w %*% x
  w2x <- w %*% wx
  fit.init <- rq(wy ~ x + wx + w2x, tau = 0.5)
  wyhat <- fitted(fit.init)

  rho_hat <- numeric()
  for(i in 1:nrho) {
    newy <- y - rho_candidates[i] * wy
    fit <- rq(newy ~ x + wyhat, tau = tau)
    rho_hat[i] <- fit$coef[length(fit$coef)]
  }

  rho_min <- rho_candidates[which.min(abs(rho_hat))]

  newy2 <- y - rho_min * wy
  fit.final <- rq(newy2 ~ x, tau = tau)
  coefs <- fit.final$coef

  b0 <- coefs[1]
  rho <- rho_min
  b <- coefs[-1]

  return(list(b0 = b0, rho = rho, b = b))
}
