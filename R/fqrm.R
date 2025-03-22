#' Functional Quantile Regression Model with Spatial Autocorrelation
#'
#' Fits a spatial scalar-on-function quantile regression model using either the Kim and Muller (2004)
#' or Chernozhukov and Hansen (2006) method.
#'
#' @param y Numeric vector of length \code{n}. Response variable.
#' @param x Numeric matrix of size \code{n x p}. Functional predictor observed on a common grid.
#' @param w Numeric matrix of size \code{n x n}. Spatial weight matrix.
#' @param tau Numeric. Quantile level to estimate (between 0 and 1).
#' @param argvals Numeric vector of length \code{p}. Grid of argument values for the functional predictor. If \code{NULL}, defaults to a regular grid on [0, 1].
#' @param nbasis Integer. Number of basis functions to use in FPCA. If \code{NULL}, defaults to \code{min(10, p / 4)}.
#' @param method Character string. Estimation method to use. One of \code{"KM"} (Kim et al.) or \code{"Ch"} (Chernozhukov et al.).
#'
#' @return A list containing:
#' \describe{
#'   \item{b}{Estimated coefficient vector in FPCA space.}
#'   \item{b0}{Estimated intercept.}
#'   \item{bhat}{Estimated coefficient function in the original function space.}
#'   \item{rho}{Estimated spatial autoregressive parameter.}
#'   \item{fitted.values}{Fitted values from the model.}
#'   \item{residuals}{Residuals from the model.}
#'   \item{tau}{Quantile level used.}
#'   \item{fpca}{Output from \code{getPCA}, including FPCA scores and basis functions.}
#' }
#'
#' @importFrom stats as.formula
#' @examples
#' \dontrun{
#' sim_data <- dgp(n = 250, rho=0.5)
#' y <- sim_data$y
#' x <- sim_data$x
#' w <- sim_data$w
#' fit_kim <- fqrm(y=y, x=x, w=w, tau=tau, method = "KM")
#' fit_ch <- fqrm(y=y, x=x, w=w, tau=tau, method = "Ch")
#' }
#'
#' @export
fqrm <- function(y, x, w, tau, argvals = NULL, nbasis = NULL,
                 method = c("KM","Ch")){

  w <- norm_wei(w)
  method <- match.arg(method)
  n <- dim(x)[1]
  p <- dim(x)[2]

  if(is.null(argvals))
    argvals <- seq(0, 1, length.out = p)
  if(is.null(nbasis))
    nbasis <- min(10, p/4)

  fpc_per <- getPCA(x, nbasis, argvals)
  fsco <- fpc_per$PCAscore
  efuns <- fpc_per$evalbase %*% fpc_per$PCAcoef$coefs

  if(method == "KM")
    qmodel <- kim_fun(y, fsco, w, tau)
  if(method == "Ch")
    qmodel <- chernozhukov_fun(y, fsco, w, tau)

  b0 <- qmodel$b0
  rho <- qmodel$rho
  b <- qmodel$b

  bhat <- efuns %*% b
  fits <- solve(diag(n) - rho * w) %*% as.matrix(rep(1, n)) * b0 +
    solve(diag(n) - rho * w) %*% fsco %*% b
  resids <- y - fits

  return(list(b = b, b0 = b0,
              bhat = bhat,
              rho = rho,
              fitted.values = fits,
              residuals = resids,
              tau = tau,
              fpca = fpc_per))
}
