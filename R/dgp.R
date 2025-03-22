#' Data Generating Process for Spatial Scalar-on-Function Model
#'
#' Generates synthetic data from a spatial scalar-on-function regression model.
#'
#' @param n Integer. Sample size.
#' @param rho Numeric. Spatial autoregressive parameter.
#' @param out.p Numeric. Proportion of outliers to introduce (default is 0).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{y}{Response vector of length \code{n}.}
#'   \item{x}{Matrix of functional predictor values.}
#'   \item{w}{Spatial weight matrix of size \code{n x n}.}
#'   \item{tcoefs}{True coefficient function used to generate the data.}
#'   \item{out.index}{Indices of outlier observations, if any.}
#' }
#'
#' @import fda.usc
#' @importFrom stats rnorm
#' @examples
#' \dontrun{
#' sim_data <- dgp(n = 250, rho=0.5)
#' y <- sim_data$y
#' x <- sim_data$x
#' w <- sim_data$w
#' tcoef <- sim_data$tcoefs
#' }
#'
#' @export
dgp <- function(n, rho, out.p = 0){

  X <- rproc2fdata(n = n, t = seq(0, 1, l = 101), sigma = "OU")

  beta0 <- fdata(mdata = cos(2*pi*seq(0,1,l=101))-(seq(0,1,l=101)-0.5)^2,
                 argvals = seq(0,1,l=101), rangeval=c(0,1))

  wei <- matrix(, n, n)

  for(i in 1:n){
    for(j in 1:n){
      if(i != j){
        wei[i,j] <- 1/(abs(i-j))
      }else{
        wei[i,j] <- 0
      }
    }
  }

  W <- matrix(0, n, n)
  for(i in 1:n)
    W[i,] <- wei[i,] / sum(wei[i,])

  err <- rnorm(n, sd = 0.1)
  argx <- inprod.fdata(X, beta0)

  out.index <- NULL
  if(out.p > 0){
    nout <- round(n * out.p)
    out.index <- sample(1:n, nout)
    err[out.index] <- rnorm(nout, mean = 5, sd = 0.1)

    fYe = solve(diag(n) - rho * W) %*% argx + solve(diag(n) - rho * W) %*% err
  }

  fYe = solve(diag(n) - rho * W) %*% argx + solve(diag(n) - rho * W) %*% err

  return(list("y" = fYe, "x" = X$data, w = W, tcoefs = beta0$data, out.index = out.index))
}
