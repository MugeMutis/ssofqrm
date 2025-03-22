#' Predict Method for Functional Quantile Regression with Spatial Effects
#'
#' Makes predictions for new observations using a fitted spatial functional quantile regression model.
#'
#' @param object A fitted model object returned by \code{\link{fqrm}}.
#' @param xnew Numeric matrix of size \code{n x p}. New functional predictor data (on the same grid used in training).
#' @param wnew Numeric matrix of size \code{n x n}. Spatial weight matrix for the new observations.
#'
#' @return A numeric vector of predicted values of length \code{n}.
#' @examples
#' \dontrun{
#' sim_data <- dgp(n = 250, rho=0.5)
#' y <- sim_data$y
#' x <- sim_data$x
#' w <- sim_data$w
#' fit_kim <- fqrm(y=y, x=x, w=w, tau=tau, method = "KM")
#' fit_ch <- fqrm(y=y, x=x, w=w, tau=tau, method = "Ch")
#' sim_test <- dgp(n = 1000, rho=rho)
#' y_test <- sim_test$y
#' x_test <- sim_test$x
#' w_test <- sim_test$w
#' predict_kim <- predict_fqrm(object = fit_kim, xnew = x_test, wnew = w_test)
#' predict_ch <- predict_fqrm(object = fit_ch, xnew = x_test, wnew = w_test)
#' }
#'
#' @export
predict_fqrm <- function(object, xnew, wnew){

  n <- dim(xnew)[1]

  sco_new <- getPCA_test(object$fpca, xnew)

  b0 <- object$b0
  b <- object$b
  rho <- object$rho

  preds <- solve(diag(n) - rho * wnew) %*% as.matrix(rep(1, n)) * b0 +
    solve(diag(n) - rho * wnew) %*% sco_new %*% b

  return(preds)
}
