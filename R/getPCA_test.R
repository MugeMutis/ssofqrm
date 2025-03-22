#' Project Test Data onto FPCA Basis
#'
#' Computes principal component scores for new functional test data using an existing FPCA model.
#'
#' @param object A list object returned from \code{\link{getPCA}} containing FPCA model components (e.g., basis, eigenfunctions, mean function, grid points).
#' @param data Numeric matrix of size \code{n x p}. Test data observed on the same grid used in training.
#'
#' @return A matrix of principal component scores for the test data.
#'
#' @import fda
getPCA_test <- function(object, data)
{
  bs_basis <- object$bs_basis
  PCAcoef <- object$PCAcoef
  gp <- object$gp
  mean.tr <- c(object$meanScore$coefs)

  n <- dim(data)[1]
  p <- dim(data)[2]
  dimnames(data) <- list(as.character(1:n), as.character(1:p))

  pcaobj <- smooth.basisPar(gp, t(data), bs_basis, Lfdobj = NULL, lambda = 0)$fd
  sdata <- scale(t(pcaobj$coefs), center = mean.tr, scale = FALSE)

  pcaobj2 <- pcaobj
  pcaobj2$coefs <- t(sdata)

  PCAscore_test <- inprod(pcaobj2, PCAcoef)
  colnames(PCAscore_test) <- paste0("Score", 1:dim(PCAcoef$coefs)[2])

  return(PCAscore_test)
}
