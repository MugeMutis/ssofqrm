#' Row-Normalize a Spatial Weight Matrix
#'
#' Checks whether a spatial weight matrix is row-normalized. If not, it normalizes the matrix by dividing each row by its sum.
#'
#' @param W Numeric matrix. Spatial weight matrix of size \code{n x n}.
#' @param tol Numeric. Tolerance for checking whether rows sum to 1. Default is \code{1e-8}.
#'
#' @return A row-normalized version of the input matrix \code{W}.

norm_wei <- function(W, tol = 1e-8) {
  row_sums <- rowSums(W)

  if (all(abs(row_sums - 1) < tol)) {
    message("The weight matrix is already row-normalized.")
    return(W)
  } else {
    message("The weight matrix is not row-normalized. Normalizing now.")
    W_normalized <- sweep(W, 1, row_sums, FUN = "/")
    return(W_normalized)
  }
}
