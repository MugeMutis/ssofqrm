#' Air Quality Data Example
#'
#' A dataset containing air quality metrics used in the spatial quantile regression models.
#'
#' @format A list including three elements.
#' \describe{
#'   \item{data_2023}{Data for the year 2023.}
#'   \item{data_2024}{Data for the year 2024.}
#'   \item{wie_mat}{Spatial weight matrix.}
#' }
#'
#' @source https://cran.r-project.org/web/packages/ARPALData/index.html
#' @usage data(air_data)
#' @examples
#' \dontrun{
#' data(air_data)
#' yy <- matrix(air_data$data_2023$PM2.5_mean, ncol = 365, byrow = T)
#' y <- apply(yy, 1, mean)
#' x <- matrix(air_data$data_2023$Ozone_max_8h, ncol = 365, byrow = T)
#' yy_test <- matrix(air_data$data_2024$PM2.5_mean, ncol = 366, byrow = T)[,1:365]
#' y_test <- apply(yy_test, 1, mean)
#' x_test <- matrix(air_data$data_2024$Ozone_max_8h, ncol = 366, byrow = T)[,1:365]
#' wei_mat <- air_data$wei_mat
#' fit_kim <- fqrm(y=y, x=x, w=wei_mat, tau=0.5, method = "KM")
#' fit_ch <- fqrm(y=y, x=x, w=wei_mat, tau=0.5, method = "Ch")
#' predict_kim <- predict_fqrm(object = fit_kim, xnew = x_test, wnew = wei_mat)
#' predict_ch <- predict_fqrm(object = fit_ch, xnew = x_test, wnew = wei_mat)
#' }
#'
"air_data"
