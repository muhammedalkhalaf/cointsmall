#' @keywords internal
"_PACKAGE"

#' cointsmall: Cointegration Tests with Structural Breaks in Small Samples
#'
#' @description
#' The cointsmall package implements cointegration tests designed for small 
#' sample sizes with potential structural breaks. It follows the methodology
#' of Trinh (2022), which provides small-sample adjusted critical values
#' for testing the null hypothesis of no cointegration.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{cointsmall}}}{Main function to perform cointegration 
#'     tests with optional structural breaks}
#'   \item{\code{\link{cointsmall_combined}}}{Combined testing procedure that 
#'     evaluates all model specifications}
#'   \item{\code{\link{cointsmall_cv}}}{Retrieve critical values for different 
#'     model specifications}
#' }
#'
#' @section Models:
#' Three model specifications are supported:
#' \describe{
#'   \item{Model "o"}{No structural break - standard cointegration test}
#'   \item{Model "c"}{Break in constant only - shift in intercept at break date}
#'   \item{Model "cs"}{Break in constant and slope - shift in both intercept 
#'     and cointegrating coefficients}
#' }
#'
#' @references
#' Trinh, H. H. (2022). Testing for cointegration with structural changes in 
#' very small sample. THEMA Working Paper n°2022-01, CY Cergy Paris Université.
#' \\url{https://ideas.repec.org/p/ema/worpap/2022-01.html}
#'
#' Gregory, A. W., & Hansen, B. E. (1996). Residual-based tests for 
#' cointegration in models with regime shifts. Journal of Econometrics, 
#' 70(1), 99-126. \doi{10.1016/0304-4076(69)41685-7}
#'
#' Hatemi-J, A. (2008). Tests for cointegration with two unknown regime 
#' shifts with an application to financial market integration. Empirical 
#' Economics, 35(3), 497-505. \doi{10.1007/s00181-007-0175-9}
#'
#' Engle, R. F., & Granger, C. W. J. (1987). Co-integration and error 
#' correction: Representation, estimation, and testing. Econometrica, 
#' 55(2), 251-276. \doi{10.2307/1913236}
#'
#' @docType package
#' @name cointsmall-package
NULL
