#' Critical Values for Cointegration Tests with Structural Breaks
#'
#' @description
#' Computes the size-corrected 5\% critical values of the ADF* test of no
#' cointegration from the response surfaces estimated by Trinh (2022).
#'
#' @param TT Sample size.
#' @param m Number of independent variables in the cointegrating regression
#'   (1, 2 or 3).
#' @param breaks Number of structural breaks (0, 1, or 2).
#' @param model Model specification ("o", "c", or "cs").
#' @param level Significance level. Only 5 is available, because Trinh (2022)
#'   reports the response surfaces for the 5\% quantile only. If NULL, the
#'   list of all levels is returned, with \code{NA} at 1\% and 10\%.
#'
#' @return If level is specified, returns the critical value. If NULL, returns
#'   a named list with cv01, cv05, and cv10 (cv01 and cv10 are \code{NA}).
#'
#' @details
#' The critical value is
#' \deqn{cv(T) = \psi_\infty + \sum_{k=1}^{K} \psi_k T^{-k},}
#' with the coefficients of Table 13 of Trinh (2022), estimated from 10,000
#' replications for T = 12 to 1,000 with the polynomial order K chosen by
#' AIC (K at most 6). They cover m = 1, 2, 3 regressors, zero breaks
#' (model "o"), and one or two breaks in the constant (model "c") or in the
#' constant and the slopes (model "cs"). The resulting values reproduce
#' Table 1 of the paper (for example -5.40 for m = 1, one break, model "cs"
#' and T = 50).
#'
#' @references
#' Trinh, J. (2022). Testing for cointegration with structural changes in
#' very small sample. THEMA Working Paper 2022-01, CY Cergy Paris
#' Universite.
#' \url{https://ideas.repec.org/p/ema/worpap/2022-01.html}
#'
#' @examples
#' # 5\% critical value for m = 1 regressor, T = 30, no breaks
#' cointsmall_cv(TT = 30, m = 1, breaks = 0, model = "o", level = 5)
#'
#' # One break in constant and slope, m = 2, T = 50
#' cointsmall_cv(TT = 50, m = 2, breaks = 1, model = "cs")
#'
#' @export
cointsmall_cv <- function(TT, m, breaks = 0, model = "o", level = NULL) {
  cv <- .get_critical_values(TT, m, breaks, model)
  if (is.null(level)) {
    return(cv)
  }
  if (!identical(as.numeric(level), 5)) {
    stop("Trinh (2022) reports 5% critical values only; use level = 5.")
  }
  cv$cv05
}

#' Get Critical Values (Internal)
#'
#' Response surface coefficients of Trinh (2022), Table 13 (5% quantile).
#'
#' @keywords internal
#' @noRd
.get_critical_values <- function(TT, m, breaks, model) {
  key <- paste(m, breaks, model)
  coefs <- list(
    "1 0 o"  = c(-3.33, -16.88, 798.01, -30818.40, 460634.58, -2279397.87),
    "2 0 o"  = c(-3.75, -10.25, 80.17, -13337.52, 302551.21, -1848305.35),
    "3 0 o"  = c(-4.10, -12.16, -321.05, 7197.98, -40759.64),
    "1 1 c"  = c(-4.62, -13.05, -1399.49, 76213.27, -1939275.51,
                 23030362.35, -99593635.82),
    "2 1 c"  = c(-4.97, -28.28, 112.01, -3338.30, 48647.86),
    "3 1 c"  = c(-5.30, -40.62, 1759.22, -94294.31, 2287166.84,
                 -25326822.40, 106995759.84),
    "1 2 c"  = c(-5.21, -279.74, 35643.95, -1963265.34, 49133408.78,
                 -564440298.77, 2411884754.22),
    "2 2 c"  = c(-5.53, -287.06, 36360.62, -2001591.63, 50089808.47,
                 -575604522.50, 2460441619.22),
    "3 2 c"  = c(-5.88, -272.48, 34788.63, -1938696.30, 48860977.69,
                 -564638501.27, 2425723155.38),
    "1 1 cs" = c(-4.96, -20.19, -64.18, -1901.05, 45903.20),
    "2 1 cs" = c(-5.55, -29.61, 205.17, -7483.02, 84068.24),
    "3 1 cs" = c(-6.09, -13.81, -2439.38, 125430.43, -2972990.17,
                 32209309.76, -126768011.40),
    "1 2 cs" = c(-5.94, -207.45, 26499.94, -1492615.06, 37796851.87,
                 -437913647.90, 1884665034.50),
    "2 2 cs" = c(-6.90, -57.52, 5352.61, -403380.76, 12234074.42,
                 -164125034.32, 802737634.37),
    "3 2 cs" = c(-7.67, 65.2, -18638.41, 1.263981e6, -4.139881e7,
                 6.391186e8, -3.757791e9)
  )
  if (!key %in% names(coefs)) {
    if (m > 3) {
      stop("Trinh (2022) provides critical values for at most 3 regressors.")
    }
    stop("Invalid combination of 'breaks' and 'model'.")
  }
  psi <- coefs[[key]]
  cv05 <- sum(psi * TT^(-(seq_along(psi) - 1)))
  list(cv01 = NA_real_, cv05 = cv05, cv10 = NA_real_)
}
