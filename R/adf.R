#' Augmented Dickey-Fuller Test for Residuals
#'
#' @description
#' Computes the ADF test statistic for residuals from a cointegrating 
#' regression. Uses BIC to select the optimal lag length.
#'
#' @param x Numeric vector of residuals.
#' @param maxlags Maximum number of lags to consider.
#'
#' @return A list with components:
#'   \item{stat}{ADF t-statistic}
#'   \item{lags}{Selected lag length}
#'
#' @keywords internal
#' @noRd
.adf_test <- function(x, maxlags) {
  n <- length(x)
  
  # Ensure we have enough observations
  maxlags <- min(maxlags, floor((n - 3) / 2))
  if (maxlags < 0) maxlags <- 0
  
  best_bic <- Inf
  best_lag <- 0
  
  # Select optimal lag using BIC
  for (p in 0:maxlags) {
    result <- .adf_regression(x, p)
    if (!is.null(result)) {
      n_eff <- result$n
      k <- p + 1  # rho + p lags
      bic <- n_eff * log(result$ssr / n_eff) + k * log(n_eff)
      if (bic < best_bic) {
        best_bic <- bic
        best_lag <- p
      }
    }
  }
  
  # Run ADF with optimal lag
  result <- .adf_regression(x, best_lag)
  
  list(stat = result$tstat, lags = best_lag)
}

#' ADF Regression (No Constant)
#'
#' @description
#' Runs the ADF regression without constant term (for cointegrating 
#' regression residuals).
#'
#' @param x Numeric vector.
#' @param p Number of lags.
#'
#' @return List with t-statistic, SSR, and effective sample size.
#'
#' @keywords internal
#' @noRd
.adf_regression <- function(x, p) {
  n <- length(x)
  
  # Need at least p + 2 observations
  if (n < p + 3) return(NULL)
  
  # First difference
  dx <- diff(x)
  
  # Lagged level
  x_lag <- x[1:(n - 1)]
  
  # Effective sample size after differencing and lagging
  n_eff <- length(dx) - p
  if (n_eff < 3) return(NULL)
  
  # Build matrices
  # y = d.x_t (from p+1 to n-1 in diff terms)
  y <- dx[(p + 1):length(dx)]
  
  # X = [x_{t-1}, d.x_{t-1}, ..., d.x_{t-p}]
  # x_{t-1} lagged level at time t (in diff index)
  X_level <- x_lag[(p + 1):length(x_lag)]
  
  # Lagged differences
  if (p > 0) {
    X_diff <- matrix(NA, nrow = n_eff, ncol = p)
    for (j in 1:p) {
      X_diff[, j] <- dx[(p + 1 - j):(length(dx) - j)]
    }
    X <- cbind(X_level, X_diff)
  } else {
    X <- matrix(X_level, ncol = 1)
  }
  
  # OLS without constant
  # (X'X)^{-1} X'y
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  
  # Check for singularity
  XtX_inv <- tryCatch(
    solve(XtX),
    error = function(e) NULL
  )
  
  if (is.null(XtX_inv)) return(NULL)
  
  beta <- XtX_inv %*% Xty
  resid <- y - X %*% beta
  ssr <- sum(resid^2)
  
  # Standard error of rho (first coefficient)
  s2 <- ssr / (n_eff - ncol(X))
  se_beta <- sqrt(diag(s2 * XtX_inv))
  
  # t-statistic for rho
  tstat <- beta[1] / se_beta[1]
  
  list(tstat = as.numeric(tstat), ssr = ssr, n = n_eff)
}
