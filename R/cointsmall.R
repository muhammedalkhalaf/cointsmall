#' @importFrom stats lm.fit
#'
#' Cointegration Test with Structural Breaks in Small Samples
#'
#' @description
#' Tests for cointegration between a dependent variable and one or more 
#' independent variables, allowing for structural breaks. The method is
#' designed specifically for small sample sizes following Trinh (2022).
#'
#' @param y Numeric vector of the dependent variable (must be I(1)).
#' @param x Numeric vector or matrix of independent variable(s) (must be I(1)).
#' @param breaks Integer specifying the number of structural breaks to test for.
#'   Must be 0, 1, or 2. Default is 1.
#' @param model Character string specifying the model type:
#'   \describe{
#'     \item{"o"}{No structural break (only valid when breaks = 0)}
#'     \item{"c"}{Break in constant only}
#'     \item{"cs"}{Break in constant and slope (default for breaks > 0)}
#'   }
#' @param criterion Character string specifying the criterion for break date 
#'   selection: \code{"adf"} (minimize ADF statistic) or \code{"ssr"} 
#'   (minimize sum of squared residuals). Default is \code{"adf"}.
#' @param trim Numeric value between 0 and 0.5 specifying the trimming 
#'   parameter for break date search. Default is 0.15.
#' @param maxlags Integer specifying the maximum number of lags for the ADF 
#'   test. If -1 (default), automatically determined using the rule 
#'   \code{floor(12*(T/100)^0.25)}.
#' @param level Numeric confidence level for critical values (1, 5, or 10). 
#'   Default is 5.
#'
#' @return An object of class \code{"cointsmall"} containing:
#'   \item{statistic}{The ADF* test statistic}
#'   \item{cv}{Critical value at the specified level}
#'   \item{cv01}{Critical value at 1\% level}
#'   \item{cv05}{Critical value at 5\% level}
#'   \item{cv10}{Critical value at 10\% level}
#'   \item{pvalue}{Approximate p-value}
#'   \item{decision}{Character string with test decision}
#'   \item{reject}{Logical indicating whether to reject null hypothesis}
#'   \item{breaks}{Number of breaks tested}
#'   \item{model}{Model specification used}
#'   \item{criterion}{Selection criterion used}
#'   \item{break_dates}{Vector of estimated break date indices (if breaks > 0)}
#'   \item{lags}{Number of lags used in ADF test}
#'   \item{ssr}{Sum of squared residuals from cointegrating regression}
#'   \item{nobs}{Number of observations}
#'   \item{nvar}{Number of independent variables}
#'   \item{coefficients}{Estimated cointegrating coefficients}
#'   \item{residuals}{Residuals from cointegrating regression}
#'
#' @details
#' The test follows the two-step Engle-Granger procedure with modifications
#' for structural breaks:
#' 
#' \enumerate{
#'   \item Estimate the cointegrating regression (with break dummies if 
#'     applicable)
#'   \item Apply an ADF test to the residuals
#' }
#' 
#' For models with breaks, the break date(s) are determined endogenously
#' by searching over all possible dates within the trimmed sample and 
#' selecting the date that minimizes the ADF statistic or SSR.
#' 
#' Critical values are computed using response surface methodology 
#' following Trinh (2022), which accounts for the small sample bias.
#'
#' @references
#' Trinh, H. H. (2022). Testing for cointegration with structural changes in 
#' very small sample. THEMA Working Paper n°2022-01, CY Cergy Paris Université.
#' \\url{https://ideas.repec.org/p/ema/worpap/2022-01.html}
#'
#' @examples
#' # Generate cointegrated series
#' set.seed(42)
#' n <- 50
#' e <- cumsum(rnorm(n))  # Common stochastic trend
#' y <- 2 + 3 * e + rnorm(n, sd = 0.5)
#' x <- e + rnorm(n, sd = 0.3)
#' 
#' # Test with no break
#' result0 <- cointsmall(y, x, breaks = 0)
#' print(result0)
#' 
#' # Test with one break (break in constant and slope)
#' result1 <- cointsmall(y, x, breaks = 1, model = "cs")
#' print(result1)
#' 
#' # Generate series with structural break
#' y_break <- c(2 + 2 * e[1:25], 5 + 4 * e[26:50]) + rnorm(n, sd = 0.3)
#' result_break <- cointsmall(y_break, x, breaks = 1)
#' print(result_break)
#'
#' @export
cointsmall <- function(y, x, breaks = 1, model = NULL, criterion = "adf",
                       trim = 0.15, maxlags = -1, level = 5) {
  
  # Input validation
  y <- as.numeric(y)
  x <- as.matrix(x)
  
  if (length(y) != nrow(x)) {
    stop("Length of 'y' must equal number of rows in 'x'")
  }
  
  T <- length(y)
  m <- ncol(x)
  
  if (T < 12) {
    stop("Sample size too small. Minimum 12 observations required.")
  }
  
  if (!breaks %in% c(0, 1, 2)) {
    stop("'breaks' must be 0, 1, or 2")
  }
  
  # Set default model
  if (is.null(model)) {
    model <- if (breaks == 0) "o" else "cs"
  }
  
  # Validate model
  if (!model %in% c("o", "c", "cs")) {
    stop("'model' must be 'o', 'c', or 'cs'")
  }
  
  if (model == "o" && breaks > 0) {
    stop("model 'o' is incompatible with breaks > 0")
  }
  
  if (model %in% c("c", "cs") && breaks == 0) {
    stop("models 'c' and 'cs' require breaks > 0")
  }
  
  # Validate criterion
  criterion <- match.arg(criterion, c("adf", "ssr"))
  
  # Validate trim
  if (trim <= 0 || trim >= 0.5) {
    stop("'trim' must be between 0 and 0.5")
  }
  
  # Validate level
  if (!level %in% c(1, 5, 10)) {
    stop("'level' must be 1, 5, or 10")
  }
  
  # Set maxlags
  if (maxlags == -1) {
    maxlags <- floor(12 * (T / 100)^0.25)
  }
  maxlags <- min(maxlags, floor(T / 4))
  
  # Run appropriate test
  if (breaks == 0) {
    result <- .cointsmall_test0(y, x, maxlags, level)
  } else if (breaks == 1) {
    result <- .cointsmall_test1(y, x, model, criterion, trim, maxlags, level)
  } else {
    result <- .cointsmall_test2(y, x, model, criterion, trim, maxlags, level)
  }
  
  # Add common elements
  result$breaks <- breaks
  result$model <- model
  result$criterion <- criterion
  result$nobs <- T
  result$nvar <- m
  result$call <- match.call()
  
  class(result) <- "cointsmall"
  return(result)
}

#' @rdname cointsmall
#' @export
print.cointsmall <- function(x, ...) {
  cat("\nCointegration Test with Structural Breaks in Small Sample\n")
  cat("Based on Trinh (2022)\n")
  cat(strrep("-", 70), "\n")
  
  cat("Sample size:", x$nobs, "\n")
  cat("Number of regressors:", x$nvar, "\n")
  cat("Number of breaks:", x$breaks, "\n")
  cat("Model:", switch(x$model, 
                       "o" = "No break",
                       "c" = "Break in constant",
                       "cs" = "Break in constant and slope"), "\n")
  
  if (x$breaks > 0) {
    cat("Selection criterion:", toupper(x$criterion), "\n")
    cat("Break date(s):", paste(x$break_dates, collapse = ", "), "\n")
  }
  
  cat(strrep("-", 70), "\n")
  cat(sprintf("ADF* statistic: %9.4f\n", x$statistic))
  cat(sprintf("Lags selected:  %9d\n", x$lags))
  cat(sprintf("Critical value (%d%%): %9.4f\n", 
              switch(as.character(which(c(1, 5, 10) == 5)), 
                     "1" = 1, "2" = 5, "3" = 10), 
              x$cv05))
  cat(sprintf("P-value:        %9.4f\n", x$pvalue))
  cat(strrep("-", 70), "\n")
  
  cat("Decision:", x$decision, "\n")
  cat(strrep("-", 70), "\n\n")
  
  invisible(x)
}

#' @rdname cointsmall
#' @export
summary.cointsmall <- function(object, ...) {
  print(object, ...)
  
  cat("\nCointegrating Regression Coefficients:\n")
  print(object$coefficients)
  
  if (object$breaks > 0) {
    cat("\nBreak fraction(s):", 
        paste(round(object$break_dates / object$nobs, 3), collapse = ", "), "\n")
  }
  
  cat("\nCritical Values:\n")
  cat(sprintf("  1%%:  %9.4f\n", object$cv01))
  cat(sprintf("  5%%:  %9.4f\n", object$cv05))
  cat(sprintf("  10%%: %9.4f\n", object$cv10))
  
  invisible(object)
}

# Internal function: Test with no breaks (model o)
.cointsmall_test0 <- function(y, x, maxlags, level) {
  T <- length(y)
  m <- ncol(x)
  
  # Cointegrating regression with constant
  X <- cbind(1, x)
  fit <- stats::lm.fit(X, y)
  resid <- fit$residuals
  coef <- fit$coefficients
  names(coef) <- c("(Intercept)", paste0("x", seq_len(m)))
  
  # ADF test on residuals
  adf_result <- .adf_test(resid, maxlags)
  
  # Get critical values
  cv <- .get_critical_values(T, m, breaks = 0, model = "o")
  
  # Compute p-value
  pval <- .interpolate_pvalue(adf_result$stat, cv)
  
  # Determine rejection
  cv_level <- cv[[paste0("cv", sprintf("%02d", level))]]
  reject <- adf_result$stat < cv_level
  
  decision <- if (reject) {
    "Reject H0 - Evidence of cointegration"
  } else {
    "Do not reject H0 - No evidence of cointegration"
  }
  
  list(
    statistic = adf_result$stat,
    cv = cv_level,
    cv01 = cv$cv01,
    cv05 = cv$cv05,
    cv10 = cv$cv10,
    pvalue = pval,
    decision = decision,
    reject = reject,
    lags = adf_result$lags,
    ssr = sum(resid^2),
    break_dates = NULL,
    coefficients = coef,
    residuals = resid
  )
}

# Internal function: Test with one break
.cointsmall_test1 <- function(y, x, model, criterion, trim, maxlags, level) {
  T <- length(y)
  m <- ncol(x)
  
  # Determine search range
  t1_min <- ceiling(T * trim)
  t1_max <- floor(T * (1 - trim))
  
  if (t1_min < 2) t1_min <- 2
  if (t1_max > T - 1) t1_max <- T - 1
  
  # Initialize search
  best_adf <- Inf
  best_t1 <- NA
  best_lags <- NA
  best_ssr <- Inf
  best_resid <- NULL
  best_coef <- NULL
  
  # Search over break dates
  for (t1 in t1_min:t1_max) {
    # Create break dummy
    D1 <- as.numeric(seq_len(T) >= t1)
    
    # Build design matrix
    if (model == "c") {
      X <- cbind(1, x, D1)
    } else {  # model == "cs"
      interactions <- x * D1
      X <- cbind(1, x, D1, interactions)
    }
    
    # Fit regression
    fit <- stats::lm.fit(X, y)
    resid <- fit$residuals
    ssr <- sum(resid^2)
    
    # ADF test
    adf_result <- .adf_test(resid, maxlags)
    
    # Update best based on criterion
    update <- if (criterion == "adf") {
      adf_result$stat < best_adf
    } else {
      ssr < best_ssr
    }
    
    if (update) {
      best_adf <- adf_result$stat
      best_t1 <- t1
      best_lags <- adf_result$lags
      best_ssr <- ssr
      best_resid <- resid
      best_coef <- fit$coefficients
    }
  }
  
  # Name coefficients
  if (model == "c") {
    names(best_coef) <- c("(Intercept)", paste0("x", seq_len(m)), "D1")
  } else {
    names(best_coef) <- c("(Intercept)", paste0("x", seq_len(m)), "D1",
                          paste0("x", seq_len(m), ":D1"))
  }
  
  # Get critical values
  cv <- .get_critical_values(T, m, breaks = 1, model = model)
  
  # Compute p-value
  pval <- .interpolate_pvalue(best_adf, cv)
  
  # Determine rejection
  cv_level <- cv[[paste0("cv", sprintf("%02d", level))]]
  reject <- best_adf < cv_level
  
  decision <- if (reject) {
    "Reject H0 - Evidence of cointegration with structural break"
  } else {
    "Do not reject H0 - No evidence of cointegration"
  }
  
  list(
    statistic = best_adf,
    cv = cv_level,
    cv01 = cv$cv01,
    cv05 = cv$cv05,
    cv10 = cv$cv10,
    pvalue = pval,
    decision = decision,
    reject = reject,
    lags = best_lags,
    ssr = best_ssr,
    break_dates = best_t1,
    coefficients = best_coef,
    residuals = best_resid
  )
}

# Internal function: Test with two breaks
.cointsmall_test2 <- function(y, x, model, criterion, trim, maxlags, level) {
  T <- length(y)
  m <- ncol(x)
  
  # Determine search range
  t_min <- ceiling(T * trim)
  t_max <- floor(T * (1 - trim))
  
  if (t_min < 2) t_min <- 2
  if (t_max > T - 1) t_max <- T - 1
  
  # Initialize search
  best_adf <- Inf
  best_t1 <- NA
  best_t2 <- NA
  best_lags <- NA
  best_ssr <- Inf
  best_resid <- NULL
  best_coef <- NULL
  
  # Search over break dates
  for (t1 in t_min:(t_max - 1)) {
    for (t2 in (t1 + 1):t_max) {
      # Create break dummies
      D1 <- as.numeric(seq_len(T) >= t1)
      D2 <- as.numeric(seq_len(T) >= t2)
      
      # Build design matrix
      if (model == "c") {
        X <- cbind(1, x, D1, D2)
      } else {  # model == "cs"
        int1 <- x * D1
        int2 <- x * D2
        X <- cbind(1, x, D1, D2, int1, int2)
      }
      
      # Fit regression
      fit <- stats::lm.fit(X, y)
      resid <- fit$residuals
      ssr <- sum(resid^2)
      
      # ADF test
      adf_result <- .adf_test(resid, maxlags)
      
      # Update best based on criterion
      update <- if (criterion == "adf") {
        adf_result$stat < best_adf
      } else {
        ssr < best_ssr
      }
      
      if (update) {
        best_adf <- adf_result$stat
        best_t1 <- t1
        best_t2 <- t2
        best_lags <- adf_result$lags
        best_ssr <- ssr
        best_resid <- resid
        best_coef <- fit$coefficients
      }
    }
  }
  
  # Name coefficients
  if (model == "c") {
    names(best_coef) <- c("(Intercept)", paste0("x", seq_len(m)), "D1", "D2")
  } else {
    names(best_coef) <- c("(Intercept)", paste0("x", seq_len(m)), "D1", "D2",
                          paste0("x", seq_len(m), ":D1"),
                          paste0("x", seq_len(m), ":D2"))
  }
  
  # Get critical values
  cv <- .get_critical_values(T, m, breaks = 2, model = model)
  
  # Compute p-value
  pval <- .interpolate_pvalue(best_adf, cv)
  
  # Determine rejection
  cv_level <- cv[[paste0("cv", sprintf("%02d", level))]]
  reject <- best_adf < cv_level
  
  decision <- if (reject) {
    "Reject H0 - Evidence of cointegration with structural breaks"
  } else {
    "Do not reject H0 - No evidence of cointegration"
  }
  
  list(
    statistic = best_adf,
    cv = cv_level,
    cv01 = cv$cv01,
    cv05 = cv$cv05,
    cv10 = cv$cv10,
    pvalue = pval,
    decision = decision,
    reject = reject,
    lags = best_lags,
    ssr = best_ssr,
    break_dates = c(best_t1, best_t2),
    coefficients = best_coef,
    residuals = best_resid
  )
}
