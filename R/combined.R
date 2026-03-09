#' Combined Cointegration Testing Procedure
#'
#' @description
#' Performs cointegration tests under all model specifications (no break,
#' break in constant, break in constant and slope) and provides model
#' selection guidance.
#'
#' @param y Numeric vector of the dependent variable (must be I(1)).
#' @param x Numeric vector or matrix of independent variable(s) (must be I(1)).
#' @param breaks Integer specifying the number of structural breaks to test for.
#'   Must be 1 or 2. Default is 1.
#' @param trim Numeric value between 0 and 0.5 specifying the trimming 
#'   parameter for break date search. Default is 0.15.
#' @param maxlags Integer specifying the maximum number of lags for the ADF 
#'   test. If -1 (default), automatically determined.
#' @param level Numeric confidence level for critical values (1, 5, or 10). 
#'   Default is 5.
#'
#' @return An object of class \code{"cointsmall_combined"} containing:
#'   \item{results}{List of cointsmall objects for each model}
#'   \item{summary}{Data frame summarizing test statistics and decisions}
#'   \item{selected_model}{Character string indicating the selected model}
#'   \item{nobs}{Number of observations}
#'   \item{nvar}{Number of independent variables}
#'   \item{breaks}{Number of breaks tested}
#'
#' @details
#' The combined procedure tests three model specifications:
#' \enumerate{
#'   \item Model "o": No structural break
#'   \item Model "c": Break in constant only
#'   \item Model "cs": Break in constant and slope
#' }
#'
#' Model selection follows these rules:
#' \itemize{
#'   \item If no model rejects H0: No evidence of cointegration
#'   \item If exactly one model rejects H0: Select that model
#'   \item If multiple models reject H0: Select the most general model that
#'     rejects H0 (cs > c > o)
#' }
#'
#' @references
#' Trinh, H. H. (2022). Testing for cointegration with structural changes in 
#' very small sample. THEMA Working Paper n°2022-01, CY Cergy Paris Université.
#' \\url{https://ideas.repec.org/p/ema/worpap/2022-01.html}
#'
#' @examples
#' # Generate cointegrated series with break
#' set.seed(123)
#' n <- 50
#' e <- cumsum(rnorm(n))
#' x <- e + rnorm(n, sd = 0.3)
#' y <- c(2 + 2 * e[1:25], 5 + 4 * e[26:50]) + rnorm(n, sd = 0.3)
#' 
#' # Combined test
#' result <- cointsmall_combined(y, x, breaks = 1)
#' print(result)
#'
#' @export
cointsmall_combined <- function(y, x, breaks = 1, trim = 0.15, 
                                 maxlags = -1, level = 5) {
  
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
  
  if (!breaks %in% c(1, 2)) {
    stop("'breaks' must be 1 or 2 for combined testing")
  }
  
  # Test all models
  result_o <- cointsmall(y, x, breaks = 0, model = "o", 
                          trim = trim, maxlags = maxlags, level = level)
  result_c <- cointsmall(y, x, breaks = breaks, model = "c", criterion = "adf",
                          trim = trim, maxlags = maxlags, level = level)
  result_cs <- cointsmall(y, x, breaks = breaks, model = "cs", criterion = "adf",
                           trim = trim, maxlags = maxlags, level = level)
  
  # Create summary
  summary_df <- data.frame(
    Model = c("o (no break)", "c (break in constant)", "cs (break in constant + slope)"),
    ADF_stat = c(result_o$statistic, result_c$statistic, result_cs$statistic),
    CV = c(result_o$cv, result_c$cv, result_cs$cv),
    Reject_H0 = c(result_o$reject, result_c$reject, result_cs$reject),
    stringsAsFactors = FALSE
  )
  
  # Model selection
  n_reject <- sum(c(result_o$reject, result_c$reject, result_cs$reject))
  
  if (n_reject == 0) {
    selected <- "none"
    conclusion <- "No evidence of cointegration under any model specification"
  } else if (n_reject == 1) {
    if (result_o$reject) {
      selected <- "o"
      conclusion <- "Evidence of cointegration without structural break"
    } else if (result_c$reject) {
      selected <- "c"
      conclusion <- "Evidence of cointegration with break in constant"
    } else {
      selected <- "cs"
      conclusion <- "Evidence of cointegration with break in constant and slope"
    }
  } else {
    # Multiple rejections - select most general model that rejects
    if (result_cs$reject) {
      selected <- "cs"
      conclusion <- "Evidence of cointegration with regime change (break in constant and slope)"
    } else if (result_c$reject) {
      selected <- "c"
      conclusion <- "Evidence of cointegration with level shift (break in constant)"
    } else {
      selected <- "o"
      conclusion <- "Evidence of cointegration without structural break"
    }
  }
  
  result <- list(
    results = list(o = result_o, c = result_c, cs = result_cs),
    summary = summary_df,
    selected_model = selected,
    conclusion = conclusion,
    nobs = T,
    nvar = m,
    breaks = breaks,
    level = level,
    call = match.call()
  )
  
  class(result) <- "cointsmall_combined"
  return(result)
}

#' @rdname cointsmall_combined
#' @export
print.cointsmall_combined <- function(x, ...) {
  cat("\nCombined Cointegration Testing Procedure\n")
  cat("Based on Trinh (2022)\n")
  cat(strrep("=", 70), "\n")
  
  cat("Sample size:", x$nobs, "\n")
  cat("Number of regressors:", x$nvar, "\n")
  cat("Number of breaks tested:", x$breaks, "\n")
  cat("Significance level:", x$level, "%\n")
  cat(strrep("-", 70), "\n\n")
  
  cat("Test Results:\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("%-35s %12s %12s %10s\n", "Model", "ADF* stat", "Crit. val", "Reject?"))
  cat(strrep("-", 70), "\n")
  
  for (i in 1:nrow(x$summary)) {
    cat(sprintf("%-35s %12.4f %12.4f %10s\n",
                x$summary$Model[i],
                x$summary$ADF_stat[i],
                x$summary$CV[i],
                ifelse(x$summary$Reject_H0[i], "Yes", "No")))
  }
  
  cat(strrep("-", 70), "\n\n")
  
  # Break dates for models with breaks
  if (x$breaks > 0) {
    cat("Estimated Break Date(s):\n")
    cat("  Model c: ", paste(x$results$c$break_dates, collapse = ", "), "\n")
    cat("  Model cs:", paste(x$results$cs$break_dates, collapse = ", "), "\n")
    cat("\n")
  }
  
  cat("Model Selection:\n")
  cat(strrep("-", 70), "\n")
  cat("Selected model:", x$selected_model, "\n")
  cat("Conclusion:", x$conclusion, "\n")
  cat(strrep("=", 70), "\n\n")
  
  invisible(x)
}

#' @rdname cointsmall_combined
#' @export
summary.cointsmall_combined <- function(object, ...) {
  print(object, ...)
  
  cat("\nDetailed Results by Model:\n")
  cat(strrep("=", 70), "\n\n")
  
  for (model_name in c("o", "c", "cs")) {
    cat("Model", model_name, ":\n")
    cat(strrep("-", 40), "\n")
    res <- object$results[[model_name]]
    cat(sprintf("  ADF statistic: %9.4f\n", res$statistic))
    cat(sprintf("  Critical values: 1%%: %.4f, 5%%: %.4f, 10%%: %.4f\n",
                res$cv01, res$cv05, res$cv10))
    cat(sprintf("  P-value: %.4f\n", res$pvalue))
    cat(sprintf("  Lags: %d\n", res$lags))
    if (!is.null(res$break_dates)) {
      cat(sprintf("  Break date(s): %s\n", paste(res$break_dates, collapse = ", ")))
    }
    cat("\n")
  }
  
  invisible(object)
}
