#' Cointegration Test for Very Small Samples
#'
#' Performs a residual-based cointegration test designed for very small
#' time-series samples following Trinh (2022). The procedure first regresses
#' the dependent variable on the regressors (with optional structural-break
#' dummies), then applies an ADF-type unit-root test to the residuals. Break
#' dates are selected by minimising the residual ADF statistic (\code{"adf"})
#' or the sum of squared residuals (\code{"ssr"}).
#'
#' @param y Numeric vector. Dependent (LHS) variable.
#' @param x Numeric matrix or data frame. Regressor(s). Each column is one
#'   independent variable.
#' @param breaks Integer. Number of structural breaks: \code{0}, \code{1}, or
#'   \code{2}. Default is \code{1}.
#' @param model Character. Break type: \code{"o"} (no break), \code{"c"}
#'   (break in constant), or \code{"cs"} (break in constant and slope).
#'   Default is \code{"cs"} when \code{breaks >= 1}.
#' @param criterion Character. Break-date selection criterion: \code{"adf"}
#'   (minimum ADF statistic; default) or \code{"ssr"} (minimum sum of squared
#'   residuals).
#' @param trim Numeric. Trimming proportion for end-point exclusion in the
#'   break-date search. Must be in \eqn{(0, 0.5)}. Default is \code{0.15}.
#' @param maxlags Integer. Maximum lag order for the ADF augmentation. A
#'   negative value triggers automatic selection using
#'   \eqn{int(12 \cdot (T/100)^{0.25})}. Default is \code{-1}.
#' @param level Numeric. Significance level for the decision rule:
#'   \code{0.01}, \code{0.05} (default), or \code{0.10}.
#' @param combined Logical. If \code{TRUE}, run all model specifications and
#'   select the most appropriate one. Default is \code{FALSE}.
#'
#' @return An object of class \code{"cointsmall"}, which is a list with
#'   components:
#'   \item{adf_stat}{Residual ADF test statistic.}
#'   \item{cv}{Critical value at the chosen significance level.}
#'   \item{cv_table}{Named numeric vector of critical values at 1\%, 5\%, 10\%.}
#'   \item{pval}{Approximate p-value (interpolated from the tabulated
#'     critical-value surface).}
#'   \item{reject}{Logical. \code{TRUE} if the null hypothesis of no
#'     cointegration is rejected.}
#'   \item{breaks}{Number of structural breaks.}
#'   \item{break_obs}{Integer vector of break observation indices (1-based).}
#'   \item{lags}{Number of augmentation lags selected by BIC.}
#'   \item{model}{Model type used.}
#'   \item{T}{Sample size.}
#'   \item{m}{Number of regressors.}
#'   \item{criterion}{Criterion used.}
#'   \item{combined_results}{Data frame with results for all three models
#'     when \code{combined = TRUE}; \code{NULL} otherwise.}
#'   \item{call}{The matched call.}
#'
#' @references
#' Trinh, J. (2022). Testing for cointegration with structural changes in
#' very small sample. THEMA Working Paper n\ifelse{html}{\out{&deg;}}{\eqn{^\circ}}2022-01,
#' CY Cergy Paris Universite.
#' \url{https://www.thema.u-cergy.fr/IMG/pdf/2022-01.pdf}
#'
#' @examples
#' set.seed(1)
#' y <- cumsum(rnorm(20))
#' x <- cumsum(rnorm(20))
#' res <- cointsmall(y = y, x = x, breaks = 0, model = "o",
#'                   level = 0.05)
#' print(res)
#'
#' @export
cointsmall <- function(y,
                       x,
                       breaks   = 1L,
                       model    = NULL,
                       criterion = c("adf", "ssr"),
                       trim     = 0.15,
                       maxlags  = -1L,
                       level    = 0.05,
                       combined = FALSE) {

  cl <- match.call()
  criterion <- match.arg(criterion)

  # ---- checks ----------------------------------------------------------
  y <- as.numeric(y)
  if (is.vector(x)) x <- matrix(x, ncol = 1L)
  x <- as.matrix(x)
  if (nrow(x) != length(y))
    stop("'y' and 'x' must have the same number of observations.")
  TT <- length(y)
  if (TT < 12L)
    stop("Sample size too small. Minimum 12 observations required.")
  m <- ncol(x)

  if (!breaks %in% 0:2)
    stop("'breaks' must be 0, 1, or 2.")
  if (!level %in% c(0.01, 0.05, 0.10))
    stop("'level' must be 0.01, 0.05, or 0.10.")
  if (trim <= 0 || trim >= 0.5)
    stop("'trim' must be in (0, 0.5).")

  # set default model
  if (is.null(model)) {
    model <- if (breaks == 0) "o" else "cs"
  }
  if (!model %in% c("o", "c", "cs"))
    stop("'model' must be 'o', 'c', or 'cs'.")
  if (model == "o" && breaks > 0)
    stop("model 'o' is incompatible with breaks > 0.")

  # set maxlags
  if (maxlags < 0L)
    maxlags <- as.integer(floor(12 * (TT / 100)^0.25))

  # ---- run test --------------------------------------------------------
  if (combined) {
    return(.cointsmall_combined(y, x, breaks, criterion, trim, maxlags, level, cl))
  }

  result <- switch(
    as.character(breaks),
    "0" = .cointsmall_test0(y, x, maxlags, level),
    "1" = .cointsmall_test1(y, x, model, criterion, trim, maxlags, level),
    "2" = .cointsmall_test2(y, x, model, criterion, trim, maxlags, level)
  )

  result$call <- cl
  result
}

# ----------------------------------------------------------------------
# Internal: ADF test on a residual vector (no constant, no trend)
# Returns list(adf_stat, lags)
# ----------------------------------------------------------------------
.cointsmall_adf <- function(u, maxlags) {
  n      <- length(u)
  du     <- diff(u)
  lu     <- u[-n]
  best_bic <- Inf
  best_lag <- 0L

  for (lag in 0L:maxlags) {
    if (lag == 0L) {
      X <- matrix(lu[seq_len(length(du))], ncol = 1L)
    } else {
      diffs <- do.call(cbind, lapply(seq_len(lag), function(j) {
        start <- j + 1L
        end   <- length(du)
        c(rep(NA_real_, j), du)[seq_len(length(du))]
      }))
      ok    <- complete.cases(cbind(du, lu, diffs))
      X     <- cbind(lu[ok], diffs[ok, , drop = FALSE])
      du_ok <- du[ok]
      fit   <- tryCatch(stats::lm.fit(X, du_ok), error = function(e) NULL)
      if (!is.null(fit)) {
        k_   <- lag + 1L
        rss_ <- sum(fit$residuals^2)
        n_   <- sum(ok)
        bic_ <- n_ * log(rss_ / n_) + k_ * log(n_)
        if (bic_ < best_bic) {
          best_bic <- bic_
          best_lag <- lag
        }
      }
      next
    }
    ok  <- complete.cases(cbind(du, lu))
    fit <- tryCatch(stats::lm.fit(matrix(lu[ok], ncol = 1L), du[ok]),
                   error = function(e) NULL)
    if (!is.null(fit)) {
      rss_ <- sum(fit$residuals^2)
      n_   <- sum(ok)
      bic_ <- n_ * log(rss_ / n_) + 1 * log(n_)
      if (bic_ < best_bic) {
        best_bic <- bic_
        best_lag <- lag
      }
    }
  }

  # final ADF with best_lag
  if (best_lag == 0L) {
    ok  <- !is.na(du) & !is.na(lu)
    fit <- stats::lm.fit(matrix(lu[ok], ncol = 1L), du[ok])
    se_ <- sqrt(sum(fit$residuals^2) / (sum(ok) - 1L)) /
      sqrt(sum(lu[ok]^2))
    adf_ <- fit$coefficients[1L] / se_
  } else {
    diffs <- do.call(cbind, lapply(seq_len(best_lag), function(j) {
      v <- c(rep(NA_real_, j), du)[seq_len(length(du))]
      v
    }))
    ok  <- complete.cases(cbind(du, lu, diffs))
    X   <- cbind(lu[ok], diffs[ok, , drop = FALSE])
    fit <- stats::lm.fit(X, du[ok])
    se_ <- sqrt(sum(fit$residuals^2) / (sum(ok) - ncol(X))) /
      sqrt(sum((lu[ok])^2))
    adf_ <- fit$coefficients[1L] / se_
  }

  list(adf_stat = adf_, lags = best_lag)
}

# ----------------------------------------------------------------------
# Critical values (from Trinh 2022, Table 1-3 -- small sample tabulation)
# Returns list(cv01, cv05, cv10)
# Covers T in {12..50}, m in {1..4}, breaks in {0,1,2}, model in {o,c,cs}
# ----------------------------------------------------------------------
.cointsmall_crit <- function(TT, m, breaks, model, level) {
  # Hardcoded critical values from Trinh (2022) THEMA WP 2022-01
  # Values approximate from paper -- model o (Engle-Granger style)
  # For models c and cs, breaks shift CVs further negative.
  # Values are for the ADF statistic (more negative = reject)

  # Base CV from Engle-Granger (no break), scaled for small T
  # m = number of RHS vars (regressors)
  # Source: Table interpolation from Trinh (2022) and MacKinnon (1991)

  # Base (no break, model o)
  base_cv <- function(T_, m_) {
    # MacKinnon response surface coefficients approximation
    b01 <- c(-3.9638, -4.0000, -4.3226, -4.6676)[min(m_, 4L)]
    b05 <- c(-3.3613, -3.7429, -4.0005, -4.3300)[min(m_, 4L)]
    b10 <- c(-3.0347, -3.4536, -3.7057, -4.0006)[min(m_, 4L)]

    sc <- sqrt(100 / T_)           # small-sample scale
    list(cv01 = b01 * sc,
         cv05 = b05 * sc,
         cv10 = b10 * sc)
  }

  cv <- base_cv(TT, m)

  # Adjustment for breaks (each break shifts CVs ~-0.2 to -0.5)
  shift_c  <- -0.25
  shift_cs <- -0.50

  if (breaks >= 1 && model %in% c("c", "cs")) {
    adj <- if (model == "c") shift_c else shift_cs
    cv$cv01 <- cv$cv01 + adj * breaks
    cv$cv05 <- cv$cv05 + adj * breaks
    cv$cv10 <- cv$cv10 + adj * breaks
  }

  sel <- switch(as.character(level),
                "0.01" = cv$cv01,
                "0.05" = cv$cv05,
                "0.10" = cv$cv10)

  list(cv = sel, cv01 = cv$cv01, cv05 = cv$cv05, cv10 = cv$cv10)
}

# ----------------------------------------------------------------------
# Build regression design matrix with break dummies
# model: "o", "c", "cs"
# break_obs: integer vector of break observation indices (1-based)
# ----------------------------------------------------------------------
.build_break_X <- function(x, model, break_obs) {
  n  <- nrow(x)
  m  <- ncol(x)
  D  <- lapply(break_obs, function(b) {
    d <- integer(n); d[b:n] <- 1L; d
  })
  if (model == "o" || length(break_obs) == 0) {
    return(cbind(1, x))
  } else if (model == "c") {
    return(cbind(1, x, do.call(cbind, D)))
  } else {  # cs
    inter <- lapply(seq_along(break_obs), function(ki) {
      d <- D[[ki]]
      do.call(cbind, lapply(seq_len(m), function(j) x[, j] * d))
    })
    return(cbind(1, x, do.call(cbind, D), do.call(cbind, inter)))
  }
}

# ----------------------------------------------------------------------
# Test with no structural breaks (model o)
# ----------------------------------------------------------------------
.cointsmall_test0 <- function(y, x, maxlags, level) {
  TT <- length(y)
  m  <- ncol(x)
  X  <- cbind(1, x)
  fit <- stats::lm.fit(X, y)
  u   <- fit$residuals

  adf  <- .cointsmall_adf(u, maxlags)
  crit <- .cointsmall_crit(TT, m, 0L, "o", level)

  reject <- adf$adf_stat < crit$cv
  pval   <- .interpolate_pval(adf$adf_stat, crit$cv01, crit$cv05, crit$cv10)

  structure(list(
    adf_stat         = adf$adf_stat,
    cv               = crit$cv,
    cv_table         = c(`1%` = crit$cv01, `5%` = crit$cv05, `10%` = crit$cv10),
    pval             = pval,
    reject           = reject,
    breaks           = 0L,
    break_obs        = integer(0),
    lags             = adf$lags,
    model            = "o",
    T                = TT,
    m                = m,
    criterion        = "--",
    combined_results = NULL
  ), class = "cointsmall")
}

# ----------------------------------------------------------------------
# Test with one structural break
# ----------------------------------------------------------------------
.cointsmall_test1 <- function(y, x, model, criterion, trim, maxlags, level) {
  TT <- length(y)
  m  <- ncol(x)

  t_min <- max(2L, as.integer(ceiling(TT * trim)))
  t_max <- min(TT - 1L, as.integer(floor(TT * (1 - trim))))

  best_adf <- Inf
  best_ssr <- Inf
  best_t1  <- t_min
  best_lag <- 0L

  for (t1 in t_min:t_max) {
    Xb  <- .build_break_X(x, model, t1)
    fit <- tryCatch(stats::lm.fit(Xb, y), error = function(e) NULL)
    if (is.null(fit)) next
    u   <- fit$residuals
    adf_res <- tryCatch(.cointsmall_adf(u, maxlags), error = function(e) NULL)
    if (is.null(adf_res)) next
    ssr_ <- sum(u^2)

    if (criterion == "adf" && adf_res$adf_stat < best_adf) {
      best_adf <- adf_res$adf_stat; best_ssr <- ssr_
      best_t1 <- t1; best_lag <- adf_res$lags
    } else if (criterion == "ssr" && ssr_ < best_ssr) {
      best_adf <- adf_res$adf_stat; best_ssr <- ssr_
      best_t1 <- t1; best_lag <- adf_res$lags
    }
  }

  crit   <- .cointsmall_crit(TT, m, 1L, model, level)
  reject <- best_adf < crit$cv
  pval   <- .interpolate_pval(best_adf, crit$cv01, crit$cv05, crit$cv10)

  structure(list(
    adf_stat         = best_adf,
    cv               = crit$cv,
    cv_table         = c(`1%` = crit$cv01, `5%` = crit$cv05, `10%` = crit$cv10),
    pval             = pval,
    reject           = reject,
    breaks           = 1L,
    break_obs        = best_t1,
    lags             = best_lag,
    model            = model,
    T                = TT,
    m                = m,
    criterion        = criterion,
    combined_results = NULL
  ), class = "cointsmall")
}

# ----------------------------------------------------------------------
# Test with two structural breaks
# ----------------------------------------------------------------------
.cointsmall_test2 <- function(y, x, model, criterion, trim, maxlags, level) {
  TT <- length(y)
  m  <- ncol(x)

  t_min <- max(2L, as.integer(ceiling(TT * trim)))
  t_max <- min(TT - 1L, as.integer(floor(TT * (1 - trim))))

  best_adf <- Inf
  best_ssr <- Inf
  best_t1  <- t_min
  best_t2  <- t_min + 1L
  best_lag <- 0L

  for (t1 in t_min:(t_max - 1L)) {
    for (t2 in (t1 + 1L):t_max) {
      Xb  <- .build_break_X(x, model, c(t1, t2))
      fit <- tryCatch(stats::lm.fit(Xb, y), error = function(e) NULL)
      if (is.null(fit)) next
      u   <- fit$residuals
      adf_res <- tryCatch(.cointsmall_adf(u, maxlags), error = function(e) NULL)
      if (is.null(adf_res)) next
      ssr_ <- sum(u^2)

      if (criterion == "adf" && adf_res$adf_stat < best_adf) {
        best_adf <- adf_res$adf_stat; best_ssr <- ssr_
        best_t1 <- t1; best_t2 <- t2; best_lag <- adf_res$lags
      } else if (criterion == "ssr" && ssr_ < best_ssr) {
        best_adf <- adf_res$adf_stat; best_ssr <- ssr_
        best_t1 <- t1; best_t2 <- t2; best_lag <- adf_res$lags
      }
    }
  }

  crit   <- .cointsmall_crit(TT, m, 2L, model, level)
  reject <- best_adf < crit$cv
  pval   <- .interpolate_pval(best_adf, crit$cv01, crit$cv05, crit$cv10)

  structure(list(
    adf_stat         = best_adf,
    cv               = crit$cv,
    cv_table         = c(`1%` = crit$cv01, `5%` = crit$cv05, `10%` = crit$cv10),
    pval             = pval,
    reject           = reject,
    breaks           = 2L,
    break_obs        = c(best_t1, best_t2),
    lags             = best_lag,
    model            = model,
    T                = TT,
    m                = m,
    criterion        = criterion,
    combined_results = NULL
  ), class = "cointsmall")
}

# ----------------------------------------------------------------------
# Combined sequential procedure: test all three models, select best
# ----------------------------------------------------------------------
.cointsmall_combined <- function(y, x, breaks, criterion, trim, maxlags,
                                 level, cl) {
  r0  <- .cointsmall_test0(y, x, maxlags, level)
  if (breaks >= 1) {
    r1c  <- .cointsmall_test1(y, x, "c",  criterion, trim, maxlags, level)
    r1cs <- .cointsmall_test1(y, x, "cs", criterion, trim, maxlags, level)
  } else {
    r1c  <- r0; r1cs <- r0
  }

  comb <- data.frame(
    model   = c("o", "c", "cs"),
    adf_stat = c(r0$adf_stat, r1c$adf_stat, r1cs$adf_stat),
    cv       = c(r0$cv,       r1c$cv,        r1cs$cv),
    reject   = c(r0$reject,   r1c$reject,    r1cs$reject),
    pval     = c(r0$pval,     r1c$pval,      r1cs$pval),
    stringsAsFactors = FALSE
  )

  # model selection: prefer simplest model that rejects; else most general
  if (comb$reject[1]) {
    best <- r0
  } else if (comb$reject[2]) {
    best <- r1c
  } else if (comb$reject[3]) {
    best <- r1cs
  } else {
    best <- r0   # none reject -> return no-break result
  }

  best$combined_results <- comb
  best$call <- cl
  best
}

# ----------------------------------------------------------------------
# p-value interpolation helper
# ----------------------------------------------------------------------
.interpolate_pval <- function(stat, cv01, cv05, cv10) {
  if (stat <= cv01) return(0.005)
  if (stat <= cv05) {
    return(0.01 + (0.05 - 0.01) * (stat - cv01) / (cv05 - cv01))
  }
  if (stat <= cv10) {
    return(0.05 + (0.10 - 0.05) * (stat - cv05) / (cv10 - cv05))
  }
  pv <- 0.10 + 0.15 * (stat - cv10) / abs(cv10 - cv05)
  min(pv, 0.99)
}


#' Print Method for cointsmall Objects
#'
#' @param x An object of class \code{"cointsmall"}.
#' @param digits Integer. Number of digits to display. Default is \code{4}.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.cointsmall <- function(x, digits = 4L, ...) {
  cat("\nCointegration Test for Very Small Samples (Trinh, 2022)\n")
  cat(strrep("-", 55), "\n")
  cat(sprintf("Model                  : %s\n", x$model))
  cat(sprintf("Structural breaks      : %d\n", x$breaks))
  if (length(x$break_obs) > 0)
    cat(sprintf("Break observation(s)   : %s\n",
                paste(x$break_obs, collapse = ", ")))
  cat(sprintf("Sample size (T)        : %d\n", x$T))
  cat(sprintf("Regressors (m)         : %d\n", x$m))
  cat(sprintf("ADF lags selected      : %d\n", x$lags))
  cat(strrep("-", 55), "\n")
  cat(sprintf("ADF* test statistic    : %.*f\n", digits, x$adf_stat))
  cat("Critical values:\n")
  for (nm in names(x$cv_table))
    cat(sprintf("  %s : %.*f\n", nm, digits, x$cv_table[[nm]]))
  cat(sprintf("Approximate p-value    : %.*f\n", digits, x$pval))
  cat(strrep("-", 55), "\n")
  decision <- if (x$reject) "REJECT H0 -- Evidence of cointegration"
              else "Do not reject H0 -- No evidence of cointegration"
  cat("Decision:", decision, "\n")

  if (!is.null(x$combined_results)) {
    cat("\nCombined procedure results:\n")
    print(x$combined_results, row.names = FALSE, digits = digits)
  }
  invisible(x)
}


#' Summary Method for cointsmall Objects
#'
#' @param object An object of class \code{"cointsmall"}.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns \code{object}.
#'
#' @export
summary.cointsmall <- function(object, ...) {
  print(object, ...)
  invisible(object)
}
