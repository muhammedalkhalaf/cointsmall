\name{cointsmall}
\alias{cointsmall}
\title{Cointegration Test for Very Small Samples}
\usage{
cointsmall(
  y,
  x,
  breaks = 1L,
  model = NULL,
  criterion = c("adf", "ssr"),
  trim = 0.15,
  maxlags = -1L,
  level = 0.05,
  combined = FALSE
)
}
\arguments{
\item{y}{Numeric vector. Dependent (LHS) variable.}

\item{x}{Numeric matrix or data frame. Regressor(s). Each column is one
independent variable.}

\item{breaks}{Integer. Number of structural breaks: \code{0}, \code{1}, or
\code{2}. Default is \code{1}.}

\item{model}{Character. Break type: \code{"o"} (no break), \code{"c"}
(break in constant), or \code{"cs"} (break in constant and slope).
Default is \code{"cs"} when \code{breaks >= 1}.}

\item{criterion}{Character. Break-date selection criterion: \code{"adf"}
(minimum ADF statistic; default) or \code{"ssr"} (minimum sum of squared
residuals).}

\item{trim}{Numeric. Trimming proportion for end-point exclusion in the
break-date search. Must be in \eqn{(0, 0.5)}. Default is \code{0.15}.}

\item{maxlags}{Integer. Maximum lag order for the ADF augmentation. A
negative value triggers automatic selection using
\eqn{int(12 \cdot (T/100)^{0.25})}. Default is \code{-1}.}

\item{level}{Numeric. Significance level for the decision rule:
\code{0.01}, \code{0.05} (default), or \code{0.10}.}

\item{combined}{Logical. If \code{TRUE}, run all model specifications and
select the most appropriate one. Default is \code{FALSE}.}
}
\value{
An object of class \code{"cointsmall"}, which is a list with
  components:
  \item{adf_stat}{Residual ADF test statistic.}
  \item{cv}{Critical value at the chosen significance level.}
  \item{cv_table}{Named numeric vector of critical values at 1\%, 5\%, 10\%.}
  \item{pval}{Approximate p-value (interpolated from the tabulated
    critical-value surface).}
  \item{reject}{Logical. \code{TRUE} if the null hypothesis of no
    cointegration is rejected.}
  \item{breaks}{Number of structural breaks.}
  \item{break_obs}{Integer vector of break observation indices (1-based).}
  \item{lags}{Number of augmentation lags selected by BIC.}
  \item{model}{Model type used.}
  \item{T}{Sample size.}
  \item{m}{Number of regressors.}
  \item{criterion}{Criterion used.}
  \item{combined_results}{Data frame with results for all three models
    when \code{combined = TRUE}; \code{NULL} otherwise.}
  \item{call}{The matched call.}
}
\description{
Performs a residual-based cointegration test designed for very small
time-series samples following Trinh (2022). The procedure first regresses
the dependent variable on the regressors (with optional structural-break
dummies), then applies an ADF-type unit-root test to the residuals. Break
dates are selected by minimising the residual ADF statistic (\code{"adf"})
or the sum of squared residuals (\code{"ssr"}).
}
\examples{
set.seed(1)
y <- cumsum(rnorm(20))
x <- cumsum(rnorm(20))
res <- cointsmall(y = y, x = x, breaks = 0, model = "o",
                  level = 0.05)
print(res)
}
\references{
Trinh, J. (2022). Testing for cointegration with structural changes in
very small sample. THEMA Working Paper n\ifelse{html}{\out{&deg;}}{\eqn{^\circ}}2022-01,
CY Cergy Paris Universite.
\url{https://www.thema.u-cergy.fr/IMG/pdf/2022-01.pdf}
}
