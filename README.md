# cointsmall

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen.svg)](https://cran.r-project.org/)
[![CRAN status](https://www.r-pkg.org/badges/version/cointsmall)](https://CRAN.R-project.org/package=cointsmall)
<!-- badges: end -->

## Overview

The `cointsmall` package implements cointegration tests with structural breaks designed specifically for small sample sizes. It follows the methodology of Trinh (2022), which provides small-sample adjusted critical values for testing the null hypothesis of no cointegration.

## Installation

```r
# Install from CRAN (when available)
install.packages("cointsmall")

# Or install the development version from GitHub
# devtools::install_github("merwanroudane/cointsmall")
```

## Features

- **Small sample focused**: Critical values adjusted for small samples using response surface methodology
- **Structural breaks**: Test for cointegration allowing for 0, 1, or 2 structural breaks
- **Multiple model specifications**:
  - Model "o": No structural break
  - Model "c": Break in constant (level shift)
  - Model "cs": Break in constant and slope (regime change)
- **Endogenous break detection**: Break dates determined by minimizing ADF statistic or SSR
- **Combined testing**: Evaluate all models simultaneously with automatic model selection

## Usage

### Basic Cointegration Test

```r
library(cointsmall)

# Generate cointegrated series
set.seed(42)
n <- 50
e <- cumsum(rnorm(n))  # Common stochastic trend
y <- 2 + 3 * e + rnorm(n, sd = 0.5)
x <- e + rnorm(n, sd = 0.3)

# Test with no structural break
result <- cointsmall(y, x, breaks = 0)
print(result)
```

### Testing with Structural Breaks

```r
# Generate series with structural break
y_break <- c(2 + 2 * e[1:25], 5 + 4 * e[26:50]) + rnorm(n, sd = 0.3)

# Test with one break (break in constant and slope)
result <- cointsmall(y_break, x, breaks = 1, model = "cs")
print(result)
summary(result)
```

### Combined Testing Procedure

```r
# Test all model specifications
combined <- cointsmall_combined(y_break, x, breaks = 1)
print(combined)
```

### Critical Values

```r
# Get critical values for specific parameters
cv <- cointsmall_cv(T = 50, m = 1, breaks = 1, model = "cs")
print(cv)
```

## Models

| Model | Description | Cointegrating Regression |
|-------|-------------|-------------------------|
| o | No break | y = α + β'x + ε |
| c | Break in constant | y = α₁ + α₂D + β'x + ε |
| cs | Break in constant and slope | y = α₁ + α₂D + β₁'x + β₂'(x·D) + ε |

Where D is a dummy variable equal to 1 after the break date.

## References

- Trinh, H. H. (2022). Testing for cointegration with structural changes in very small sample. THEMA Working Paper n°2022-01, CY Cergy Paris Université. [doi:10.2139/ssrn.4037745](https://doi.org/10.2139/ssrn.4037745)

- Gregory, A. W., & Hansen, B. E. (1996). Residual-based tests for cointegration in models with regime shifts. *Journal of Econometrics*, 70(1), 99-126. [doi:10.1016/0304-4076(69)41685-7](https://doi.org/10.1016/0304-4076(69)41685-7)

- Hatemi-J, A. (2008). Tests for cointegration with two unknown regime shifts with an application to financial market integration. *Empirical Economics*, 35(3), 497-505. [doi:10.1007/s00181-007-0175-9](https://doi.org/10.1007/s00181-007-0175-9)

- Engle, R. F., & Granger, C. W. J. (1987). Co-integration and error correction: Representation, estimation, and testing. *Econometrica*, 55(2), 251-276. [doi:10.2307/1913236](https://doi.org/10.2307/1913236)

## Author

Dr. Merwan Roudane (merwanroudane920@gmail.com)

## License

GPL (>= 3)
