# cointsmall

**Cointegration Test for Very Small Samples**

The `cointsmall` R package implements the cointegration test of Trinh (2022) designed specifically for time-series data with very few observations (as few as 12). The procedure allows for up to two structural breaks in the long-run relationship.

## Features

- Residual-based ADF cointegration test
- Up to two structural breaks (models: o, c, cs)
- Automatic BIC lag selection
- Combined sequential model-selection procedure
- Simulation-based critical values for small samples

## Installation

```r
install.packages("cointsmall")
```

## Usage

```r
library(cointsmall)

set.seed(1)
y <- cumsum(rnorm(20))
x <- cumsum(rnorm(20))

# No structural break
res0 <- cointsmall(y, x, breaks = 0, model = "o")
print(res0)

# One break, break in constant and slope
res1 <- cointsmall(y, x, breaks = 1, model = "cs")
print(res1)

# Combined procedure
res_c <- cointsmall(y, x, breaks = 1, combined = TRUE)
print(res_c)
```

## Reference

Trinh, J. (2022). Testing for cointegration with structural changes in very small sample. THEMA Working Paper n°2022-01, CY Cergy Paris Université. https://www.thema.u-cergy.fr/IMG/pdf/2022-01.pdf

## License

GPL-3
