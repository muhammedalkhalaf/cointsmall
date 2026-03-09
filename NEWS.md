# cointsmall 1.0.0

* Initial CRAN release.

## Features

* `cointsmall()`: Main function for cointegration tests with structural breaks in small samples
* `cointsmall_combined()`: Combined testing procedure evaluating all model specifications
* `cointsmall_cv()`: Function to retrieve critical values for different model configurations

## Model Specifications

* Model "o": No structural break (standard Engle-Granger test)
* Model "c": Break in constant only (Gregory-Hansen style level shift)
* Model "cs": Break in constant and slope (regime change model)

## Supported Options

* 0, 1, or 2 structural breaks
* Break date selection via minimum ADF statistic or minimum SSR
* Adjustable trimming parameter for break date search
* Automatic lag selection for ADF test using BIC
* Small-sample adjusted critical values via response surface methodology

## References

* Based on Trinh (2022) "Testing for cointegration with structural changes in very small sample"
