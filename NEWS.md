# cointsmall 1.0.4

* Bug fix: the critical values did not come from Trinh (2022). The model "o" table was shifted by one regressor (the m = 1 row held the Dickey-Fuller values of a single series), and the tables for one and two breaks were not the published response surfaces. The 5% critical values are now computed from the response surfaces in Table 13 of Trinh (2022), for m = 1, 2, 3 regressors; they reproduce Table 1 of the paper.
* Bug fix: the p-value was extrapolated from these tables and could exceed 1 (16.29 in one test case). Trinh (2022) publishes the 5% quantile only, so `pvalue`, `cv01` and `cv10` are now `NA`, `level` must be 5, and `cointsmall_cv()` stops for other levels or more than three regressors.
* The author of the methodology is Jerome Trinh; the reference is corrected, and the contributor entry in Authors@R, which named another person and did not correspond to any contributed code, is removed.
* The ADF* statistic and break dates are unchanged; they agree with the Stata command cointsmall (SSC) on the same data.

# cointsmall 1.0.3

* Removed a DOI wrongly attached to the MacKinnon (2010) working paper; the citation text is unchanged. No changes to code.

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
