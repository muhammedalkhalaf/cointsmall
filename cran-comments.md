## cointsmall 1.0.4

This release corrects the critical values; the 1.0.3 submission (DOI metadata only) should be discarded in favour of this one.

* Bug fix: the critical values did not come from Trinh (2022). The model "o" table was shifted by one regressor (the m = 1 row held the Dickey-Fuller values of a single series), and the tables for one and two breaks were not the published response surfaces. The 5% critical values are now computed from the response surfaces in Table 13 of Trinh (2022), for m = 1, 2, 3 regressors; they reproduce Table 1 of the paper.
* Bug fix: the p-value was extrapolated from these tables and could exceed 1 (16.29 in one test case). Trinh (2022) publishes the 5% quantile only, so `pvalue`, `cv01` and `cv10` are now `NA`, `level` must be 5, and `cointsmall_cv()` stops for other levels or more than three regressors.
* The author of the methodology is Jerome Trinh; the reference is corrected, and the contributor entry in Authors@R, which named another person and did not correspond to any contributed code, is removed.
* The ADF* statistic and break dates are unchanged; they agree with the Stata command cointsmall (SSC) on the same data.

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
