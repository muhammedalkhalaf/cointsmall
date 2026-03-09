test_that("cointsmall works with no breaks", {
  set.seed(42)
  n <- 50
  e <- cumsum(rnorm(n))
  y <- 2 + 3 * e + rnorm(n, sd = 0.5)
  x <- e + rnorm(n, sd = 0.3)
  
  result <- cointsmall(y, x, breaks = 0)
  
  expect_s3_class(result, "cointsmall")
  expect_equal(result$breaks, 0)
  expect_equal(result$model, "o")
  expect_true(is.numeric(result$statistic))
  expect_true(is.numeric(result$cv))
  expect_true(is.logical(result$reject))
  expect_true(length(result$residuals) == n)
})

test_that("cointsmall works with one break", {
  set.seed(42)
  n <- 50
  e <- cumsum(rnorm(n))
  x <- e + rnorm(n, sd = 0.3)
  # Create series with structural break
  y <- c(2 + 2 * e[1:25], 5 + 4 * e[26:50]) + rnorm(n, sd = 0.3)
  
  # Test model c
  result_c <- cointsmall(y, x, breaks = 1, model = "c")
  expect_s3_class(result_c, "cointsmall")
  expect_equal(result_c$breaks, 1)
  expect_equal(result_c$model, "c")
  expect_true(length(result_c$break_dates) == 1)
  
  # Test model cs
  result_cs <- cointsmall(y, x, breaks = 1, model = "cs")
  expect_s3_class(result_cs, "cointsmall")
  expect_equal(result_cs$model, "cs")
  expect_true(length(result_cs$break_dates) == 1)
})

test_that("cointsmall works with two breaks", {
  set.seed(42)
  n <- 60
  e <- cumsum(rnorm(n))
  x <- e + rnorm(n, sd = 0.3)
  # Create series with two structural breaks
  y <- c(2 + 2 * e[1:20], 5 + 4 * e[21:40], 3 + 3 * e[41:60]) + rnorm(n, sd = 0.3)
  
  result <- cointsmall(y, x, breaks = 2, model = "cs")
  
  expect_s3_class(result, "cointsmall")
  expect_equal(result$breaks, 2)
  expect_true(length(result$break_dates) == 2)
  expect_true(result$break_dates[1] < result$break_dates[2])
})

test_that("cointsmall_combined works", {
  set.seed(42)
  n <- 50
  e <- cumsum(rnorm(n))
  x <- e + rnorm(n, sd = 0.3)
  y <- c(2 + 2 * e[1:25], 5 + 4 * e[26:50]) + rnorm(n, sd = 0.3)
  
  result <- cointsmall_combined(y, x, breaks = 1)
  
  expect_s3_class(result, "cointsmall_combined")
  expect_true("o" %in% names(result$results))
  expect_true("c" %in% names(result$results))
  expect_true("cs" %in% names(result$results))
  expect_true(result$selected_model %in% c("none", "o", "c", "cs"))
})

test_that("critical values are correct", {
  # Test that critical values become less negative with sample size
  # (small samples need more stringent critical values)
  cv_small <- cointsmall_cv(T = 20, m = 1, breaks = 0, model = "o")
  cv_large <- cointsmall_cv(T = 100, m = 1, breaks = 0, model = "o")
  
  # Small sample critical values should be more negative (more stringent)
  expect_true(cv_small$cv05 < cv_large$cv05)
  
  # Test ordering: cv01 < cv05 < cv10
  cv <- cointsmall_cv(T = 50, m = 1, breaks = 0, model = "o")
  expect_true(cv$cv01 < cv$cv05)
  expect_true(cv$cv05 < cv$cv10)
})

test_that("input validation works", {
  # Use 20 observations to avoid sample size errors
  x <- 1:20
  y <- 1:20
  
  # Test breaks validation
  expect_error(cointsmall(y, x, breaks = 3), "'breaks' must be 0, 1, or 2")
  
  # Test model validation
  expect_error(cointsmall(y, x, breaks = 0, model = "c"), 
               "models 'c' and 'cs' require breaks > 0")
  
  # Test sample size
  expect_error(cointsmall(1:5, matrix(1:5, ncol = 1), breaks = 0), 
               "Sample size too small")
  
  # Test dimension mismatch
  expect_error(cointsmall(1:20, matrix(1:15, ncol = 1), breaks = 0),
               "Length of 'y' must equal number of rows in 'x'")
})

test_that("ADF test selects reasonable lags", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 2 * x + rnorm(n)
  
  result <- cointsmall(y, matrix(x, ncol = 1), breaks = 0)
  
  # Lags should be reasonable
  expect_true(result$lags >= 0)
  expect_true(result$lags <= floor(12 * (n / 100)^0.25))
})

test_that("print and summary methods work", {
  set.seed(42)
  n <- 30
  e <- cumsum(rnorm(n))
  y <- 2 + 3 * e + rnorm(n, sd = 0.5)
  x <- e + rnorm(n, sd = 0.3)
  
  result <- cointsmall(y, x, breaks = 0)
  
  # These should not error
  expect_output(print(result))
  expect_output(summary(result))
  
  # Combined
  result_comb <- cointsmall_combined(y, x, breaks = 1)
  expect_output(print(result_comb))
  expect_output(summary(result_comb))
})

test_that("multiple regressors work", {
  set.seed(42)
  n <- 50
  e <- cumsum(rnorm(n))
  x1 <- e + rnorm(n, sd = 0.3)
  x2 <- 0.5 * e + rnorm(n, sd = 0.2)
  x <- cbind(x1, x2)
  y <- 1 + 2 * e + rnorm(n, sd = 0.5)
  
  result <- cointsmall(y, x, breaks = 0)
  
  expect_s3_class(result, "cointsmall")
  expect_equal(result$nvar, 2)
  expect_true(length(result$coefficients) == 3)  # intercept + 2 slopes
})

test_that("criterion selection works", {
  set.seed(42)
  n <- 50
  e <- cumsum(rnorm(n))
  x <- e + rnorm(n, sd = 0.3)
  y <- c(2 + 2 * e[1:25], 5 + 4 * e[26:50]) + rnorm(n, sd = 0.3)
  
  result_adf <- cointsmall(y, x, breaks = 1, criterion = "adf")
  result_ssr <- cointsmall(y, x, breaks = 1, criterion = "ssr")
  
  expect_equal(result_adf$criterion, "adf")
  expect_equal(result_ssr$criterion, "ssr")
  
  # Both should find a break date
  expect_true(!is.na(result_adf$break_dates))
  expect_true(!is.na(result_ssr$break_dates))
})
