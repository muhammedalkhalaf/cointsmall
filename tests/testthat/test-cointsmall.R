test_that("cointsmall model o works", {
  set.seed(1)
  y <- cumsum(rnorm(20))
  x <- cumsum(rnorm(20))
  res <- cointsmall(y = y, x = x, breaks = 0, model = "o")
  expect_s3_class(res, "cointsmall")
  expect_true(is.numeric(res$adf_stat))
  expect_true(is.logical(res$reject))
  expect_true(res$breaks == 0)
})

test_that("cointsmall model c with 1 break works", {
  set.seed(2)
  y <- cumsum(rnorm(25))
  x <- cumsum(rnorm(25))
  res <- cointsmall(y = y, x = x, breaks = 1, model = "c")
  expect_s3_class(res, "cointsmall")
  expect_length(res$break_obs, 1)
  expect_true(res$break_obs >= 1 && res$break_obs <= 25)
})

test_that("cointsmall model cs with 2 breaks works", {
  set.seed(3)
  y <- cumsum(rnorm(30))
  x <- cumsum(rnorm(30))
  res <- cointsmall(y = y, x = x, breaks = 2, model = "cs")
  expect_s3_class(res, "cointsmall")
  expect_length(res$break_obs, 2)
  expect_true(res$break_obs[1] < res$break_obs[2])
})

test_that("cointsmall combined procedure works", {
  set.seed(4)
  y <- cumsum(rnorm(25))
  x <- cumsum(rnorm(25))
  res <- cointsmall(y = y, x = x, breaks = 1, combined = TRUE)
  expect_false(is.null(res$combined_results))
  expect_true(nrow(res$combined_results) == 3)
})

test_that("print.cointsmall does not error", {
  set.seed(5)
  y <- cumsum(rnorm(20))
  x <- cumsum(rnorm(20))
  res <- cointsmall(y = y, x = x, breaks = 0, model = "o")
  expect_output(print(res))
})

test_that("error for too small sample", {
  expect_error(
    cointsmall(y = rnorm(10), x = rnorm(10), breaks = 0, model = "o"),
    "too small"
  )
})
