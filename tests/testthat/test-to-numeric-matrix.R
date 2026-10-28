test_that("to_numeric_matrix converts numeric, logical, and POSIXct columns", {
  df <- data.frame(
    a = c(1.5, 2.5),
    b = c(TRUE, FALSE),
    d = as.POSIXct(c("2020-01-01", "2020-01-02"), tz = "UTC")
  )
  m <- to_numeric_matrix(df, c("a", "b", "d"))
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(2L, 3L))
  expect_false(anyNA(m))
  expect_equal(m[, "a"], c(1.5, 2.5), ignore_attr = TRUE)
  expect_equal(m[, 2], c(1, 0), ignore_attr = TRUE)
})

test_that("to_numeric_matrix refuses character columns instead of returning all-NA", {
  df <- data.frame(x = c(1, 2), sex = c("F", "M"), race = c("a", "b"))
  expect_error(to_numeric_matrix(df, c("x", "sex", "race")), "sex, race")
  # Numeric columns of the same frame still convert.
  expect_silent(to_numeric_matrix(df, "x"))
})
