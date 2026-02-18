test_that("reordering of data frames works", {
  # Generate random data
  set.seed(42)
  data <- data.frame(matrix(rnorm(400), ncol = 4))
  colnames(data) <- c("var1.loc1", "var1.loc2", "var2.loc1", "var2.loc2")

  # Test correct reordering
  expect_equal(
    reorder_dataset(data, direction = "variable-major"),
    data[, c("var1.loc1", "var1.loc2", "var2.loc1", "var2.loc2")]
  )
  expect_equal(
    reorder_dataset(data, direction = "variable-major", order = c(2, 1)),
    data[, c("var2.loc1", "var2.loc2", "var1.loc1", "var1.loc2")]
  )
  expect_equal(
    reorder_dataset(data, direction = "location-major"),
    data[, c("var1.loc1", "var2.loc1", "var1.loc2", "var2.loc2")]
  )
  expect_equal(
    reorder_dataset(data, direction = "location-major", order = c(2, 1)),
    data[, c("var2.loc1", "var1.loc1", "var2.loc2", "var1.loc2")]
  )

  expect_error(
    reorder_dataset(data, direction = "invalid-direction"),
    "Please select either 'variable-major' or 'location-major' ordering."
  )

  expect_error(reorder_dataset(
    data,
    direction = "variable-major",
    order = c("var1", "var2")
  ))
})


test_that("transform_to_wide_format returns correct structure and values", {
  data <- data.frame(
    var1.1 = c(1, 2, 3),
    var1.2 = c(4, 5, 6),
    var2.1 = c(7, 8, 9),
    var2.2 = c(10, 11, 12)
  )

  time <- as.Date("2000-01-01") + 0:2

  locs <- data.frame(
    Id = c(1, 2),
    Lat = c(10, 20),
    Lon = c(30, 40)
  )

  vars <- c("var1", "var2")

  out <- transform_to_wide_format(data, locs, vars, time)

  # ---- structure ----
  expect_s3_class(out, "data.frame")

  expect_true(all(
    c("time", "Id", "t", "Lat", "Lon", "var1", "var2") %in% names(out)
  ))

  # ---- dimensions ----
  expect_equal(nrow(out), nrow(data) * nrow(locs))

  # ---- types ----
  expect_s3_class(out$time, "Date")
  expect_type(out$Id, "double")

  # ---- temporal covariate ----
  expect_equal(out$t, lubridate::yday(out$time))

  # ---- spatial join ----
  expect_equal(
    out$Lat[out$Id == 1],
    rep(10, nrow(data))
  )

  # ---- reshaping correctness (spot checks) ----
  expect_equal(
    out$var1[out$Id == 1],
    c(1, 2, 3)
  )

  expect_equal(
    out$var2[out$Id == 2],
    c(10, 11, 12)
  )
})
