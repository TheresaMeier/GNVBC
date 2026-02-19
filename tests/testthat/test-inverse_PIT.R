test_that("inverse_PITs returns correct structure", {
  set.seed(1)
  n_locs <- 10
  # Simulate some data for demonstration
  mp <- data.frame(
    cbind(matrix(rnorm(500), ncol = n_locs)),
    matrix(rgamma(500, shape = 2), ncol = n_locs)
  )
  mc <- data.frame(
    cbind(matrix(rnorm(1000), ncol = n_locs)),
    matrix(rgamma(1000, shape = 2), ncol = n_locs)
  )
  rc <- data.frame(
    cbind(matrix(rnorm(1000), ncol = n_locs)),
    matrix(rgamma(1000, shape = 2), ncol = n_locs)
  )

  colnames(mp) <- colnames(mc) <- colnames(rc) <- paste0(
    rep(c("tas.", "pr."), each = n_locs),
    c(1:n_locs)
  )

  time_c <- as.Date("2000-01-01") + 0:99
  time_p <- as.Date("2020-01-01") + 0:49

  var_names <- c("tas", "pr")

  families <- list("tas" = gaussian(), "pr" = Gamma(link = "log"))

  locs <- data.frame(
    Id = 1:n_locs,
    Lon = runif(n_locs, -180, 180),
    Lat = runif(n_locs, -90, 90)
  )

  gam_fit <- get_GAMs(
    mp,
    mc,
    rc,
    locs,
    time_c,
    time_p,
    var_names,
    families
  )

  components_rc <- extract_components(gam_fit$rc, rc, locs, time_c)

  u_mph <- rvinecopulib::pseudo_obs(mp)

  u_mph_wide <- transform_to_wide_format(u_mph, locs, var_names, time_p)

  result <- inverse_PITs(
    gam_fit,
    u_mph_wide,
    u_mph,
    components_rc
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), n_locs * length(time_p))

  # No NAs produced
  expect_equal(sum(is.na(result)), 0)

  # Correct columns
  expect_setequal(
    colnames(result),
    c("time", "Id", "tas", "pr", "t", "Lon", "Lat")
  )
  expect_equal(as.Date(result$time), result$time)

  for (col in c("Id", "tas", "pr", "t", "Lon", "Lat")) {
    expect_true(
      is.numeric(result[[col]]),
      info = paste("Column", col, "is not numeric")
    )
  }
})
