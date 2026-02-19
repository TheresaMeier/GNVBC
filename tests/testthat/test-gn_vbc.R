test_that("gn_vbc works", {
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

  # Simulate 10 locations
  locs <- data.frame(
    Id = 1:n_locs,
    Lon = runif(n_locs, -180, 180),
    Lat = runif(n_locs, -90, 90)
  )

  out <- gn_vbc(
    mp = mp,
    mc = mc,
    rc = rc,
    time_p = time_p,
    time_c = time_c,
    locs = locs,
    nrows = 5,
    ncols = 2,
    families = list("tas" = gaussian(), "pr" = Gamma(link = "log")),
    fixed = FALSE,
    mask = TRUE,
    bridge_var = 3,
    trunc_lvl = NULL,
    seed = 1,
    cores = 5
  )

  expect_type(out, "list")
  expect_named(out, c("corrected_mp", "rvine_mp", "rvine_rc"))

  expect_s3_class(out$rvine_mp, "vinecop")
  expect_s3_class(out$rvine_rc, "vinecop")

  expect_true(is.data.frame(out$corrected_mp))

  expect_equal(
    dim(out$corrected_mp),
    dim(mp)
  )

  expect_equal(
    colnames(out$corrected_mp),
    colnames(mp)
  )
  expect_false(anyNA(out$corrected_mp))

  # seed ensures reproducibility

  out1 <- gn_vbc(
    mp = mp,
    mc = mc,
    rc = rc,
    time_p = time_p,
    time_c = time_c,
    locs = locs,
    nrows = 5,
    ncols = 2,
    families = list("tas" = gaussian(), "pr" = Gamma(link = "log")),
    fixed = FALSE,
    mask = TRUE,
    bridge_var = 3,
    trunc_lvl = NULL,
    seed = 12,
    cores = 5
  )
  out2 <- gn_vbc(
    mp = mp,
    mc = mc,
    rc = rc,
    time_p = time_p,
    time_c = time_c,
    locs = locs,
    nrows = 5,
    ncols = 2,
    families = list("tas" = gaussian(), "pr" = Gamma(link = "log")),
    fixed = FALSE,
    mask = TRUE,
    bridge_var = 3,
    trunc_lvl = NULL,
    seed = 12,
    cores = 5
  )

  expect_equal(out1$corrected_mp, out2$corrected_mp)

  # truncation level is respected
  out <- gn_vbc(
    mp = mp,
    mc = mc,
    rc = rc,
    time_p = time_p,
    time_c = time_c,
    locs = locs,
    nrows = 5,
    ncols = 2,
    families = list("tas" = gaussian(), "pr" = Gamma(link = "log")),
    fixed = FALSE,
    mask = TRUE,
    bridge_var = 3,
    trunc_lvl = 2,
    seed = 12,
    cores = 5
  )

  expect_equal(
    out$rvine_mp$structure$trunc_lvl,
    2
  )
})
