

test_that("get_GAMs returns correct nested structure", {

  set.seed(1)
  # Simulate some data for demonstration
  mp = data.frame(matrix(rnorm(500), ncol=10))
  mc = data.frame(matrix(rnorm(1000), ncol=10))
  rc = data.frame(matrix(rnorm(1000), ncol=10))
  colnames(mp) = colnames(mc) = colnames(rc) = paste0(rep(c("tas.", "pr."), each = 5), c(1:5))

  time_c = as.Date("2000-01-01") + 0:99
  time_p = as.Date("2020-01-01") + 0:49

  # Simulate 5 locations
  locs = data.frame(Lon = runif(5, -180, 180), Lat = runif(5, -90, 90), Id = 1:5)

  var_names = c("tas", "pr")

  families = list("tas" = gaussian(), "pr" = gaussian())

  res = get_GAMs(
    mp = mp,
    mc = mc,
    rc = rc,
    time_p = time_p,
    time_c = time_c,
    locs = locs,
    var_names = var_names,
    families = families,
    cores = 5
  )

  # ---- top-level structure ----
  expect_type(res, "list")
  expect_setequal(names(res), c("mc", "mp", "rc"))

  # ---- second level ----
  for (ds in names(res)) {
    expect_setequal(names(res[[ds]]), var_names)
  }

  # ---- objects returned ----
  expect_s3_class(res$mc$tas, "bam")
  expect_s3_class(res$rc$pr,  "bam")

  # ------ constructs correct GAM formulas ------
  expect_equal(
    deparse(res$mc$tas$formula),
    "tas ~ te(t, Lat, Lon, bs = c(\"cc\", \"tp\", \"tp\"))"
  )
  expect_equal(
    deparse(res$mp$pr$formula),
    "pr ~ te(t, Lat, Lon, bs = c(\"cc\", \"tp\", \"tp\"))"
  )


  expect_silent(
    get_GAMs(
      mp, mc, rc,
      locs,
      time_c, time_p,
      var_names,
      families,
      cores = 1
    )
  )

})
