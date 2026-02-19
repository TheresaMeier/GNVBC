test_that("compute_pit returns valid PIT values for all supported families", {
  set.seed(123)

  # ------------------------
  # Gaussian
  # ------------------------
  dat_g <- data.frame(
    y = rnorm(100, mean = 5, sd = 2),
    x0 = runif(100)
  )
  fit_g <- mgcv::gam(y ~ s(x0), data = dat_g, family = gaussian())
  pit_g <- compute_pit(dat_g$y, fit_g)
  expect_type(pit_g, "double")
  expect_true(all(pit_g >= 0 & pit_g <= 1))

  # ------------------------
  # Gamma
  # ------------------------
  dat_gamma <- data.frame(
    y = rgamma(100, shape = 2, rate = 0.5),
    x0 = runif(100)
  )
  fit_gamma <- mgcv::gam(
    y ~ s(x0),
    data = dat_gamma,
    family = Gamma(link = "log")
  )
  pit_gamma <- compute_pit(dat_gamma$y, fit_gamma)
  expect_type(pit_gamma, "double")
  expect_true(all(pit_gamma >= 0 & pit_gamma <= 1))

  # ------------------------
  # Tweedie
  # ------------------------

  dat_tw <- data.frame(
    y = tweedie::rtweedie(100, mu = 1, phi = 1, power = 1.5),
    x0 = runif(100)
  )

  fit_tw <- mgcv::gam(y ~ s(x0), data = dat_tw, family = mgcv::tw(link = "log"))
  pit_tw <- compute_pit(dat_tw$y, fit_tw)

  expect_type(pit_tw, "double")
  expect_true(all(pit_tw >= 0 & pit_tw <= 1))

  # ------------------------
  # Beta
  # ------------------------
  dat_beta <- data.frame(
    y = rbeta(100, 2, 5),
    x0 = runif(100)
  )
  fit_beta <- mgcv::gam(y ~ s(x0), data = dat_beta, family = mgcv::betar())
  pit_beta <- compute_pit(dat_beta$y, fit_beta)
  expect_type(pit_beta, "double")
  expect_true(all(pit_beta >= 0 & pit_beta <= 1))

  # ------------------------
  # Poisson
  # ------------------------
  dat_pois <- data.frame(
    y = rpois(100, lambda = 5),
    x0 = runif(100)
  )
  fit_pois <- mgcv::gam(y ~ s(x0), data = dat_pois, family = poisson())
  pit_pois <- compute_pit(dat_pois$y, fit_pois)
  expect_type(pit_pois, "double")
  expect_true(all(pit_pois >= 0 & pit_pois <= 1))
  expect_false(any(duplicated(pit_pois))) # randomized PIT, no ties

  # ------------------------
  # Binomial
  # ------------------------
  dat_bin <- data.frame(
    y = rbinom(100, size = 1, prob = 0.3),
    x0 = runif(100)
  )
  fit_bin <- mgcv::gam(y ~ s(x0), data = dat_bin, family = binomial())
  pit_bin <- compute_pit(dat_bin$y, fit_bin)
  expect_type(pit_bin, "double")
  expect_true(all(pit_bin >= 0 & pit_bin <= 1))
  expect_false(any(duplicated(pit_bin))) # randomized PIT, no ties

  # ------------------------
  # Unsupported family triggers error
  # ------------------------
  fit_dummy <- mgcv::gam(y ~ s(x0), data = dat_bin, family = quasibinomial())
  expect_error(compute_pit(dat_g$y, fit_dummy))
})

make_test_data <- function(n_locs = 10, n_time = 50) {
  set.seed(42)

  data <- data.frame(
    cbind(matrix(rnorm(n_locs * n_time), ncol = n_locs)),
    matrix(rgamma(500, shape = 2), ncol = n_locs)
  )
  colnames(data) <- paste0(rep(c("tas.", "pr."), each = n_locs), c(1:n_locs))

  time <- as.Date("2000-01-01") + seq_len(n_time) - 1

  locs <- data.frame(
    Id = 1:n_locs,
    Lon = runif(n_locs, -10, 10),
    Lat = runif(n_locs, 40, 50)
  )

  families <- list(
    tas = gaussian(),
    pr = Gamma(link = "log")
  )

  gam_list <- get_GAMs(
    mp = data,
    mc = data,
    rc = data,
    time_p = time,
    time_c = time,
    locs = locs,
    var_names = c("tas", "pr"),
    families = families,
    cores = 1
  )$mc

  list(
    data = data,
    time = time,
    locs = locs,
    gam_list = gam_list
  )
}

test_that("extract_components returns expected structure", {
  x <- make_test_data()

  res <- extract_components(
    gam_list = x$gam_list,
    data = x$data,
    locs = x$locs,
    time = x$time
  )

  expect_type(res, "list")
  expect_named(res, c("seasonality", "remainder", "remainder_orig"))
})

test_that("extract_components returns consistent dimensions", {
  x <- make_test_data()

  res <- extract_components(
    x$gam_list,
    x$data,
    x$locs,
    x$time
  )

  expect_equal(
    dim(res$seasonality),
    dim(res$remainder_orig)
  )

  expect_equal(
    dim(res$remainder),
    dim(res$remainder_orig)
  )
})

test_that("output columns follow var.Id naming convention", {
  x <- make_test_data()

  res <- extract_components(
    x$gam_list,
    x$data,
    x$locs,
    x$time
  )

  expect_true(all(grepl("^(tas|pr)\\.\\d+$", colnames(res$seasonality))))
})

test_that("PIT values are in valid ranges", {
  x <- make_test_data()

  res <- extract_components(
    x$gam_list,
    x$data,
    x$locs,
    x$time
  )

  expect_true(all(res$remainder_orig >= 0))
  expect_true(all(res$remainder_orig <= 1))

  expect_true(all(res$remainder > 0))
  expect_true(all(res$remainder < 1))
})

test_that("extract_components is deterministic for fixed seed", {
  x <- make_test_data()

  res1 <- extract_components(
    x$gam_list,
    x$data,
    x$locs,
    x$time
  )

  res2 <- extract_components(
    x$gam_list,
    x$data,
    x$locs,
    x$time
  )

  expect_equal(res1$seasonality, res2$seasonality)
  expect_equal(res1$remainder_orig, res2$remainder_orig)
})

test_that("wrong time length triggers error", {
  x <- make_test_data()

  bad_time <- x$time[-1]

  expect_error(
    extract_components(x$gam_list, x$data, x$locs, bad_time)
  )
})


test_that("inverse_cdf_from_pit works for all supported families", {
  p <- runif(100)

  expect_type(
    inverse_cdf_from_pit(p, mu = 0, phi = 1, family_name = "gaussian"),
    "double"
  )

  expect_type(
    inverse_cdf_from_pit(p, mu = 2, phi = 0.5, family_name = "Gamma"),
    "double"
  )

  expect_type(
    inverse_cdf_from_pit(p, mu = 0.5, phi = 20, family_name = "Beta"),
    "double"
  )

  expect_type(
    inverse_cdf_from_pit(p, mu = 5, phi = 1, family_name = "poisson"),
    "double"
  )

  expect_type(
    inverse_cdf_from_pit(
      p,
      mu = 0.3,
      phi = NA,
      family_name = "binomial",
      size = 1
    ),
    "double"
  )
})
