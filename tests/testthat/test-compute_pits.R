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
  fit_gamma <- mgcv::gam(y ~ s(x0), data = dat_gamma, family = Gamma(link = "log"))
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

