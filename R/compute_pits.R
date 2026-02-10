#' @title Compute Probability Integral Transform (PIT) values from a fitted GAM
#' @description Transforms observed responses into PIT values using the conditional
#' distribution implied by a fitted GAM. The distributional form is
#' selected based on the variable name prefix.
#'
#' This function is typically used to extract model residuals on the
#' copula scale (uniform on `[0,1]`) for dependence modeling.
#' @param y Numeric vector of observed responses
#' @param fit Fitted GAM object (mgcv), providing mean and dispersion estimates
#'
#' @returns Numeric vector of PIT values in `[0,1]`
#' @export
#'
#' @details The function supports the following families: Gaussian, Gamma, Tweedie, Beta, Poisson, Binomial.
#' The family is inferred from the fitted model's family component. For discrete families (Poisson, Binomial), randomization is applied to ensure PIT values are continuous on `[0,1]`. For continuous families, the appropriate CDF is used directly.
#'
#' @examples
#'
#' # Gaussian
#' set.seed(2)
#' dat <- data.frame(y = rnorm(100), x0 = runif(100))
#' b <- mgcv::gam(y ~ s(x0), data = dat)
#' pits_norm <- compute_pit(dat$y, b)
#'
#' # Tweedie
#' set.seed(2)
#' dat <- data.frame(
#'   y = tweedie::rtweedie(200, mu = 1, phi = 1, power = 1.5),
#'   x0 = runif(200)
#' )
#' b <- mgcv::gam(y ~ s(x0), data = dat, family = mgcv::tw(link = "log"))
#' pits_tw <- compute_pit(dat$y, b)
#'
#' # Beta
#' set.seed(2)
#' # simulate some data...
#' dat <- data.frame(matrix(rbeta(200, 1, 1), ncol = 2))
#' colnames(dat) = c("y", "x0")
#' b <- mgcv::gam(y~s(x0),data=dat, family = betar())
#' pits_beta = compute_pit(dat$y, b)
#'
#' # Poisson
#' set.seed(2)
#' # simulate some data...
#' dat <- data.frame(matrix(rpois(200, 5), ncol = 2))
#' colnames(dat) = c("y", "x0")
#' b <- mgcv::gam(y~s(x0),data=dat, family = poisson())
#' pits_pois = compute_pit(dat$y, b)
#'
#' # Binomial
#' set.seed(2)
#' dat <- data.frame(y  = rbinom(100, 1, 0.3), x0 = stats::runif(100))
#' b <- mgcv::gam(y ~ s(x0), data = dat, family = binomial())
#' pits_binom = compute_pit(dat$y, b)
#'


compute_pit <- function(y, fit) {

  family_name <- fit$family$family
  if (grepl("Tweedie", family_name)) {
    family_name <- "Tweedie"
  } else if (grepl("Beta", family_name)) {
    family_name <- "betar"
  }
  mu_hat      <- mgcv::predict.gam(fit, type = "response")
  phi_hat     <- fit$sig2

  u <- switch(
    family_name,

    # --------------------------------------------------
    # Gaussian
    # --------------------------------------------------
    gaussian = {
      sigma_hat <- sqrt(phi_hat)
      stats::pnorm(y, mean = mu_hat, sd = sigma_hat)
    },

    # --------------------------------------------------
    # Gamma
    # --------------------------------------------------
    Gamma = {
      shape <- 1 / phi_hat
      scale <- phi_hat * mu_hat
      stats::pgamma(y, shape = shape, scale = scale)
    },

    # --------------------------------------------------
    # Tweedie
    # --------------------------------------------------
    Tweedie = {
      power_hat <- fit$family$getTheta(TRUE)

      u <- tweedie::ptweedie(
        q     = y,
        mu    = mu_hat,
        phi   = phi_hat,
        power = power_hat
      )

      zero_idx <- which(y == 0)
      if (length(zero_idx) > 0) {
        u[zero_idx] <- stats::runif(length(zero_idx), 0, u[zero_idx])
      }

      u
    },

    # --------------------------------------------------
    # Beta regression
    # --------------------------------------------------
    betar = {
      phi_hat <- fit$family$getTheta(TRUE)

      stats::pbeta(
        y,
        shape1 = mu_hat * phi_hat,
        shape2 = (1 - mu_hat) * phi_hat
      )
    },

    # --------------------------------------------------
    # Poisson
    # --------------------------------------------------
    poisson = {
      stats::ppois(y - 1, lambda = mu_hat) +
        stats::runif(length(y)) * stats::dpois(y, lambda = mu_hat)
    },

    # --------------------------------------------------
    # Binomial (Bernoulli)
    # --------------------------------------------------
    binomial = {
      stats::pbinom(y - 1, size = 1, prob = mu_hat) +
        stats::runif(length(y)) * stats::dbinom(y, size = 1, prob = mu_hat)
    },

    # --------------------------------------------------
    # Fallback
    # --------------------------------------------------
    stop(
      "Unsupported family: ", family_name,
      call. = FALSE
    )
  )

  return(u)
}
