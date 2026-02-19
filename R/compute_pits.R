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
  mu_hat <- mgcv::predict.gam(fit, type = "response")
  phi_hat <- fit$sig2

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
        q = y,
        mu = mu_hat,
        phi = phi_hat,
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
      "Unsupported family: ",
      family_name,
      call. = FALSE
    )
  )

  return(u)
}


#' @title Extract seasonality and remainder components from fitted GAMs
#' @description Decomposes multivariate time series data into:
#' - Seasonality component: the fitted mean from the GAMs, capturing temporal and spatial patterns (GAM predictions)
#' - Remainder component: the Probability Integral Transform (PIT) values of the original data relative to the fitted GAMs, representing residuals on the copula scale.
#'
#' This function transforms the input data into a long format, applies the fitted GAMs to extract the mean and compute PIT values for each variable, and then reshapes the results back into a wide format suitable for dependence modeling.
#' @param gam_list Named list of fitted GAM objects (one per variable); obtained from `get_GAMs()`.
#' @param data Original data frame containing the multivariate time series, with columns corresponding to variables and spatial locations.
#' @param locs Data frame of spatial locations with columns: Id, Lat, Lon
#' @param time Date vector corresponding to rows of the input data frame, used to align with GAM predictions.
#'
#' @returns A list containing:
#' - `seasonality`: Wide-format data frame of fitted mean values from the GAMs
#' - `remainder`: Wide-format data frame of pseudo-observations of PIT values for dependence modeling
#' - `remainder_orig`: Wide-format data frame of original PIT values before transformation to pseudo-observations
#' @export
#'
#' @examples
#' set.seed(1)
#' n_locs = 10
#' # Simulate some data for demonstration
#' mp = data.frame(cbind(matrix(rnorm(500), ncol=n_locs)),
#'                       matrix(rgamma(500, shape = 2), ncol=n_locs))
#' mc = data.frame(cbind(matrix(rnorm(1000), ncol=n_locs)),
#'                       matrix(rgamma(1000, shape = 2), ncol=n_locs))
#' rc = data.frame(cbind(matrix(rnorm(1000), ncol=n_locs)),
#'                       matrix(rgamma(1000, shape = 2), ncol=n_locs))
#'
#' colnames(mp) = colnames(mc) = colnames(rc) = paste0(rep(c("tas.", "pr."),
#'                                                 each = n_locs), c(1:n_locs))
#'
#' time_c = as.Date("2000-01-01") + 0:99
#' time_p = as.Date("2020-01-01") + 0:49
#'
#' # Simulate 10 locations
#' locs <- data.frame(
#' Id = 1:n_locs,
#' Lon = runif(n_locs, -180, 180),
#' Lat = runif(n_locs, -90, 90),
#' Altitude = runif(n_locs, 0, 3000)
#' )
#' var_names = c("tas", "pr")
#'
#' families = list("tas" = gaussian(), "pr" = Gamma(link = "log"))
#'
#' fit = get_GAMs(
#' mp = mp,
#' mc = mc,
#' rc = rc,
#' time_p = time_p,
#' time_c = time_c,
#' locs = locs,
#' var_names = var_names,
#' families = families,
#' cores = 5
#' )
#'
#' comp_mc = extract_components(fit$mc, mc, locs, time_c)
#'
extract_components <- function(gam_list, data, locs, time) {
  stopifnot(length(time) == nrow(data))

  # Variable names inferred from GAM list
  vars <- names(gam_list)

  # Convert data into long format with explicit Id and time columns
  data_long = transform_to_wide_format(data, locs, vars, time)

  # Unique spatial/location identifiers
  id_locs = unique(data_long$Id)

  # Number of observations in long format
  n <- nrow(data_long)

  # Initialize data frames for components
  seasonality_df <- data.frame(matrix(NA, nrow = n, ncol = length(vars)))
  colnames(seasonality_df) <- vars
  remainder_orig_df = seasonality_df

  # Loop over variables and extract components
  for (var in vars) {
    # Corresponding GAM fit
    fit <- gam_list[[var]]

    # Extract fitted mean (seasonal component)
    mu_hat <- mgcv::predict.gam(fit, type = "response")
    seasonality_df[[var]] <- mu_hat

    # Compute PIT values (remainder on copula scale)
    pits = compute_pit(data_long[[var]], fit)
    remainder_orig_df[[var]] <- pits
  }

  # Convert remainder PITs to wide format (variable.location)
  remainder_orig_wide = remainder_orig_df %>%
    dplyr::mutate(
      Id = data_long$Id,
      time = data_long$time
    ) %>%
    tidyr::pivot_wider(
      id_cols = "time",
      names_from = "Id",
      values_from = tidyselect::all_of(vars),
      names_sep = "."
    ) %>%
    dplyr::select(-"time")

  # Convert seasonality component to wide format
  seasonality_wide = seasonality_df %>%
    dplyr::mutate(
      Id = data_long$Id,
      time = data_long$time
    ) %>%
    tidyr::pivot_wider(
      id_cols = "time",
      names_from = "Id",
      values_from = tidyselect::all_of(vars),
      names_sep = "."
    ) %>%
    dplyr::select(-"time")

  return(list(
    seasonality = seasonality_wide,
    remainder = rvinecopulib::pseudo_obs(remainder_orig_wide),
    remainder_orig = remainder_orig_wide
  ))
}

#' @title Inverse CDF from PIT values for various distributions
#' @description Maps PIT values back to the original data scale using the inverse
#' cumulative distribution function (quantile function) of the specified
#' distribution family.
#'
#' @param p Numeric vector of PIT values in `[0, 1]`
#' @param mu Numeric vector of mean parameters from the fitted GAM
#' @param phi Numeric vector of dispersion parameters from the fitted GAM
#' @param family_name Character string indicating the distribution family
#' (e.g., "gaussian", "Gamma", "Tweedie", "Beta", "Poisson", "Binomial")
#' @param power Numeric value of the power parameter for the Tweedie distribution (required if family is Tweedie)
#' @param size Numeric value of the size parameter for the Binomial distribution (required if family is Binomial)
#'
#' @returns Numeric vector of values on the original data scale corresponding to the input PIT values
#' @export
#'
#' @examples
#' p = runif(100)
#' phi = 4
#'
#' # Gaussian
#' mu = rnorm(100, mean = 5, sd = 2)
#' inv_gaussian = inverse_cdf_from_pit(p, mu, phi, "gaussian")
#'
#' # Gamma
#' mu = rgamma(100, shape = 2, rate = 0.5)
#' inv_gamma = inverse_cdf_from_pit(p, mu, phi, "Gamma")
#'
#' # Tweedie
#' mu = runif(100, 0.1, 10)
#' inv_tweedie = inverse_cdf_from_pit(p, mu, phi, "Tweedie", power = 1.5)
#'
#' # Beta
#' mu = runif(100, 0.1, 0.9)
#' inv_beta = inverse_cdf_from_pit(p, mu, phi, "Beta")
#'
#'# Poisson
#' mu = rpois(100, lambda = 5)
#' inv_poisson = inverse_cdf_from_pit(p, mu, phi, "Poisson")
#'
#' # Binomial
#' mu = rbinom(100, size = 1, prob = 0.3)
#' inv_binomial = inverse_cdf_from_pit(p, mu, phi, "Binomial", size = 1)

inverse_cdf_from_pit <- function(
  p,
  mu,
  phi,
  family_name,
  power = NULL,
  size = 1
) {
  # safety
  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("p must be in [0, 1]")
  }

  # --- Gaussian distribution ---
  if (grepl("^gaussian", family_name, ignore.case = TRUE)) {
    return(stats::qnorm(p, mean = mu, sd = sqrt(phi)))

    # --- Gamma distribution ---
  } else if (grepl("^Gamma", family_name, ignore.case = TRUE)) {
    shape <- 1 / phi
    scale <- phi * mu
    return(stats::qgamma(p, shape = shape, scale = scale))

    # --- Inverse Gaussian distribution ---
  } else if (grepl("^inverse.gaussian", family_name, ignore.case = TRUE)) {
    return(statmod::qinvgauss(p, mean = mu, shape = 1 / phi))

    # --- Beta distribution ---
  } else if (grepl("^Beta", family_name, ignore.case = TRUE)) {
    shape1 <- mu * phi
    shape2 <- (1 - mu) * phi
    return(stats::qbeta(p, shape1 = shape1, shape2 = shape2))

    # --- Tweedie distribution ---
  } else if (grepl("^Tweedie", family_name, ignore.case = TRUE)) {
    if (is.null(power)) {
      stop("Tweedie inverse CDF requires 'power'")
    }

    return(tweedie::qtweedie(p, mu = mu, phi = phi, power = power))

    # --- Poisson distribution (discrete) ---
  } else if (grepl("^poisson", family_name, ignore.case = TRUE)) {
    return(stats::qpois(p, lambda = mu))

    # --- Binomial distribution (discrete) ---
  } else if (grepl("^binomial", family_name, ignore.case = TRUE)) {
    return(stats::qbinom(p, size = size, prob = mu))

    # --- Unsupported distribution ---
  } else {
    stop(paste("Inverse CDF not implemented for family:", family_name))
  }
}
