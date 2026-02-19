#' @title Back-transformation of bias-corrected PIT values to the original scale
#' @description
#' Transforms bias-corrected Probability Integral Transform (PIT) values
#' back to the original data scale using the fitted marginal GAMs.
#'
#' For each variable:
#' - The corrected PIT values are mapped to empirical reference PITs.
#' - The inverse CDF of the fitted marginal distribution is applied.
#' - The mean structure is adjusted using:
#'   \deqn{\mu_{rc} + \mu_{mp} - \mu_{mc}}
#'
#' This step reconstructs bias-corrected projections on the physical scale.
#'
#' @param gam_fit A list containing fitted GAM objects for `rc`, `mc` and `mp`
#'
#' @param u_mph_wide Data frame of bias-corrected PIT values in wide format,
#'   including columns `t`, `Lon`, `Lat`, and variable-location columns.
#'
#' @param u_mph Long-format version of bias-corrected PIT values.
#'
#' @param components_rc Output of `extract_components()` for calibration data,
#'   containing seasonality and remainder components.
#'
#' @return
#' A wide-format data frame of reconstructed, bias-corrected values
#' on the original physical scale.
#'
#' @details
#' The inverse transformation depends on the marginal family:
#' Gaussian, Gamma, Tweedie, Beta, Poisson, or Binomial.
#'
#' For Tweedie and Beta models, the additional parameters
#' (power or precision) are extracted from the fitted GAM object.
#'
inverse_PITs <- function(gam_fit, u_mph_wide, u_mph, components_rc) {
  vars_unique <- names(gam_fit$rc)

  # Initialize corrected data container
  x_mph_wide <- u_mph_wide

  for (var in vars_unique) {
    # Retrieve GAM fits
    fit <- gam_fit$rc[[var]]

    # Predict mean components
    mu_rc <- mgcv::predict.gam(
      fit,
      newdata = data.frame(u_mph_wide[, c("t", "Lon", "Lat")]),
      type = "response"
    )
    mu_mc <- mgcv::predict.gam(
      gam_fit$mc[[var]],
      newdata = data.frame(u_mph_wide[, c("t", "Lon", "Lat")]),
      type = "response"
    )
    mu_mp <- mgcv::predict.gam(
      gam_fit$mp[[var]],
      newdata = data.frame(u_mph_wide[, c("t", "Lon", "Lat")]),
      type = "response"
    )

    # Identify marginal distribution
    family_name <- fit$family$family

    # Extract Tweedie power parameter if needed
    if (grepl("Tweedie", family_name)) {
      p_rc <- fit$family$getTheta(TRUE)
    } else {
      p_rc <- NULL
    }

    # Extract dispersion / precision parameter
    if (grepl("Beta", family_name)) {
      phi_rc <- fit$family$getTheta(TRUE)
    } else {
      phi_rc <- fit$sig2
    }

    # Identify location-specific PIT columns
    cols <- grep(
      paste0("^", var, "\\."),
      names(components_rc$remainder_orig),
      value = TRUE
    )

    # Compute empirical quantiles of reference PITs
    q_df <- sapply(
      cols,
      function(col) {
        stats::quantile(
          x = components_rc$remainder_orig[[col]],
          probs = data.frame(u_mph)[col]
        )
      },
      simplify = "data.frame"
    )

    p_var <- c(t(q_df))

    # Inverse PIT mapping with mean adjustment
    x_mph_wide[[var]] <- inverse_cdf_from_pit(
      p = p_var,
      mu = mu_rc + mu_mp - mu_mc,
      phi = phi_rc,
      family_name = family_name,
      power = p_rc
    )
  }
  x_mph_wide
}
