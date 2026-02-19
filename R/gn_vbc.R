#' @title Multivariate bias correction using GAMs and nested vine copulas (NVC)
#' @description
#' Performs multivariate bias correction of climate model simulations using
#' Probability Integral Transforms (PIT) derived from Generalized additive models
#' (GAMs) and nested vine copulas (NVC). The method is designed
#' for potentially zero-inflated climate variables and preserves inter-variable,
#' spatial and temporal dependence structures.
#'
#' The algorithm consists of the following steps:
#'   1. Decompose modeled and observed data into seasonal (mean) and remainder
#'      components using GAMs and transform remainder components to PIT space
#'   2. Fit spatial vine copulas derived by NVC merging algorithm to PIT
#'      residuals of model and reference data
#'   3. Map model PITs to reference PITs via (inverse) Rosenblatt transforms
#'   4. Backtransform corrected PITs to data space using GAM parameters and
#'      climate adjustment
#'
#' @param mp Model simulations during projection period (data frame)
#' @param mc Model simulations during calibration period (data frame)
#' @param rc Observed reference data during calibration period (data frame)
#' @param var_names Variable names to be corrected (default: colnames(rc))
#' @param time_c Time vector for calibration period
#' @param time_p Time vector for projection period
#' @param locs Data frame with spatial coordinates (Lon, Lat) and Id per location,
#'             and potential other spatial covariates (e.g., altitude)
#' @param nrows Number of grid rows
#' @param ncols Number of grid columns
#' @param families Distribution families for GAM fitting as list: e.g., `list(tas = gaussian(), pr = tw(link = "log"))`
#' @param fixed Logical; whether to assume fixed local vine structure
#' @param mask Logical; whether to apply spatial adjacency mask (see Details)
#' @param bridge_var Optional variable used to bridge vine structures (see Details)
#' @param trunc_lvl Optional truncation level for vine copulas
#' @param seed Random seed
#' @param cores Number of cores for parallel computation
#' @param ... Additional arguments passed to `rvinecopulib::vinecop()`
#'
#' @returns  A list containing:
#'  - `corrected_mp` : Bias-corrected projection data
#'  - `rvine_mp`     : Nested vine copula fitted to model PITs
#'  - `rvine_rc`     : Nested vine copula fitted to reference PITs
#' @export
#' @details
#'  - `mask`: if `TRUE`, an adjacency mask is
#'   constructed based on the provided grid structure (`nrows`, `ncols`),
#'   and candidate edges in the first tree are limited to neighboring
#'   grid cells (including diagonal neighbors). This enforces spatial
#'   locality in the maximum spanning tree selection and reduces the
#'   number of admissible pair-copulas. If `FALSE`, the vine structure
#'   is selected without spatial constraints.
#'  - `bridge_var`: The bridging location, i.e., an integer between 1 and the
#'  number of locations, acts as a shared conditioning node that ensures the resulting
#'  global vine structure satisfies the proximity condition and
#'  remains a valid R-vine. If `NULL`, the bridging location is chosen randomly.
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
#' Lat = runif(n_locs, -90, 90)
#' )
#'
#' mp_corrected = gn_vbc(
#'       mp = mp,
#'       mc = mc,
#'       rc = rc,
#'       time_p = time_p,
#'       time_c = time_c,
#'       locs = locs,
#'       nrows = 5,
#'       ncols = 2,
#'       families = list("tas" = gaussian(), "pr" = Gamma(link = "log")),
#'       fixed = FALSE,
#'       mask = TRUE,
#'       bridge_var = 3,
#'       trunc_lvl = NULL,
#'       seed = 1,
#'       cores = 5,
#'       family_set = "tll" # for non-parametric vine copula estimation
#' )
#'
#' # plot(mp_corrected$rvine_mp, var_names = "use")
#'
gn_vbc <- function(
  mp,
  mc,
  rc,
  var_names = colnames(rc),
  time_c, # time vector of calibration period
  time_p, # time vector of projection period
  locs, # data frame with Lon, Lat, and Id
  nrows,
  ncols,
  families,
  fixed = TRUE, # indicator for fixed local vine structure
  mask = TRUE, # apply spatial adjacency mask
  bridge_var = NULL,
  trunc_lvl = NULL,
  seed = 123,
  cores = 11,
  ...
) {
  # Extract unique variable names (strip location suffixes)
  vars_unique = unique(sub("\\..*$", "", var_names))

  # ---------------------------------------------------------------------------
  # Step 1: GAM-based decomposition into seasonal and remainder components
  # ---------------------------------------------------------------------------

  gam_fit = get_GAMs(
    mp,
    mc,
    rc,
    locs,
    time_c,
    time_p,
    vars_unique,
    families,
    cores = cores
  )

  # Extract seasonal means and PIT-based remainders
  components_mc = extract_components(gam_fit$mc, mc, locs, time_c)
  components_mp = extract_components(gam_fit$mp, mp, locs, time_p)
  components_rc = extract_components(gam_fit$rc, rc, locs, time_c)

  # ---------------------------------------------------------------------------
  # Step 2: Fit spatial vine copula to model PIT residuals using NVC
  # ---------------------------------------------------------------------------

  mpu <- get_nested_vine(
    components_mp$remainder,
    nrows = nrows,
    ncols = ncols,
    fixed = fixed,
    mask = mask,
    bridge_var = bridge_var,
    ids = locs$Id,
    trunc_lvl = trunc_lvl,
    seed = seed,
    cores = cores,
    ...
  )

  rcu <- get_nested_vine(
    components_rc$remainder,
    nrows = nrows,
    ncols = ncols,
    fixed = fixed,
    mask = mask,
    bridge_var = bridge_var,
    ids = locs$Id,
    trunc_lvl = trunc_lvl,
    seed = seed,
    cores = cores,
    ...
  )

  # ---------------------------------------------------------------------------
  # Step 3: (Inverse) Rosenblatt transform
  # ---------------------------------------------------------------------------

  #  Rosenblatt transform model PITs using model vine
  u = rvinecopulib::rosenblatt(
    as.matrix(components_mp$remainder),
    mpu$vine_level3
  )

  # Inverse Rosenblatt transform to reference dependence structure
  u_mph = rvinecopulib::inverse_rosenblatt(
    u,
    rcu$vine_level3
  )

  # ---------------------------------------------------------------------------
  # Step 4: Backtransform corrected PITs to data space via inverse CDFs
  # ---------------------------------------------------------------------------

  # Prepare corrected PITs for backtransformation
  u_mph_wide = transform_to_wide_format(
    data.frame(u_mph),
    locs,
    vars_unique,
    time_p
  )

  x_mph_wide = inverse_PITs(gam_fit, u_mph_wide, u_mph, components_rc)

  # Reshape corrected data back to final wide format
  x_mph = x_mph_wide |>
    tidyr::pivot_wider(
      id_cols = "time",
      names_from = "Id",
      values_from = tidyselect::all_of(vars_unique),
      names_sep = "."
    ) |>
    dplyr::select(-"time")

  return(list(
    corrected_mp = x_mph,
    rvine_mp = mpu$vine_level3,
    rvine_rc = rcu$vine_level3
  ))
}
