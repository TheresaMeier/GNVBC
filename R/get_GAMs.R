#' @title Fit tensor-product GAMs for multiple climate datasets
#' @description Fits generalized additive models (GAMs) using tensor-product smooths to
#' capture joint temporal and spatial variation in climate variables.
#'
#' Temporal variation is modeled using a cyclic cubic spline on day-of-year,
#' while spatial variation is modeled using thin plate splines over latitude
#' and longitude.
#'
#' Models are fitted separately for multiple datasets (e.g., modeled
#' calibration mc, modeled projection mp, reference calibration rc) and for multiple
#' climate variables, using parallel processing if requested.
#'
#' @param mp Data frame with modeled projection data with columns named "variable.location_id"
#' @param mc Data frame with modeled calibration data with columns named "variable.location_id"
#' @param rc Data frame with reference calibration data with columns named "variable.location_id"
#' @param locs Data frame of spatial locations with columns: Id, Lat, Lon
#' @param time_c Date vector corresponding to rows of calibration datasets (mc, rc)
#' @param time_p Date vector corresponding to rows of projection dataset (mp)
#' @param var_names Character vector of climate variable names (e.g., c("tas", "pr"))
#' @param families Named list of family objects for each variable (e.g., list(tas = gaussian(), pr = tw(link = "log")))
#' @param cores Number of parallel workers used for model fitting; if NULL or 1, runs sequentially
#' @param extra_smooths Optional character vector of additional predictor variable names to include as smooth terms (e.g., c("elev") or c("s(elev, k=5)"));
#' if provided, these will be added as separate smooth terms in the GAM formula.
#'
#' @importFrom mgcv bam
#'
#' @returns Nested list of fitted GAM objects with structure:`[[dataset]][[variable]]`
#' @export
#'
#' @examples
#' set.seed(1)
#' # Simulate some data for demonstration
#' mp = data.frame(matrix(rnorm(500), ncol=10))
#' mc = data.frame(matrix(rnorm(1000), ncol=10))
#' rc = data.frame(matrix(rnorm(1000), ncol=10))
#' colnames(mp) = colnames(mc) = colnames(rc) = paste0(rep(c("tas.", "pr."), each = 5), c(1:5))
#'
#' time_c = as.Date("2000-01-01") + 0:99
#' time_p = as.Date("2020-01-01") + 0:49
#'
#' # Simulate 5 locations
#' locs = data.frame(Lon = runif(5, -180, 180), Lat = runif(5, -90, 90), Id = 1:5,
#'                   Altitude = runif(5, 0, 3000))
#'
#' var_names = c("tas", "pr")
#'
#' families = list("tas" = gaussian(), "pr" = gaussian())
#'
#' test = get_GAMs(
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
#' # With extra smooths
#' test = get_GAMs(
#' mp = mp,
#' mc = mc,
#' rc = rc,
#' time_p = time_p,
#' time_c = time_c,
#' locs = locs,
#' var_names = var_names,
#' families = families,
#' cores = 5,
#' extra_smooths = c("s(Altitude, k=3)")
#' )
#'
get_GAMs <- function(
  mp,
  mc,
  rc,
  locs,
  time_c,
  time_p,
  var_names,
  families,
  cores = NULL,
  extra_smooths = NULL
) {
  # ---------------------------------------------------------------------------
  # Step 1: Transform datasets
  # ---------------------------------------------------------------------------

  dfs <- list(
    mc = transform_to_wide_format(mc, locs, var_names, time_c),
    mp = transform_to_wide_format(mp, locs, var_names, time_p),
    rc = transform_to_wide_format(rc, locs, var_names, time_c)
  )

  # ---------------------------------------------------------------------------
  # Step 2: Prepare model combinations
  # ---------------------------------------------------------------------------

  inputs <- expand.grid(
    dataset = names(dfs),
    var = var_names,
    stringsAsFactors = FALSE
  )

  # Save & restore future plan
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (!is.null(cores) && cores > 1) {
    future::plan(future::multisession, workers = cores)
  }

  # ---------------------------------------------------------------------------
  # Step 3: Fit GAMs
  # ---------------------------------------------------------------------------

  fit_one <- function(dataset, var) {
    form <- build_gam_formula(var, extra_smooths)

    mgcv::bam(
      formula = form,
      data = dfs[[dataset]],
      family = families[[var]]
    )
  }

  if (!is.null(cores) && cores > 1) {
    fits <- furrr::future_pmap(
      inputs,
      fit_one,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    fits <- purrr::pmap(
      inputs,
      fit_one
    )
  }

  # ---------------------------------------------------------------------------
  # Step 4: Return nested list
  # ---------------------------------------------------------------------------

  split(fits, inputs$dataset) |>
    purrr::map2(
      split(inputs$var, inputs$dataset),
      purrr::set_names
    )
}

#' Build GAM formula with tensor-product smooths for time and space, plus optional extra smooths.
#'
#' @param response Name of the response variable (e.g., "tas", "pr")
#' @param extra_smooths Optional character vector of additional predictor variable names to include as smooth terms (e.g., c("elev", "dist_to_coast"))
#'
#' @returns A formula object for use in mgcv::bam(), with a base tensor-product smooth for time and space, plus any specified extra smooths.
#' @export
build_gam_formula <- function(response, extra_smooths = NULL) {
  # Base smooth: time × space
  base_term <- "te(t, Lat, Lon, bs = c('cc', 'tp', 'tp'))"

  rhs_terms <- base_term

  if (!is.null(extra_smooths)) {
    extra_terms <- purrr::map_chr(
      extra_smooths,
      function(x) {
        if (grepl("\\(", x)) {
          x # full mgcv term supplied by user
        } else {
          paste0("s(", x, ")") # shorthand
        }
      }
    )

    rhs_terms <- paste(c(rhs_terms, extra_terms), collapse = " + ")
  }

  stats::as.formula(paste(response, "~", rhs_terms))
}
