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
#' @param families Named list of family objects for each variable (e.g., list(tas = gaussian(), pr = Tweedie()))
#' @param cores Number of parallel workers used for model fitting; if NULL or 1, runs sequentially
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
#' locs = data.frame(Lon = runif(5, -180, 180), Lat = runif(5, -90, 90), Id = 1:5)
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
get_GAMs <- function(
    mp, mc, rc,
    locs,
    time_c, time_p,
    var_names,
    families,
    cores = NULL
) {

  # ---------------------------------------------------------------------------
  # Step 1: Transform all datasets into GAM-ready wide format
  # ---------------------------------------------------------------------------

  # Each dataset is expanded to include spatial coordinates and
  # day-of-year for cyclic temporal smoothing

  dfs <- list(
    mc = transform_to_wide_format(mc, locs, var_names, time_c),
    mp = transform_to_wide_format(mp, locs, var_names, time_p),
    rc = transform_to_wide_format(rc, locs, var_names, time_c)
  )

  # ---------------------------------------------------------------------------
  # Step 2: Define parallel execution plan and model combinations
  # ---------------------------------------------------------------------------

  inputs <- expand.grid(
    dataset = names(dfs),
    var = var_names,
    stringsAsFactors = FALSE
  )

  # Save and restore future plan
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (!is.null(cores) && cores > 1) {
    future::plan(future::multisession, workers = cores)
  }

  # ---------------------------------------------------------------------------
  # Step 3: Fit GAMs in parallel
  # ---------------------------------------------------------------------------

  # Model structure:
  #   response ~ te(t, Lat, Lon)
  # with:
  #   - cyclic cubic spline for day-of-year (t)
  #   - thin plate splines for spatial dimensions

  if (!is.null(cores) && cores > 1) {

    fits <- furrr::future_pmap(
      inputs,
      function(dataset, var) {
        mgcv::bam(
          stats::as.formula(
            paste0(var, " ~ te(t, Lat, Lon, bs = c('cc', 'tp', 'tp'))")
          ),
          data = dfs[[dataset]],
          family = families[[var]]
        )
      },
      .options = furrr::furrr_options(seed = TRUE)
    )

  } else {

    fits <- purrr::pmap(
      inputs,
      function(dataset, var) {
        mgcv::bam(
          stats::as.formula(
            paste0(var, " ~ te(t, Lat, Lon, bs = c('cc', 'tp', 'tp'))")
          ),
          data = dfs[[dataset]],
          family = families[[var]]
        )
      }
    )
  }


  # ---------------------------------------------------------------------------
  # Step 4: Organize fitted models into nested list
  # ---------------------------------------------------------------------------

  split(fits, inputs$dataset) |>
    purrr::map2(split(inputs$var, inputs$dataset), purrr::set_names)
}
