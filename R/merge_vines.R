.validate_merge_structure <- function(x, name) {
  if (!inherits(x, "rvine_structure")) {
    stop(name, " must be an rvine_structure.", call. = FALSE)
  }
  if (length(x$d) != 1L || !is.finite(x$d) || x$d < 2L ||
      length(x$trunc_lvl) != 1L || !is.finite(x$trunc_lvl) ||
      x$trunc_lvl < 1L) {
    stop(name, " must be a non-degenerate, non-empty vine structure.",
      call. = FALSE
    )
  }
  invisible(x)
}

.validate_bridge <- function(bridge_var, spatial_dimension) {
  if (length(bridge_var) != 1L || is.na(bridge_var) ||
      bridge_var != as.integer(bridge_var) || bridge_var < 1L ||
      bridge_var > spatial_dimension) {
    stop("bridge_var must be an integer between 1 and ", spatial_dimension, ".",
      call. = FALSE
    )
  }
  as.integer(bridge_var)
}

.merge_with_maps <- function(components, maps, global_dimension) {
  merge_rvine_structures(components, maps, as.integer(global_dimension))
}

#' @title Merge vine structures for hierarchical copula
#'   (fixed spatial structure)
#' @description Combines a spatial R-vine structure replicated for each
#'   variable with an inter-variable R-vine structure at one bridge location.
#' @param rvs_level1 Spatial R-vine structure to replicate.
#' @param rvs_level2 Inter-variable R-vine structure.
#' @param bridge_var Integer location index used to connect the hierarchies.
#' @returns A merged `rvine_structure`.
#' @export
merge_edges_fixed_full <- function(rvs_level1, rvs_level2, bridge_var) {
  .validate_merge_structure(rvs_level1, "rvs_level1")
  .validate_merge_structure(rvs_level2, "rvs_level2")

  spatial_dimension <- as.integer(rvs_level1$d)
  variable_dimension <- as.integer(rvs_level2$d)
  bridge_var <- .validate_bridge(bridge_var, spatial_dimension)

  # Each component remains valid in its local labels. C++ obtains its
  # original-label edge list and applies these maps directly.
  spatial_maps <- lapply(seq_len(variable_dimension), function(variable) {
    (variable - 1L) * spatial_dimension + seq_len(spatial_dimension)
  })
  bridge_map <- bridge_var +
    (seq_len(variable_dimension) - 1L) * spatial_dimension

  .merge_with_maps(
    c(rep(list(rvs_level1), variable_dimension), list(rvs_level2)),
    c(spatial_maps, list(bridge_map)),
    spatial_dimension * variable_dimension
  )
}

#' @title Merge vine structures for hierarchical copula
#'   (individual spatial structures)
#' @description Combines one spatial R-vine structure per variable with an
#'   inter-variable R-vine structure at one bridge location.
#' @param rvs_level1 List of spatial R-vine structures, one per variable.
#' @param rvs_level2 Inter-variable R-vine structure.
#' @param bridge_var Integer location index used to connect the hierarchies.
#' @returns A merged `rvine_structure`.
#' @export
merge_edges_individual_full <- function(rvs_level1, rvs_level2, bridge_var) {
  .validate_merge_structure(rvs_level2, "rvs_level2")
  variable_dimension <- as.integer(rvs_level2$d)
  if (!is.list(rvs_level1) || length(rvs_level1) != variable_dimension) {
    stop("rvs_level1 must contain one spatial structure per level-2 variable.",
      call. = FALSE
    )
  }
  for (index in seq_along(rvs_level1)) {
    .validate_merge_structure(rvs_level1[[index]],
      paste0("rvs_level1[[", index, "]]"))
  }

  spatial_dimension <- as.integer(rvs_level1[[1]]$d)
  if (!all(vapply(rvs_level1, function(x) x$d == spatial_dimension, logical(1)))) {
    stop("All individual spatial structures must have the same dimension.",
      call. = FALSE
    )
  }
  bridge_var <- .validate_bridge(bridge_var, spatial_dimension)

  spatial_maps <- lapply(seq_len(variable_dimension), function(variable) {
    (variable - 1L) * spatial_dimension + seq_len(spatial_dimension)
  })
  bridge_map <- bridge_var +
    (seq_len(variable_dimension) - 1L) * spatial_dimension

  .merge_with_maps(
    c(rvs_level1, list(rvs_level2)),
    c(spatial_maps, list(bridge_map)),
    spatial_dimension * variable_dimension
  )
}
