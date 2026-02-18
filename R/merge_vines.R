#' @title Merge vine structures for hierarchical copula
#'   (fixed spatial structure)
#' @description Combines level 1 (spatial dependence within variables)
#'   and level 2 (inter-variable dependence) R-vine structures into a
#'   single unified structure using a fixed spatial structure approach.
#'   The level 2 structure is replicated across variables, then merged
#'   with the level 1 structure using a bridging variable.
#' @param rvs_level1 R-vine structure object for spatial dependence
#'   (to be replicated)
#' @param rvs_level2 R-vine structure object for inter-variable
#'   dependence (single structure)
#' @param bridge_var Integer. Location index used to connect the two
#'   hierarchies.
#'
#' @returns Merged R-vine structure (rvs_level3) combining spatial and
#'   inter-variable dependence with variables properly ordered.
#' @details This function replicates the level 2 structure across all
#'   spatial locations, then transforms both structures to a common
#'   natural ordering before merging. The bridge variable anchors the
#'   connection between the two levels.
#' @export
merge_edges_fixed_full <- function(rvs_level1, rvs_level2, bridge_var) {
  # Step 1: Replicate level 2 structure for each spatial location
  list_rvs_level1 <- list(rvs_level1)

  for (i in 2:rvs_level2$d) {
    rvs_level1_tmp <- rvs_level1
    # Offset structure indices to correspond to location i
    rvs_level1_tmp$order <- rvs_level1_tmp$order + (i - 1) * rvs_level1$d
    rvs_level1_tmp$struct_array <- lapply(
      rvs_level1_tmp$struct_array,
      function(x) x + (i - 1) * rvs_level1$d
    )
    list_rvs_level1[[i]] <- rvs_level1_tmp
  }

  lev1_tmp <- list_rvs_level1

  # Step 2: Extract combined ordering from all level 2 structures
  order_level1 <- unlist(lapply(list_rvs_level1, function(x) x$order))

  # Step 3: Relabel all level 2 structures to natural order (1, 2, 3, ...)
  list_rvs_level1 <- lapply(lev1_tmp, function(x) {
    x$order <- x$order[order(x$order)] # Convert to natural ordering
    x
  })

  # Step 4: Prepare level 1 structure for merging
  # Identify bridge variable indices for each spatial location
  bridges <- bridge_var + rvs_level1_tmp$d * seq(0, rvs_level2$d - 1)

  # Reindex level 1 structure relative to bridge variables
  order_tmp <- bridges[rvs_level2$order]

  rvs_level2$order <- order_tmp
  rvs_level2$struct_array <- lapply(
    rvs_level2$struct_array,
    function(x) order_tmp[x]
  )
  rvs_level2_orig <- rvs_level2

  # Step 5: Transform level 1 to natural ordering
  inverse_map <- integer(length(order_level1))
  inverse_map[order_level1] <- seq_along(order_level1)

  rvs_level2$order <- inverse_map[rvs_level2$order]
  rvs_level2$struct_array <- lapply(
    rvs_level2$struct_array,
    function(x) inverse_map[x]
  )

  # Step 6: Merge level 1 and level 2 structures
  rvs_level3_tmp <- merge_rvine_structures(c(
    list_rvs_level1,
    list(rvs_level2)
  ))

  # Step 7: Relabel merged structure back to original variable-location order
  rvs_level3 <- rvs_level3_tmp
  rvs_level3$order <- order_level1[rvs_level3_tmp$order]

  rvs_level3
}

#' @title Merge vine structures for hierarchical copula
#'   (individual spatial structures)
#' @description Combines level 1 (spatial dependence within variables)
#'   and level 2 (inter-variable dependence) R-vine structures into a
#'   single unified structure using an individual spatial structure
#'   approach. Unlike the fixed version, each variable can have its own
#'   level 1 structure. Level 2 structures are replicated and merged with
#'   their corresponding level 1 structures.
#' @param rvs_level1 List of R-vine structure objects for spatial dependence
#' (one structure per variable)
#' @param rvs_level2 R-vine structure object for inter-variable dependence
#' @param bridge_var Integer. Location index used to connect the hierarchies.
#'
#' @returns Merged R-vine structure (rvs_level3) combining spatial and
#'   inter-variable dependence with variables properly ordered.
#' @details This function handles variable-specific level 1 structures by
#'   replicating each corresponding level 2 structure, then transforming
#'   both to natural ordering before merging. The bridge variable anchors
#'   the connection between the two levels.
#' @export
#'
merge_edges_individual_full <- function(rvs_level1, rvs_level2, bridge_var) {
  # Step 1: Replicate level 2 structures, one for each variable
  list_rvs_level1 <- list(rvs_level1[[1]])

  for (i in 2:rvs_level2$d) {
    rvs_level1_tmp <- rvs_level1[[i]]
    # Offset structure indices to correspond to variable i at location i
    rvs_level1_tmp$order <- rvs_level1_tmp$order + (i - 1) * rvs_level1_tmp$d
    rvs_level1_tmp$struct_array <- lapply(
      rvs_level1_tmp$struct_array,
      function(x) x + (i - 1) * rvs_level1_tmp$d
    )
    list_rvs_level1[[i]] <- rvs_level1_tmp
  }

  lev1_tmp <- list_rvs_level1

  # Step 2: Extract combined ordering from all level 2 structures
  order_level1 <- unlist(lapply(list_rvs_level1, function(x) x$order))

  # Step 3: Relabel all level 2 structures to natural order (1, 2, 3, ...)
  list_rvs_level1 <- lapply(lev1_tmp, function(x) {
    x$order <- x$order[order(x$order)] # Convert to natural ordering
    x
  })

  # Step 4: Prepare level 1 structure for merging
  # Identify bridge variable indices for each spatial location
  bridges <- bridge_var + rvs_level1_tmp$d * seq(0, rvs_level2$d - 1)

  # Reindex level 1 structure relative to bridge variables
  order_tmp <- bridges[rvs_level2$order]

  rvs_level2$order <- order_tmp
  rvs_level2$struct_array <- lapply(
    rvs_level2$struct_array,
    function(x) order_tmp[x]
  )
  rvs_level2_orig <- rvs_level2

  # Step 5: Transform level 1 to natural ordering
  inverse_map <- integer(length(order_level1))
  inverse_map[order_level1] <- seq_along(order_level1)

  rvs_level2$order <- inverse_map[rvs_level2$order]
  rvs_level2$struct_array <- lapply(
    rvs_level2$struct_array,
    function(x) inverse_map[x]
  )

  # Step 6: Merge level 1 and level 2 structures
  rvs_level3_tmp <- merge_rvine_structures(c(
    list_rvs_level1,
    list(rvs_level2)
  ))

  # Step 7: Relabel merged structure back to original variable-location order
  rvs_level3 <- rvs_level3_tmp
  rvs_level3$order <- order_level1[rvs_level3_tmp$order]

  rvs_level3
}
