#' @title Fit nested hierarchical vine copula structure for spatial climate data
#' @description Constructs a 3-level hierarchical R-vine copula model designed
#' for multivariate climate variables observed on a spatial grid.
#'
#' The hierarchy consists of:
#'  - Level 1: Spatial dependence within each climate variable
#'  - Level 2: Inter-variable dependence at a selected bridging location
#'  - Level 3: Merged vine capturing full dependence structure
#' @param data Data frame of pseudo-observations in `[0, 1]`, with columns
#' named as "variable.location"
#' @param nrows Number of rows in the spatial grid
#' @param ncols Number of columns in the spatial grid
#' @param var_types Character vector specifying marginal variable types
#' (e.g., "c" for continuous)
#' @param fixed Logical; if TRUE, fit one spatial vine structure and reuse it
#' for all variables
#' @param mask Logical; if TRUE, restrict spatial edges to neighboring locations
#' @param bridge_var Integer index of the spatial location used to connect
#' variables across levels; randomly chosen if NULL
#' @param ids Optional vector of location identifiers; defaults to full grid
#' @param trunc_lvl Integer truncation level for vine copulas
#' @param seed Integer random seed for reproducibility
#' @param ... Additional arguments passed to `rvinecopulib::vinecop()`
#'
#' @returns A list containing:
#'  - bridge_var  : Selected bridge location
#'  - rvs_level1  : Spatial R-vine structure(s)
#'  - rvs_level2  : Inter-variable R-vine structure
#'  - rvs_level3  : Final merged R-vine structure
#'  - vine_level1 : Fitted spatial vine(s), if applicable
#'  - vine_level2 : Fitted Level 2 vine, if requested
#'  - vine_level3 : Final fitted vine copula, if requested
#' @export
#'
#' @examples
#' data = matrix(runif(600), ncol = 6)
#' colnames(data) = c("tas.1", "tas.2", "tas.3", "pr.1", "pr.2", "pr.3")
#' test = get_nested_vine(data, nrows = 2, ncols = 2, fixed = FALSE,
#'                         ids = c(1,2,3), family_set = "tll", cores = 12)
#'
get_nested_vine <- function(
    data,
    nrows, ncols,
    var_types = rep("c", ncol(data)),
    fixed = TRUE,
    mask = TRUE,
    bridge_var = NULL,
    ids = NULL,
    trunc_lvl = NULL,
    seed = 123,
    ...
) {

  # Ensure reproducibility
  set.seed(seed)
  # ---------------------------------------------------------------------------
  # Step 1: Initialize dimensions and variable structure
  # ---------------------------------------------------------------------------

  # Number of spatial locations
  nlocs <- if (!is.null(ids)) length(ids) else nrows * ncols

  # Set truncation level if not provided
  trunc_lvl <- if (is.null(trunc_lvl)) nlocs else trunc_lvl

  # Extract unique climate variable names
  vars <- unique(sub("\\..*$", "", colnames(data)))
  nvars <- length(vars)

  # Total dimensionality of the joint copula
  d <- nlocs * nvars

  # ---------------------------------------------------------------------------
  # Step 2: Select or validate bridge location
  # ---------------------------------------------------------------------------

  # Bridge location connects spatial and inter-variable dependence structures
  if (is.null(bridge_var)) {
    bridge_var <- sample(1:nlocs, 1)
  } else if (!(bridge_var %in% seq_len(nlocs))) {
    stop(paste("Bridging variable must be between 1 and", nlocs))
  }

  # Initialize output container
  output <- list(bridge_var = bridge_var)

  # Convert data to matrix for efficient access
  data_matrix <- as.matrix(data)

  # ---------------------------------------------------------------------------
  # Step 3: Pre-compute spatial adjacency mask (optional)
  # ---------------------------------------------------------------------------

  # Spatial mask restricts edges to neighboring grid cells
  global_mask <- if (mask) get_spatial_mask(nrows, ncols, ids) else NULL

  # ---------------------------------------------------------------------------
  # Helper: Map variable index to column indices in data matrix
  # ---------------------------------------------------------------------------
  get_var_indices <- function(iVar) {
    seq(iVar + (iVar - 1) * (nlocs - 1), iVar + iVar * (nlocs - 1))
  }

  # ---------------------------------------------------------------------------
  # Helper: Fit spatial vine copula for a single variable
  # ---------------------------------------------------------------------------
  fit_var_vine <- function(iVar, var, ...) {

    # Column indices corresponding to the variable across locations
    idx <- get_var_indices(iVar)

    # Extract data for the current variable
    var_data <- as.data.frame(data_matrix[, grep(var, colnames(data))])

    # Initialize structure (NULL = automatic selection)
    structure <- NA

    # -------------------------------------------------------------------------
    # Optional spatial structure construction via minimum spanning tree
    # -------------------------------------------------------------------------
    if (mask) {

      # Compute Kendall's tau matrix as dependence measure
      ktau_matrix <- wdm::wdm(data_matrix[, idx], method = "kendall")

      # Convert dependence to edge weights (low weight = strong dependence)
      weight_matrix <- global_mask * (1 - abs(ktau_matrix))

      # Build weighted graph and extract minimum spanning tree
      g <- igraph::graph_from_adjacency_matrix(
        as.matrix(weight_matrix),
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
      )
      tree <- igraph::mst(g)

      # Convert spanning tree to R-vine structure
      structure <- spanning_tree_to_rvine_structure(tree)
    }

    # Fit vine copula with optional fixed structure
    rvinecopulib::vinecop(
      var_data,
      trunc_lvl = trunc_lvl,
      structure = structure,
      var_types = var_types[idx],
      ...
    )
  }

  # ---------------------------------------------------------------------------
  # Level 1: Spatial dependence within variables
  # ---------------------------------------------------------------------------
  rvs_level1 <- if (fixed) {

    # Fixed mode: fit spatial structure for one random variable
    iVar <- sample(1:nvars, 1)
    vine <- fit_var_vine(iVar, vars[iVar], ...)

    output[["vine_level1"]] = vine
    output[["rvs_level1"]] = vine$structure

    vine$structure

  } else {

    # Flexible mode: fit separate spatial structure for each variable
    result_vines <- lapply(seq_len(nvars), function(iVar) {
      vine <- fit_var_vine(iVar, vars[iVar], ...)
      output[[paste0("vine_level1_", vars[iVar])]] <<- vine
      output[[paste0("rvs_level1_", vars[iVar])]] <<- vine$structure
      vine$structure
    })

    result_vines
  }

  # ---------------------------------------------------------------------------
  # Level 2: Inter-variable dependence at bridge location
  # ---------------------------------------------------------------------------
  vine_level2 <- rvinecopulib::vinecop(
    data[, seq(bridge_var, ncol(data), nlocs)],
    trunc_lvl = trunc_lvl,
    var_types = var_types[seq(bridge_var, d, nlocs)],
    ...
  )

  rvs_level2 <- vine_level2$structure
  output[["rvs_level2"]] <- rvs_level2
  output[["vine_level2"]] <- vine_level2

  # ---------------------------------------------------------------------------
  # Level 3: Merge spatial and inter-variable vine structures
  # ---------------------------------------------------------------------------
  rvs_level3 <- if (fixed) {
    merge_edges_fixed_full(rvs_level1, rvs_level2, bridge_var)
  } else {
    merge_edges_individual_full(rvs_level1, rvs_level2, bridge_var)
  }

  output[["rvs_level3"]] <- rvs_level3

  # ---------------------------------------------------------------------------
  # Final step: Fit copula parameters for merged structure
  # ---------------------------------------------------------------------------
  vine_level3 <- rvinecopulib::vinecop(
    data,
    trunc_lvl = trunc_lvl,
    structure = rvs_level3,
    var_types = var_types,
    ...
  )
  output[["vine_level3"]] <- vine_level3
  output[["rvs_level3"]] <- vine_level3$structure


  return(output)
}


