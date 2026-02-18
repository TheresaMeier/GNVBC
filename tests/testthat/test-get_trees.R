test_that("spatial grid derivation works", {
  mask <- get_spatial_mask(4, 6)

  # Correct class and dimensions
  expect_s4_class(mask, "dgCMatrix")
  expect_equal(dim(mask), c(4 * 6, 4 * 6))

  # Symmetry and no self-loops
  expect_equal(mask, Matrix::t(mask))
  expect_equal(Matrix::diag(mask), rep(0, nrow(mask)))

  # Correct neighbor counts
  # Helper to compute 1D index
  idx <- function(r, c) (r - 1) * 4 + c

  # Corner cell (1,1)
  expect_equal(sum(mask[idx(1, 1), ]), 3)
  # Edge cell (1,3)
  expect_equal(sum(mask[idx(1, 3), ]), 5)
  # Interior cell (3,3)
  expect_equal(sum(mask[idx(3, 3), ]), 8)

  # Correct for 2x2 grid
  mask <- get_spatial_mask(2, 2)

  expected <- Matrix::sparseMatrix(
    i = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4),
    j = c(2, 3, 4, 1, 3, 4, 1, 2, 4, 1, 2, 3),
    x = rep(1, 12),
    dims = c(4, 4)
  )

  expect_equal(mask, expected)

  # Subsetting by ids preserves adjacency structure
  ids <- c(1, 5, 6, 10, 15)
  full_mask <- get_spatial_mask(5, 3)
  sub_mask <- get_spatial_mask(5, 3, ids)

  expect_equal(dim(sub_mask), c(length(ids), length(ids)))
  expect_equal(sub_mask, full_mask[ids, ids])
})

test_that("rvine_from_spanning_tree returns valid R-vine structure", {
  set.seed(1)
  g <- igraph::make_tree(5, children = 2, mode = "undirected")

  rvs <- rvine_from_spanning_tree(g)

  expect_s3_class(rvs, "rvine_structure")
  expect_equal(rvs$d, igraph::vcount(g))

  # Function errors for disconnected graph
  g <- igraph::make_graph(c(1, 2, 2, 3, 3, 4, 5, 6), directed = FALSE)

  expect_error(
    rvine_from_spanning_tree(g),
    "Input spanning_tree must be a connected graph."
  )

  # R-vine order is a permutation of graph nodes
  g <- igraph::make_tree(6, children = 2, mode = "undirected")

  rvs <- rvine_from_spanning_tree(g)

  expect_equal(
    sort(rvs$order),
    seq_len(igraph::vcount(g))
  )

  # First-tree structure has correct size and nodes
  rvs <- rvine_from_spanning_tree(g)

  first_tree <- rvs$struct_array[[1]]

  expect_equal(length(first_tree), igraph::vcount(g) - 1)
  expect_true(all(first_tree %in% seq_len(igraph::vcount(g))))

  # Known small tree produces expected R-vine structure
  g <- igraph::make_graph(
    edges = c(1, 2, 2, 3, 2, 4),
    directed = FALSE
  )

  rvs <- rvine_from_spanning_tree(g)

  expect_equal(sort(rvs$order), 1:4)
  expect_equal(length(rvs$struct_array[[1]]), 3)
})
