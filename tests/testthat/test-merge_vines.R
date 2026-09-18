edge_sets <- function(rvs, local_to_global = seq_len(rvs$d)) {
  stopifnot(length(local_to_global) == rvs$d)
  lapply(seq_len(rvs$trunc_lvl), function(tree) {
    sort(vapply(seq_len(rvs$d - tree), function(edge) {
      conditioning_indices <- if (tree == 1L) {
        integer()
      } else {
        vapply(rvs$struct_array[seq_len(tree - 1L)], `[[`, numeric(1), edge)
      }
      conditioned <- local_to_global[rvs$order[c(
        edge, rvs$struct_array[[tree]][edge]
      )]]
      conditioning <- local_to_global[rvs$order[conditioning_indices]]
      paste(
        paste(sort(conditioned), collapse = ","),
        paste(sort(conditioning), collapse = ","),
        sep = "|"
      )
    }, character(1)))
  })
}

canonical_edge_sets <- edge_sets

expect_inherited_edges <- function(merged, components, maps) {
  merged_edges <- canonical_edge_sets(merged)
  for (component in seq_along(components)) {
    inherited_edges <- edge_sets(components[[component]], maps[[component]])
    for (tree in seq_along(inherited_edges)) {
      expect_true(
        all(inherited_edges[[tree]] %in% merged_edges[[tree]]),
        info = paste("component", component, "tree", tree)
      )
    }
  }
}

test_that("canonical edge sets decode equivalent matrix orientations", {
  # These are stored natural-order rows, as returned by the upstream wrapper.
  # Construct the lists directly because rvine_structure() accepts original-
  # label rows and converts them before storing them in this form.
  first_orientation <- structure(
    list(order = c(1, 3, 2), struct_array = list(c(3, 3), 2), d = 3, trunc_lvl = 2),
    class = c("rvine_structure", "list")
  )
  second_orientation <- structure(
    list(order = c(3, 1, 2), struct_array = list(c(3, 3), 2), d = 3, trunc_lvl = 2),
    class = c("rvine_structure", "list")
  )

  expected <- list(c("1,2|", "2,3|"), "1,3|2")
  expect_equal(canonical_edge_sets(first_orientation), expected)
  expect_equal(canonical_edge_sets(second_orientation), expected)
})

test_that("Merging of vine copulas works - fixed structure", {
  # Generate vine copulas for testing
  rvs_level1 <- rvinecopulib::rvine_structure(
    order = c(4, 2, 5, 3, 1),
    struct_array = list(
      c(5, 3, 1, 1),
      c(1, 1, 3),
      c(3, 5),
      2
    )
  )

  rvs_level2 <- rvinecopulib::rvine_structure(
    order = c(1, 3, 2),
    struct_array = list(
      c(2, 2),
      3
    )
  )

  bridge_var <- 2
  rvs_level3 <- merge_edges_fixed_full(rvs_level1, rvs_level2, bridge_var)

  rvs_level3_true <- rvinecopulib::rvine_structure(
    order = c(4, 9, 14, 15, 11, 13, 10, 6, 8, 12, 5, 1, 3, 2, 7),
    struct_array = list(
      c(5, 10, 15, 11, 13, 12, 6, 8, 7, 7, 1, 3, 2, 7),
      c(1, 6, 11, 13, 12, 7, 8, 7, 12, 2, 3, 2, 7),
      c(3, 8, 13, 12, 7, 8, 7, 12, 2, 3, 2, 7),
      c(2, 7, 12, 7, 8, 6, 12, 2, 3, 1, 7)
    )
  )

  expect_s3_class(rvs_level3, "rvine_structure")
  expect_equal(rvs_level3$d, rvs_level1$d * rvs_level2$d)
  expect_equal(sort(rvs_level3$order), seq_len(rvs_level3$d))
  expect_true(all(rvs_level3$order >= 1))
  expect_true(all(rvs_level3$order <= rvs_level3$d))

  # The upstream converter is free to use another matrix orientation. Compare
  # trees, not matrix cells, so this remains a structural regression fixture.
  expect_equal(canonical_edge_sets(rvs_level3), canonical_edge_sets(rvs_level3_true))
  expect_inherited_edges(
    rvs_level3,
    c(rep(list(rvs_level1), rvs_level2$d), list(rvs_level2)),
    c(
      lapply(seq_len(rvs_level2$d), function(i) {
        (i - 1L) * rvs_level1$d + seq_len(rvs_level1$d)
      }),
      list(bridge_var + (seq_len(rvs_level2$d) - 1L) * rvs_level1$d)
    )
  )
})

test_that("Merge works for different bridge variables - fixed structure", {
  # Generate vine copulas for testing
  rvs_level1 <- rvinecopulib::rvine_structure(
    order = c(2, 1, 3, 4),
    struct_array = list(
      c(3, 4, 4),
      c(4, 3),
      1
    )
  )

  rvs_level2 <- rvinecopulib::rvine_structure(
    order = c(2, 4, 1, 3, 5, 6),
    struct_array = list(
      c(3, 5, 6, 6, 6),
      c(6, 6, 5, 5),
      c(5, 3, 3),
      c(1, 1),
      4
    )
  )

  for (bridge_var in seq_len(rvs_level1$d)) {
    rvs_level3 <- merge_edges_fixed_full(rvs_level1, rvs_level2, bridge_var)

    expect_s3_class(rvs_level3, "rvine_structure")
    expect_equal(rvs_level3$d, rvs_level1$d * rvs_level2$d)
  }
})

test_that("Merging of vine copulas works - flexible structure", {
  # Generate vine copulas for testing
  rvs_level1 <- list(
    rvinecopulib::rvine_structure(
      order = c(4, 2, 5, 3, 1),
      struct_array = list(
        c(5, 3, 1, 1),
        c(1, 1, 3),
        c(3, 5),
        2
      )
    ),
    rvinecopulib::rvine_structure(
      order = c(2, 3, 4, 1, 5),
      struct_array = list(
        c(4, 1, 5, 5),
        c(5, 5, 1),
        c(1, 4),
        3
      )
    ),
    rvinecopulib::rvine_structure(
      order = c(1, 3, 2, 5, 4),
      struct_array = list(
        c(2, 5, 4, 4),
        c(4, 4, 5),
        c(5, 2),
        3
      )
    )
  )

  rvs_level2 <- rvinecopulib::rvine_structure(
    order = c(2, 1, 3),
    struct_array = list(
      c(3, 3),
      1
    )
  )

  bridge_var <- 3
  rvs_level3 <- merge_edges_individual_full(rvs_level1, rvs_level2, bridge_var)

  rvs_level3_true <- rvinecopulib::rvine_structure(
    order = c(4, 7, 11, 12, 14, 15, 9, 10, 6, 8, 2, 5, 1, 3, 13),
    struct_array = list(
      c(5, 9, 12, 14, 15, 13, 10, 6, 8, 13, 3, 1, 3, 13),
      c(1, 10, 14, 15, 13, 3, 6, 8, 13, 3, 1, 3, 13),
      c(3, 6, 15, 13, 3, 1, 8, 13, 3, 1, 5, 13),
      c(2, 8, 13, 3, 1, 5, 13, 3, 1, 5, 13)
    )
  )

  expect_s3_class(rvs_level3, "rvine_structure")
  expect_equal(rvs_level3$d, rvs_level1[[1]]$d * rvs_level2$d)
  expect_equal(sort(rvs_level3$order), seq_len(rvs_level3$d))
  expect_true(all(rvs_level3$order >= 1))
  expect_true(all(rvs_level3$order <= rvs_level3$d))

  expect_equal(canonical_edge_sets(rvs_level3), canonical_edge_sets(rvs_level3_true))
  expect_inherited_edges(
    rvs_level3,
    c(rvs_level1, list(rvs_level2)),
    c(
      lapply(seq_len(rvs_level2$d), function(i) {
        (i - 1L) * rvs_level1[[1]]$d + seq_len(rvs_level1[[1]]$d)
      }),
      list(bridge_var + (seq_len(rvs_level2$d) - 1L) * rvs_level1[[1]]$d)
    )
  )
})

test_that("different truncation levels produce a valid complete structure", {
  spatial <- rvinecopulib::rvine_structure(
    order = c(3, 1, 4, 2),
    struct_array = list(c(2, 4, 2))
  )
  inter_variable <- rvinecopulib::rvine_structure(
    order = c(3, 1, 2),
    struct_array = list(c(2, 2), 1)
  )

  merged <- merge_edges_fixed_full(spatial, inter_variable, bridge_var = 4)

  expect_equal(merged$trunc_lvl, 2)
  expect_equal(sort(merged$order), seq_len(merged$d))
  # Sending the output through the upstream edge conversion again validates the
  # completed forest without assuming a particular matrix orientation.
  expect_no_error(.merge_with_maps(
    list(merged), list(seq_len(merged$d)), merged$d
  ))
})

test_that("degenerate structures and invalid component lists fail explicitly", {
  expect_error(
    merge_edges_fixed_full(rvinecopulib::rvine_structure(1), rvinecopulib::rvine_structure(2), 1),
    "non-degenerate"
  )
  expect_error(
    merge_edges_individual_full(
      list(rvinecopulib::rvine_structure(
        order = 1:2, struct_array = list(2)
      )),
      rvinecopulib::rvine_structure(
        order = 1:3, struct_array = list(c(2, 3))
      ),
      1
    ),
    "one spatial structure"
  )
})

test_that("invalid local-to-global maps fail explicitly", {
  component <- rvinecopulib::rvine_structure(
    order = 1:2, struct_array = list(2)
  )

  expect_error(
    merge_rvine_structures(list(component), list(c(1, 1)), 2),
    "Invalid local-to-global"
  )
  expect_error(
    merge_rvine_structures(list(component), list(c(1, 3)), 3),
    "omit global variable 2"
  )
})

test_that("a complete inherited first tree needs no completion candidates", {
  dimension <- 1000L
  complete_first_tree <- rvinecopulib::rvine_structure(
    order = seq_len(dimension),
    struct_array = list(rep(dimension, dimension - 1L))
  )

  expect_no_error(merge_rvine_structures(
    list(complete_first_tree), list(seq_len(dimension)), dimension
  ))
})

test_that("Merge works for different bridge variables - flexible structure", {
  # Generate vine copulas for testing
  set.seed(1)
  rvs_level1 <- list(
    rvinecopulib::rvine_structure_sim(10),
    rvinecopulib::rvine_structure_sim(10),
    rvinecopulib::rvine_structure_sim(10),
    rvinecopulib::rvine_structure_sim(10),
    rvinecopulib::rvine_structure_sim(10)
  )

  rvs_level2 <- rvinecopulib::rvine_structure_sim(5)

  for (bridge_var in seq_len(rvs_level1[[1]]$d)) {
    rvs_level3 <- merge_edges_individual_full(
      rvs_level1,
      rvs_level2,
      bridge_var
    )

    expect_s3_class(rvs_level3, "rvine_structure")
    expect_equal(rvs_level3$d, rvs_level1[[1]]$d * rvs_level2$d)
  }
})
