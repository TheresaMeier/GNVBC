test_that("get_nested_vine returns expected components", {
  set.seed(1)

  data <- matrix(runif(120), ncol = 6)
  colnames(data) <- c("tas.1", "tas.2", "tas.3",
                      "pr.1",  "pr.2",  "pr.3")

  res <- get_nested_vine(
    data,
    nrows = 1,
    ncols = 3,
    fixed = TRUE,
    mask = FALSE
  )

  expect_type(res, "list")
  expect_named(
    res,
    c(
      "bridge_var",
      "rvs_level1",
      "rvs_level2",
      "rvs_level3",
      "vine_level1",
      "vine_level2",
      "vine_level3"
    ),
    ignore.order = TRUE
  )
})

test_that("get_nested_vine returns expected components", {
  set.seed(1)

  data <- matrix(runif(120), ncol = 6)
  colnames(data) <- c("tas.1", "tas.2", "tas.3",
                      "pr.1",  "pr.2",  "pr.3")

  res <- get_nested_vine(
    data,
    nrows = 1,
    ncols = 3,
    fixed = FALSE,
    mask = TRUE
  )

  expect_type(res, "list")
  expect_named(
    res,
    c(
      "bridge_var",
      "rvs_level1_tas",
      "rvs_level1_pr",
      "rvs_level2",
      "rvs_level3",
      "vine_level1_tas",
      "vine_level1_pr",
      "vine_level2",
      "vine_level3"
    ),
    ignore.order = TRUE
  )
})

test_that("results are reproducible given seed", {
  data <- matrix(runif(120), ncol = 6)
  colnames(data) <- c("tas.1", "tas.2", "tas.3",
                      "pr.1",  "pr.2",  "pr.3")

  res1 <- get_nested_vine(
    data, nrows = 1, ncols = 3,
    fixed = TRUE, seed = 42, mask = FALSE
  )

  res2 <- get_nested_vine(
    data, nrows = 1, ncols = 3,
    fixed = TRUE, seed = 42, mask = FALSE
  )

  expect_identical(
    res1$rvs_level3,
    res2$rvs_level3
  )
})

test_that("invalid bridge_var throws error", {
  data <- matrix(runif(120), ncol = 6)
  colnames(data) <- c("tas.1", "tas.2", "tas.3",
                      "pr.1",  "pr.2",  "pr.3")

  expect_error(
    get_nested_vine(
      data,
      nrows = 1,
      ncols = 3,
      bridge_var = 99
    ),
    "Bridging variable must be between"
  )
})

test_that("full model fits on minimal data", {
  data <- matrix(runif(60), ncol = 4)
  colnames(data) <- c("tas.1", "tas.2", "pr.1", "pr.2")

  res <- get_nested_vine(
    data,
    nrows = 1,
    ncols = 2,
    fixed = FALSE,
    mask = FALSE
  )

  expect_s3_class(res$vine_level3, "vinecop")
})

