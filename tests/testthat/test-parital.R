context("Parital Decomp")

# construct a low rank matrix in the natural parameter scale
rows = 100
cols = 10
k = 1
set.seed(1)
mat_np = outer(rnorm(rows), rnorm(cols))

mat = (matrix(runif(rows * cols), rows, cols) <= inv.logit.mat(mat_np)) * 1.0

gpca = generalizedPCA(mat, k = k, M = 4, family = "binomial", main_effects = TRUE, partial_decomp = FALSE)
gmf = generalizedMF(mat, k = k, family = "binomial", main_effects = TRUE, partial_decomp = FALSE, method = "svd")

ppca = generalizedPCA(mat, k = k, M = 4, family = "poisson", main_effects = TRUE, partial_decomp = FALSE)
pmf = generalizedMF(mat, k = k, family = "poisson", main_effects = TRUE, partial_decomp = FALSE, method = "svd")

gpca_part = generalizedPCA(mat, k = k, M = 4, family = "binomial", main_effects = TRUE, partial_decomp = TRUE)
gmf_part = generalizedMF(mat, k = k, family = "binomial", main_effects = TRUE, partial_decomp = TRUE, method = "svd")

ppca_part = generalizedPCA(mat, k = k, M = 4, family = "poisson", main_effects = TRUE, partial_decomp = TRUE)
pmf_part = generalizedMF(mat, k = k, family = "poisson", main_effects = TRUE, partial_decomp = TRUE, method = "svd")

test_that("partial_decomp = full decomp", {
  expect_equal(gpca$iters, gpca_part$iters)
  expect_equal(gmf$iters, gmf_part$iters)
  expect_equal(ppca$iters, ppca_part$iters)
  expect_equal(pmf$iters, pmf_part$iters)

  expect_equal(gpca$loss_trace, gpca_part$loss_trace)
  expect_equal(gmf$loss_trace, gmf_part$loss_trace)
  expect_equal(ppca$loss_trace, ppca_part$loss_trace)
  expect_equal(pmf$loss_trace, pmf_part$loss_trace)
})
