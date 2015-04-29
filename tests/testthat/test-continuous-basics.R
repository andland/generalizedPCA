context("Continuous")

# construct a low rank matrix in the natural parameter scale
rows = 100
cols = 10
k = 1
set.seed(1)
mat_np = outer(rnorm(rows), rnorm(cols))

mat = matrix(rnorm(rows * cols, c(mat_np)), rows, cols)

gpca = generalizedPCA(mat, k = k, M = 4, family = "gaussian", main_effects = TRUE)

pca_mu = colMeans(mat)
pca = svd(scale(mat, center = pca_mu, scale = FALSE))
gpca2 = gpca
gpca$mu = pca_mu
gpca$U = matrix(pca$v[, 1:k], cols, k)


pred1 = predict(gpca, mat)
pred1l = predict(gpca, mat, type = "link")
pred1r = predict(gpca, mat, type = "response")
fit1l = fitted(gpca, type = "link")
fit1r = fitted(gpca, type = "response")

pred2l = predict(gpca2, mat, type = "link")
pred2r = predict(gpca2, mat, type = "response")

test_that("correct classes", {
  expect_is(gpca, "gpca")

  expect_is(pred1, "matrix")
  expect_is(pred1l, "matrix")
  expect_is(pred1r, "matrix")
  expect_is(fit1l, "matrix")
  expect_is(fit1r, "matrix")
})

test_that("k = 1 dimensions", {
  expect_equal(dim(gpca$U), c(cols, 1))
  expect_equal(dim(gpca$PCs), c(rows, 1))
  expect_equal(length(gpca$mu), cols)

  expect_equal(dim(pred1), c(rows, 1))
  expect_equal(dim(pred1l), c(rows, cols))
  expect_equal(dim(pred1r), c(rows, cols))
  expect_equal(dim(fit1l), c(rows, cols))
  expect_equal(dim(fit1r), c(rows, cols))
})

test_that("k = 1 same fits as standard PCA", {
  expect_equal(pred1r, pred2r)
  expect_equal(pred1l, pred2l)
  # Don't compare fitted, because they rely on $PCs
  expect_equal(gpca$mu, pca_mu)
  expect_equal(pca$v[, 1], gpca$U[, 1])
  expect_equal(pca$u[, 1] * pca$d[1], gpca$PCs[, 1])
})

rm(gpca, pred1, pred1l, pred1r, fit1l, fit1r)

k = 2
gpca = generalizedPCA(mat, k = k, M = 4, family = "gaussian", main_effects = TRUE)

pred1 = predict(gpca, mat)
pred1l = predict(gpca, mat, type = "link")
pred1r = predict(gpca, mat, type = "response")
fit1l = fitted(gpca, type = "link")
fit1r = fitted(gpca, type = "response")

test_that("k = 2 dimensions", {
  expect_equal(dim(gpca$U), c(cols, 2))
  expect_equal(dim(gpca$PCs), c(rows, 2))
  expect_equal(length(gpca$mu), cols)

  expect_equal(dim(pred1), c(rows, 2))
  expect_equal(dim(pred1l), c(rows, cols))
  expect_equal(dim(pred1r), c(rows, cols))
  expect_equal(dim(fit1l), c(rows, cols))
  expect_equal(dim(fit1r), c(rows, cols))
})
