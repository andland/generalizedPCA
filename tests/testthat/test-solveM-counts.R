context("Solve M")

# construct a low rank matrix in the natural parameter scale
rows = 100
cols = 10
k = 2
set.seed(1)

loadings = rnorm(cols)
mat_np1 = outer(rnorm(rows), loadings)
mat_np2 = outer(rnorm(rows), loadings)

mat = matrix(rpois(rows * cols, c(exp(mat_np1))), rows, cols)
mat_test = matrix(rpois(rows * cols, c(exp(mat_np1))), rows, cols)

gpca = generalizedPCA(mat, k = k, M = 0, family = "poisson", main_effects = TRUE)
gpca_val = generalizedPCA(mat, k = k, M = 0, family = "poisson", main_effects = TRUE, validation = mat_test)

pred1 = predict(gpca, mat)
pred1l = predict(gpca, mat, type = "link")
pred1r = predict(gpca, mat, type = "response")
fit1l = fitted(gpca, type = "link")
fit1r = fitted(gpca, type = "response")

test_that("k = 2 poisson", {
  expect_equal(dim(gpca$U), c(cols, 2))
  expect_equal(dim(gpca$PCs), c(rows, 2))
  expect_equal(length(gpca$mu), cols)

  expect_equal(dim(pred1), c(rows, 2))
  expect_equal(dim(pred1l), c(rows, cols))
  expect_equal(dim(pred1r), c(rows, cols))
  expect_equal(dim(fit1l), c(rows, cols))
  expect_equal(dim(fit1r), c(rows, cols))
})

rm(gpca, gpca_val)

tot = rowSums(mat) + rbinom(rows, 1, 0.5)
tot_test = rowSums(mat_test) + rbinom(rows, 1, 0.5)
matp = sweep(mat, 1, tot, "/")
matp_test = sweep(mat_test, 1, tot_test, "/")
wghts = outer(tot, rep(1, cols))
wghts_test = outer(tot_test, rep(1, cols))

gpca = generalizedPCA(matp, k = k, M = 0, family = "multinomial", weights = wghts)
gpca_val = generalizedPCA(matp, k = k, M = 0, family = "multinomial", weights = wghts,
                          validation = matp_test, val_weights = wghts_test)

pred1 = predict(gpca, matp)
pred1l = predict(gpca, matp, type = "link")
pred1r = predict(gpca, matp, type = "response")
fit1l = fitted(gpca, type = "link")
fit1r = fitted(gpca, type = "response")

test_that("k = 2 multinomial", {
  expect_equal(dim(gpca$U), c(cols, 2))
  expect_equal(dim(gpca$PCs), c(rows, 2))
  expect_equal(length(gpca$mu), cols)

  expect_equal(dim(pred1), c(rows, 2))
  expect_equal(dim(pred1l), c(rows, cols))
  expect_equal(dim(pred1r), c(rows, cols))
  expect_equal(dim(fit1l), c(rows, cols))
  expect_equal(dim(fit1r), c(rows, cols))
})
