#' Exponential Family Matrix Factorization
#'
#' @description Collins et al. (2001)'s Exponential Family PCA
#'
#' @param x matrix of either binary, proportions, count, or continuous data
#' @param k dimension
#' @param family exponential family distribution of data
#' @param weights an optional matrix of the same size as the \code{x} with data weights
#' @param quiet logical; whether the calculation should give feedback
#' @param max_iters maximum number of iterations
#' @param conv_criteria convergence criteria
#' @param random_start whether to randomly initialize \code{A} and \code{B}
#' @param start_A initial value for \code{A}
#' @param start_B initial value for \code{B}
#' @param mu specific value for \code{mu}, the mean vector of \code{x}
#' @param main_effects logical; whether to include main effects in the model
#' @param method which algorithm to use. \code{"als"} uses alternating least squares.
#'   It has the benefit of majozing row-wise and column-wise for each of the updates.
#'   \code{"svd"} uses singular value decomposition (similar to de Leeuw, 2006). It has to
#'   a more gereral majorization, which may not work well for heterogeneous matrices.
#'
#' @return An S3 object of class \code{gsmf} which is a list with the
#' following components:
#' \item{mu}{the main effects for dimensionality reduction}
#' \item{A}{the \code{n}x\code{k}-dimentional matrix with the scores}
#' \item{B}{the \code{d}x\code{k}-dimentional matrix with the loadings}
#' \item{beta}{the \code{k + 1} length vector of the coefficients}
#' \item{family_x}{the exponential family of covariates}
#' \item{family_y}{the exponential family of response}
#' \item{iters}{number of iterations required for convergence}
#' \item{loss_trace}{the trace of the average deviance of the algorithm.
#'    Should be non-increasing}
#' \item{prop_deviance_expl}{the proportion of deviance explained by this model.
#'    If \code{main_effects = TRUE}, the null model is just the main effects, otherwise
#'    the null model estimates 0 for all natural parameters.}
#'
#' @export
#' @importFrom RSpectra svds
#'
#' @references
#' de Leeuw, Jan, 2006. Principal component analysis of binary data
#' by iterated singular value decomposition. Computational Statistics & Data Analysis
#' 50 (1), 21--39.
#'
#' Collins, M., Dasgupta, S., & Schapire, R. E., 2001. A generalization of principal
#' components analysis to the exponential family. In NIPS, 617--624.
#'
#' @examples
#' rows = 100
#' cols = 10
#' set.seed(1)
#' mat_np = outer(rnorm(rows), rnorm(cols))
#'
#' # generate a count matrix and binary response
#' mat = matrix(rpois(rows * cols, c(exp(mat_np))), rows, cols)
#'
#' mod = generalizedMF(mat, k = 1, family = "poisson", quiet = FALSE)
#'
generalizedMF <- function(x, k = 2, family = c("gaussian", "binomial", "poisson"),
                      weights, quiet = TRUE, max_iters = 1000, conv_criteria = 1e-5,
                      random_start = FALSE, start_A,  start_B, mu, main_effects = TRUE,
                      method = c("als", "svd")) {
  family = match.arg(family)
  check_family(x, family)
  method = match.arg(method)

  if (method == "als") {
    # used in ALS to make sure the estimates do not diverge
    # a penalization like this was recommended in Rish et al. (2008), ICML
    L2_eps = 1e-5
  }

  x = as.matrix(x)
  missing_mat = is.na(x)
  n = nrow(x)
  d = ncol(x)
  ones = rep(1, n)

  if (missing(weights)) {
    weights = 1.0
    sum_weights = sum(!is.na(x))
  } else {
    weights[is.na(x)] <- 0
    if (any(is.na(weights))) {
      stop("Can't have NA in weights")
    }
    if (any(weights < 0)) {
      stop("weights must be non-negative")
    }
    if (!all(dim(weights) == dim(x))) {
      stop("x and weights are not the same dimension")
    }
    sum_weights = sum(weights)
  }

  # calculate the null log likelihood for % deviance explained and normalization
  if (main_effects) {
    if (length(weights) == 1) {
      weighted_col_means = colMeans(x, na.rm = TRUE)
    } else {
      weighted_col_means = colSums(x * weights, na.rm = TRUE) / colSums(weights)
    }
    null_theta = as.numeric(saturated_natural_parameters(matrix(weighted_col_means, 1), family, M = Inf))
  } else {
    null_theta = rep(0, d)
  }
  null_deviance = exp_fam_deviance(x, outer(ones, null_theta), family, weights) / sum_weights

  # Initialize #
  ##################
  if (main_effects) {
    if (missing(mu)) {
      mu = saturated_natural_parameters(weighted_col_means, family, Inf)
    } else {
      mu = as.numeric(mu)
      stopifnot(length(mu) == d)
    }
  } else {
    mu = rep(0, d)
  }

  if (!missing(start_B)) {
    stopifnot(dim(start_B) == c(d, k))
    B = as.matrix(start_B)
  } else if (random_start) {
    B = matrix(rnorm(d * k), d, k)
  } else {
    udv = svd(scale(saturated_natural_parameters(x, family, 4), TRUE, FALSE))
    B = udv$v[, 1:k, drop = FALSE]
  }

  if (!missing(start_A)) {
    stopifnot(dim(start_A) == c(n, k))
    A = as.matrix(start_A)
  } else if (random_start) {
    A = matrix(rnorm(n * k), n, k)
  } else {
    if (!missing(start_B)) {
      udv = svd(scale(saturated_natural_parameters(x, family, 4), TRUE, FALSE))
    }
    A = udv$u[, 1:k, drop = FALSE] %*% diag(udv$d, k, k)
  }

  loss_trace = numeric(max_iters)
  theta = outer(ones, mu) + tcrossprod(A, B)
  loss_trace[1] = exp_fam_deviance(x, theta, family, weights) / sum_weights

  ptm <- proc.time()

  if (!quiet) {
    cat(0, "  ", loss_trace[1], "")
    cat("0 hours elapsed\n")
  }

  for (ii in seq_len(max_iters)) {
    if (method == "als") {
      # update A
      theta = outer(ones, mu) + tcrossprod(A, B)
      first_dir = exp_fam_mean(theta, family)
      second_dir = exp_fam_variance(theta, family, weights)

      W = apply(second_dir, 2, max)
      Z = as.matrix(theta + weights * (x - first_dir) / outer(ones, W))
      Z[is.na(x)] <- theta[is.na(x)]

      A = t(solve(t(B) %*% diag(W) %*% B + diag(L2_eps, k, k),
                  t(B) %*% diag(W) %*% t(scale(Z, center = mu, scale = FALSE))))

      # update B
      theta = outer(ones, mu) + tcrossprod(A, B)
      first_dir = exp_fam_mean(theta, family)
      second_dir = exp_fam_variance(theta, family, weights)

      W = apply(second_dir, 1, max)
      Z = as.matrix(theta + weights * (x - first_dir) / outer(W, rep(1, d)))
      Z[is.na(x)] <- theta[is.na(x)]

      B = t(solve(t(A) %*% diag(W, n, n) %*% A + diag(L2_eps, k, k),
                  t(A) %*% diag(W, n, n) %*% scale(Z, center = mu, scale = FALSE)))
    } else if (method == "svd") {
      theta = outer(ones, mu) + tcrossprod(A, B)
      first_dir = exp_fam_mean(theta, family)
      second_dir = exp_fam_variance(theta, family, weights)

      W = max(second_dir)
      Z = as.matrix(theta + weights * (x - first_dir) / W)
      Z[is.na(x)] <- theta[is.na(x)]

      udv = svd(scale(Z, center = mu, scale = FALSE))
      A = udv$u[, 1:k, drop = FALSE] %*% diag(udv$d, k, k)
      B = udv$v[, 1:k, drop = FALSE]
    }

    # Calc Deviance
    theta = outer(ones, mu) + tcrossprod(A, B)
    loss_trace[ii] <- exp_fam_deviance(x, theta, family, weights) / sum_weights

    if (!quiet) {
      time_elapsed = as.numeric(proc.time() - ptm)[3]
      tot_time = max_iters / ii * time_elapsed
      time_remain = tot_time - time_elapsed
      cat(ii, "  ", loss_trace[ii], "")
      cat(round(time_elapsed / 3600, 1), "hours elapsed. Max", round(time_remain / 3600, 1), "hours remain.\n")
    }

    if (ii > 1 && abs(loss_trace[ii] - loss_trace[ii - 1]) < 2 * conv_criteria) {
      break
    } else {
      B_lag = B
    }
  }

  object <- list(
    mu = mu,
    A = A,
    B = B,
    family = family,
    iters = ii,
    loss_trace = loss_trace[1:ii],
    prop_deviance_expl = 1 - loss_trace[ii] / null_deviance
  )
  class(object) <- "gmf"
  object
}

#' @title Predict generalized PCA scores or reconstruction on new data
#'
#' @description Predict generalized PCA scores or reconstruction on new data
#'
#' @param object generalized MF object
#' @param newdata matrix of the same exponential family as covariates in \code{object}.
#'  If missing, will use the data that \code{object} was fit on
#' @param type the type of fitting required.
#'  \code{type = "PCs"} gives matrix of principal components of \code{x},
#'  \code{type = "link"} gives a matrix on the natural parameter scale, and
#'  \code{type = "response"} gives a matrix on the response scale
#' @param quiet logical; whether the calculation should show progress
#' @param max_iters maximum number of iterations
#' @param conv_criteria convergence criteria
#' @param start_A initial value for \code{A}
#' @param ... Additional arguments
#' @examples
#' # construct a low rank matrices in the natural parameter space
#' rows = 100
#' cols = 10
#' set.seed(1)
#' loadings = rnorm(cols)
#' mat_np = outer(rnorm(rows), rnorm(cols))
#' mat_np_new = outer(rnorm(rows), loadings)
#'
#' # generate a count matrices
#' mat = matrix(rpois(rows * cols, c(exp(mat_np))), rows, cols)
#' mat_new = matrix(rpois(rows * cols, c(exp(mat_np_new))), rows, cols)
#'
#' # run Poisson PCA on it
#' gmf = generalizedMF(mat, k = 1, family = "poisson")
#'
#' A = predict(gmf, mat_new)
#'
#' @export
predict.gmf <- function(object, newdata, type = c("PCs", "link", "response"), quiet = TRUE,
                         max_iters = 1000, conv_criteria = 1e-5, start_A,...) {
  type = match.arg(type)

  if (missing(newdata)) {
    A = object$A
  } else {
    n = nrow(newdata)
    d = ncol(newdata)
    k = ncol(object$B)
    stopifnot(d == nrow(object$B))
    check_family(newdata, object$family)
    L2_eps = 1e-5

    ones = rep(1, n)

    # solve for A
    if (missing(start_A)) {
      theta = outer(ones, object$mu)
    } else {
      stopifnot(dim(start_A) == c(n, k))
      theta = outer(ones, object$mu) + tcrossprod(start_A, object$B)
    }

    last_deviance = exp_fam_deviance(newdata, theta, object$family) / sum(!is.na(newdata))
    if (!quiet) {
      cat(0, " ", last_deviance, "\n")
    }
    for (ii in seq_len(max_iters)) {
      first_dir = exp_fam_mean(theta, object$family)
      second_dir = exp_fam_variance(theta, object$family)

      multiplier = 1
      # while (TRUE) {
        # W = apply(second_dir, 2, max)
      W = rep(max(second_dir), d)
      Z = as.matrix(theta + (newdata - first_dir) / outer(ones, W))
      Z[is.na(newdata)] <- theta[is.na(newdata)]

      A = t(solve(t(object$B) %*% diag(W, d, d) %*% object$B + diag(L2_eps, k, k),
                  t(object$B) %*% diag(W, d, d) %*% t(scale(Z, object$mu, FALSE))))

      theta = outer(ones, object$mu) + tcrossprod(A, object$B)
      this_deviance = exp_fam_deviance(newdata, theta, object$family) / sum(!is.na(newdata))

        # if (this_deviance < last_deviance) {
        #   break
        # } else {
        #   multiplier = multiplier * 2
        # }

      if (!quiet) {
        cat(ii ," ", this_deviance, "\n")
      }

      if (abs(last_deviance - this_deviance) < 2 * conv_criteria)
        break
      last_deviance = this_deviance
    }
  }

  if (type == "PCs") {
    return(A)
  } else {
    theta = outer(rep(1, nrow(A)), object$mu) + tcrossprod(A, object$B)

    if (type == "link") {
      return(theta)
    } else if (type == "response") {
      return(exp_fam_mean(theta, object$family))
    }
  }
}

#' @title Plot generalized MF
#'
#' @description
#' Plots the results of a generalized MF
#'
#' @param x generalized MF object
#' @param type the type of plot \code{type = "trace"} plots the algorithms progress by
#' iteration, \code{type = "loadings"} plots the first 2 principal component
#' loadings, \code{type = "scores"} plots the loadings first 2 principal component scores
#' @param ... Additional arguments
#' @examples
#' # construct a low rank matrix in the logit scale
#' rows = 100
#' cols = 10
#' set.seed(1)
#' mat_logit = outer(rnorm(rows), rnorm(cols))
#'
#' # generate a binary matrix
#' mat = (matrix(runif(rows * cols), rows, cols) <= inv.logit.mat(mat_logit)) * 1.0
#'
#' # run logistic SVD on it
#' gmf = generalizedMF(mat, k = 2, family = "binomial")
#'
#' \dontrun{
#' plot(gmf)
#' }
#' @export
plot.gmf <- function(x, type = c("trace", "loadings", "scores"), ...) {
  type = match.arg(type)

  if (type == "trace") {
    df = data.frame(Iteration = 1:x$iters,
                    Deviance = x$loss_trace)
    p <- ggplot2::ggplot(df, ggplot2::aes_string("Iteration", "Deviance")) +
      ggplot2::geom_line()
  } else if (type == "loadings") {
    df = data.frame(x$B)
    colnames(df) <- paste0("PC", 1:ncol(df))
    if (ncol(df) == 1) {
      df$PC2 = 0
      p <- ggplot2::ggplot(df, ggplot2::aes_string("PC1", "PC2")) + ggplot2::geom_point() +
        ggplot2::labs(y = NULL)
    } else {
      p <- ggplot2::ggplot(df, ggplot2::aes_string("PC1", "PC2")) + ggplot2::geom_point()
    }
  } else if (type == "scores") {
    df = data.frame(x$A)
    colnames(df) <- paste0("PC", 1:ncol(df))
    if (ncol(df) == 1) {
      df$PC2 = 0
      p <- ggplot2::ggplot(df, ggplot2::aes_string("PC1", "PC2")) + ggplot2::geom_point() +
        ggplot2::labs(y = NULL)
    } else {
      p <- ggplot2::ggplot(df, ggplot2::aes_string("PC1", "PC2")) + ggplot2::geom_point()
    }
  }

  return(p)
}

#' @export
print.gmf <- function(x, ...) {
  cat(nrow(x$A), "rows and ")
  cat(nrow(x$B), "columns\n")
  cat("Rank", ncol(x$B), "solution\n")
  cat("\n")
  cat(round(x$prop_deviance_expl * 100, 1), "% of deviance explained\n", sep = "")
  cat(x$iters, "iterations to converge\n")

  invisible(x)
}


#' @title CV for generalized MF
#'
#' @description
#' Run cross validation on dimension for generalized MF
#'
#' @param x matrix of either binary, count, or continuous data
#' @param ks the different dimensions \code{k} to try
#' @param family exponential family distribution of data
#' @param folds if \code{folds} is a scalar, then it is the number of folds. If
#'  it is a vector, it should be the same length as the number of rows in \code{x}
#' @param quiet logical; whether the function should display progress
#' @param ... Additional arguments passed to generalizedMF
#'
#' @return A matrix of the CV deviance with \code{k} in rows
#'
#' @examples
#' # construct a low rank matrix in the logit scale
#' rows = 100
#' cols = 10
#' set.seed(1)
#' mat_logit = outer(rnorm(rows), rnorm(cols))
#'
#' # generate a binary matrix
#' mat = (matrix(runif(rows * cols), rows, cols) <= inv.logit.mat(mat_logit)) * 1.0
#'
#' \dontrun{
#' deviances = cv.gmf(mat, ks = 1:9, family = "binomial")
#' plot(deviances)
#' }
#' @export
cv.gmf <- function(x, ks, family = c("gaussian", "binomial", "poisson", "multinomial"),
                    folds = 5, quiet = TRUE, ...) {
  family = match.arg(family)
  check_family(x, family)

  if (length(folds) > 1) {
    # does this work if factor?
    if (length(unique(folds)) <= 1) {
      stop("If inputing CV split, must be more than one level")
    }
    if (length(folds) != nrow(x)) {
      stop("if folds is a vector, it should be of same length as nrow(x)")
    }
    cv = folds
  } else {
    cv = sample(1:folds, nrow(q), replace = TRUE)
  }

  log_likes = matrix(0, length(ks), 1,
                     dimnames = list(k = ks, M = "GMF"))
  for (k in ks) {
    if (!quiet) {
      cat("k =", k, "\n")
    }
    for (c in unique(cv)) {
      gmf = generalizedMF(x[c != cv, ], k = k, family = family, ...)
      pred_theta = predict(gmf, newdat = x[c == cv, ], type = "link")
      log_likes[k == ks] = log_likes[k == ks] +
        exp_fam_deviance(x[c == cv, ], theta = pred_theta, family = family)
    }
  }
  class(log_likes) <- c("matrix", "cv.gpca")
  which_max = which.max(log_likes)
  if (!quiet) {
    cat("Best: k =", ks[which_max], "\n")
  }

  return(-log_likes)
}
