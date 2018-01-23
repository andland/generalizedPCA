#' @title Generalized Principal Component Analysis
#'
#' @description
#' Dimension reduction for exponential family data by extending Pearson's
#' PCA formulation
#'
#' @param x matrix of either binary, proportions, count, or continuous data
#' @param k number of principal components to return
#' @param M value to approximate the saturated model
#' @param family exponential family distribution of data
#' @param weights an optional matrix of the same size as the \code{x} with data weights
#' @param quiet logical; whether the calculation should give feedback
#' @param majorizer how to majorize the deviance. \code{"row"} gives
#'  tighter majorization, but may take longer to calculate each iteration.
#'  \code{"all"} may be faster per iteration, but take more iterations
#' @param partial_decomp logical; if \code{TRUE}, the function uses the RSpectra package
#'   to more quickly calculate the SVD. When the number of columns is small,
#'   the approximation may be less accurate and slower
#' @param max_iters number of maximum iterations
#' @param conv_criteria convergence criteria. The difference between average deviance
#'   in successive iterations
#' @param random_start logical; whether to randomly inititalize the parameters. If \code{FALSE},
#'   function will use an eigen-decomposition as starting value
#' @param start_U starting value for the orthogonal matrix
#' @param start_mu starting value for mu. Only used if \code{main_effects = TRUE}
#' @param main_effects logical; whether to include main effects in the model
#' @param normalize logical; whether to weight the variables to they all have equal influence
#' @param validation a validation dataset to select \code{m} with
#' @param val_weights weights associated with validation data
#'
#' @return An S3 object of class \code{gpca} which is a list with the
#' following components:
#' \item{mu}{the main effects}
#' \item{U}{a \code{k}-dimentional orthonormal matrix with the loadings}
#' \item{PCs}{the princial component scores}
#' \item{M}{the parameter inputed}
#' \item{family}{the exponential family used}
#' \item{iters}{number of iterations required for convergence}
#' \item{loss_trace}{the trace of the average deviance of the algorithm.
#'    Should be non-increasing}
#' \item{prop_deviance_expl}{the proportion of deviance explained by this model.
#'    If \code{main_effects = TRUE}, the null model is just the main effects, otherwise
#'    the null model estimates 0 for all natural parameters.}
#' @examples
#' # construct a low rank matrix in the natural parameter space
#' rows = 100
#' cols = 10
#' set.seed(1)
#' mat_np = outer(rnorm(rows), rnorm(cols))
#'
#' # generate a count matrix
#' mat = matrix(rpois(rows * cols, c(exp(mat_np))), rows, cols)
#'
#' # run Poisson PCA on it
#' gpca = generalizedPCA(mat, k = 1, M = 4, family = "poisson")
#' @export
#' @importFrom stats var
generalizedPCA <- function(x, k = 2, M = 4, family = c("gaussian", "binomial", "poisson", "multinomial"),
                           weights, quiet = TRUE, majorizer = c("row", "all"),
                           partial_decomp = FALSE, max_iters = 1000, conv_criteria = 1e-5,
                           random_start = FALSE, start_U, start_mu, main_effects = TRUE,
                           normalize = FALSE, validation, val_weights) {
  family = match.arg(family)
  check_family(x, family)

  majorizer = match.arg(majorizer)

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
    null_theta = as.numeric(saturated_natural_parameters(matrix(weighted_col_means, 1), family, M))
  } else {
    null_theta = rep(0, d)
  }

  if (normalize) {
    if (any(apply(x, 2, stats::var, na.rm = TRUE) == 0)) {
      stop("At least one variable has variance of 0. Cannot normalize")
    }

    eta_sat_nat = saturated_natural_parameters(x, family, M = Inf)

    norms = sapply(1:d, function(j) {
      2 * (exp_fam_log_like(x[, j], eta_sat_nat[, j], family, weights) -
             exp_fam_log_like(x[, j], rep(null_theta[j], n), family, weights))
      }) / n
    if (any(norms <= 0)) {
      stop("Normalization caused weights to be <= 0")
    }

    if (length(weights) == 1) {
      weights = outer(ones, 1 / norms)
    } else {
      weights = sweep(weights, 2, 1 / norms, "*")
    }
  }

  if (M == 0) {
    if (any(is.na(x))) {
      stop("Cannot solve for M with missing weights")
    }
    M = 4
    solve_M = TRUE
    if (!missing(validation)) {
      if (ncol(validation) != ncol(x)) {
        stop("validation does not have the same variables as x")
      }
      if (missing(val_weights)) {
        val_weights = 1.0
      } else {
        if (!all(dim(val_weights) == dim(validation))) {
          stop("validation and val_weights are not the same dimension")
        }
      }
      validation = as.matrix(validation)
      M_mat = exp_fam_sat_ind_mat(validation, family)
    } else {
      M_mat = exp_fam_sat_ind_mat(x, family)
    }
  } else {
    solve_M = FALSE
  }

  # if it is standard PCA, only need 1 iteration
  if (family == "gaussian" & all(weights == 1) & sum(missing_mat) == 0) {
    max_iters = 0
  }

  # Initialize #
  ##################
  if (main_effects) {
    if (!missing(start_mu)) {
      mu = start_mu
    } else {
      eta = saturated_natural_parameters(x, family, M = M)
      is.na(eta[is.na(x)]) <- TRUE
      mu = colMeans(eta, na.rm = TRUE)
      # mu = saturated_natural_parameters(colMeans(x, na.rm = TRUE), family, M)
    }
  } else {
    mu = rep(0, d)
  }

  eta = saturated_natural_parameters(x, family, M = M) + missing_mat * outer(ones, mu)
  eta_centered = scale(eta, center = mu, scale = FALSE)

  if (!missing(start_U)) {
    U = sweep(start_U, 2, sqrt(colSums(start_U^2)), "/")
  } else if (random_start) {
    U = matrix(rnorm(d * k), d, k)
    U = qr.Q(qr(U))
  } else {
    if (partial_decomp) {
      udv = RSpectra::svds(scale(eta, center = mu, scale = normalize), k, nu = k, nv = k)
    } else {
      udv = svd(scale(eta, center = mu, scale = normalize))
    }
    U = matrix(udv$v[, 1:k], d, k)
  }

  eta_sat_nat = saturated_natural_parameters(x, family, M = Inf)
  sat_loglike = exp_fam_log_like(x, eta_sat_nat, family, weights)

  loss_trace = numeric(max_iters + 1)
  theta = outer(ones, mu) + eta_centered %*% tcrossprod(U)
  loglike <- exp_fam_log_like(x, theta, family, weights)
  loss_trace[1] = 2 * (sat_loglike - loglike) / sum_weights
  ptm <- proc.time()

  if (!quiet) {
    cat(0, "  ", loss_trace[1], "")
    cat("0 hours elapsed\n")
  }

  for (m in seq_len(max_iters)) {
    last_U = U
    last_M = M
    last_mu = mu

    # TODO: incorporate missing data
    if (solve_M) {
      gpca_obj = structure(list(mu = mu, U = U, M = M, family = family),
                             class = "gpca")
      if (missing(validation)) {
        fitted_theta = predict(gpca_obj, newdata = x, type = "link")
      } else {
        fitted_theta = predict(gpca_obj, newdata = validation, type = "link")
      }
      fitted_mean = exp_fam_mean(fitted_theta, family)

      if (missing(validation)) {
        M_slope = sum(((fitted_mean - x) * weights * (M_mat %*% tcrossprod(U)))[!is.na(M_mat)])
        fitted_var = exp_fam_variance(fitted_theta, family, weights)
      } else {
        M_slope = sum(((fitted_mean - validation) * val_weights * (M_mat %*% tcrossprod(U)))[!is.na(M_mat)])
        fitted_var = exp_fam_variance(fitted_theta, family, val_weights)
      }
      M_curve = sum((fitted_var * (M_mat %*% tcrossprod(U))^2)[!is.na(M_mat)])

      M = max(M - M_slope / M_curve, 0)

      eta = saturated_natural_parameters(x, family, M = M) + missing_mat * outer(ones, mu)
      eta_centered = scale(eta, center = mu, scale = FALSE)
      theta = outer(ones, mu) + eta_centered %*% tcrossprod(U)
    }

    first_dir = exp_fam_mean(theta, family)
    second_dir = exp_fam_variance(theta, family, weights)
    if (majorizer == "row") {
      W = apply(second_dir, 1, max)
    } else if (majorizer == "all") {
      W = rep(max(second_dir), n)
    }

    # EM style estimate of Z with theta when missing data
    Z = as.matrix(theta + weights * (x - first_dir) / outer(W, rep(1, d)))
    Z[is.na(x)] <- theta[is.na(x)]
    if (main_effects) {
      mu = as.numeric(colSums((Z - eta %*% tcrossprod(U)) * W) / sum(W))
    }

    eta = saturated_natural_parameters(x, family, M = M) + missing_mat * outer(ones, mu)
    eta_centered = scale(eta, center = mu, scale = FALSE)

    mat_temp = t(eta_centered * W) %*% scale(Z, center = mu, scale = FALSE)
    mat_temp = mat_temp + t(mat_temp) -
      t(eta_centered * W) %*% eta_centered

    # RSpectra could give poor estimates of e-vectors
    # so I switch to standard eigen if it does
    repeat {
      if (partial_decomp) {
        eig = RSpectra::eigs_sym(mat_temp, k = min(k + 2, d))
      } else {
        eig = eigen(mat_temp, symmetric = TRUE)
      }
      U = matrix(eig$vectors[, 1:k], d, k)

      theta = outer(ones, mu) + eta_centered %*% tcrossprod(U)
      this_loglike <- exp_fam_log_like(x, theta, family, weights)

      if (!partial_decomp | this_loglike>=loglike) {
        loglike = this_loglike
        break
      } else {
        partial_decomp = FALSE
        if (!quiet) {
          cat("RSpectra::eigs_sym was too inaccurate in iteration ", m ,
              ". Switched to base::eigen")
        }
      }
    }

    loss_trace[m + 1] = 2 * (sat_loglike - loglike) / sum_weights

    if (!quiet) {
      time_elapsed = as.numeric(proc.time() - ptm)[3]
      tot_time = max_iters / m * time_elapsed
      time_remain = tot_time - time_elapsed
      cat(m, "  ", loss_trace[m + 1], "")
      cat(round(time_elapsed / 3600, 1), "hours elapsed. Max", round(time_remain / 3600, 1), "hours remain.\n")
    }
    if (m > 4) {
      if (abs(loss_trace[m] - loss_trace[m+1]) < 2 * conv_criteria) {
        break
      }
    }
  }

  # test if loss function increases
  if (max_iters > 0 && (loss_trace[m + 1] - loss_trace[m]) > (1e-10)) {
    U = last_U
    mu = last_mu
    M = last_M
    m = m - 1

    eta = saturated_natural_parameters(x, family, M = M) + missing_mat * outer(ones, mu)
    eta_centered = scale(eta, center = mu, scale = FALSE)

    if (family != "poisson") {
      # maybe possible with missing data? TODO: look into
      warning("Deviance increased in last iteration.\nThis should not happen!")
    } else {
      message("Deviance increased in last iteration.")
    }
  } else if (max_iters == 0) {
    m = 0
  }

  null_loglike = exp_fam_log_like(x, outer(ones, null_theta), family, weights)

  object <- list(mu = mu,
                 U = U,
                 PCs = eta_centered %*% U,
                 M = M,
                 family = family,
                 iters = m,
                 loss_trace = loss_trace[1:(m + 1)],
                 prop_deviance_expl = 1 - (loglike - sat_loglike) / (null_loglike - sat_loglike)
  )
  class(object) <- "gpca"
  object
}


#' @title Predict generalized PCA scores or reconstruction on new data
#'
#' @description Predict generalized PCA scores or reconstruction on new data
#'
#' @param object generalized PCA object
#' @param newdata matrix of the same exponential family as in \code{object}.
#'  If missing, will use the data that \code{object} was fit on
#' @param type the type of fitting required. \code{type = "PCs"} gives the PC scores,
#'  \code{type = "link"} gives matrix on the natural parameter scale and
#'  \code{type = "response"} gives matrix on the response scale
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
#' gpca = generalizedPCA(mat, k = 1, M = 4, family = "poisson")
#'
#' PCs = predict(gpca, mat_new)
#' @export
predict.gpca <- function(object, newdata, type = c("PCs", "link", "response"), ...) {
  type = match.arg(type)

  if (missing(newdata)) {
    PCs = object$PCs
  } else {
    check_family(newdata, object$family)

    eta = saturated_natural_parameters(newdata, object$family, object$M) +
      is.na(newdata) * outer(rep(1, nrow(newdata)), object$mu)
    PCs = scale(eta, center = object$mu, scale = FALSE) %*% object$U
  }

  if (type == "PCs") {
    PCs
  } else {
    object$PCs = PCs
    fitted(object, type, ...)
  }
}

#' @title Fitted values using generalized PCA
#'
#' @description
#' Fit a lower dimentional representation of the exponential family matrix using generalized PCA
#'
#' @param object generalized PCA object
#' @param type the type of fitting required. \code{type = "link"} gives output on the natural
#'  parameter scale and \code{type = "response"} gives output on the response scale
#' @param ... Additional arguments
#' @examples
#' # construct a low rank matrix in the natural parameter space
#' rows = 100
#' cols = 10
#' set.seed(1)
#' mat_np = outer(rnorm(rows), rnorm(cols))
#'
#' # generate a count matrix
#' mat = matrix(rpois(rows * cols, c(exp(mat_np))), rows, cols)
#'
#' # run Poisson PCA on it
#' gpca = generalizedPCA(mat, k = 1, M = 4, family = "poisson")
#'
#' # construct fitted expected value of counts matrix
#' fit = fitted(gpca, type = "response")
#' @export
fitted.gpca <- function(object, type = c("link", "response"), ...) {
  type = match.arg(type)
  n = nrow(object$PCs)

  theta = outer(rep(1, n), object$mu) + tcrossprod(object$PCs, object$U)

  if (type == "link") {
    return(theta)
  } else if (type == "response") {
    return(exp_fam_mean(theta, object$family))
  }
}

#' @title Plot generalized PCA
#'
#' @description
#' Plots the results of a generalized PCA
#'
#' @param x generalized PCA object
#' @param type the type of plot \code{type = "trace"} plots the algorithms progress by
#' iteration, \code{type = "loadings"} plots the first 2 principal component
#' loadings, \code{type = "scores"} plots the loadings first 2 principal component scores
#' @param ... Additional arguments
#' @examples
#' # construct a low rank matrix in the natural parameter space
#' rows = 100
#' cols = 10
#' set.seed(1)
#' mat_np = outer(rnorm(rows), rnorm(cols))
#'
#' # generate a count matrix
#' mat = matrix(rpois(rows * cols, c(exp(mat_np))), rows, cols)
#'
#' # run logistic PCA on it
#' gpca = generalizedPCA(mat, k = 2, M = 4, family = "poisson")
#'
#' \dontrun{
#' plot(gpca)
#' }
#' @export
plot.gpca <- function(x, type = c("trace", "loadings", "scores"), ...) {
  type = match.arg(type)

  if (type == "trace") {
    df = data.frame(Iteration = 0:x$iters,
                    Deviance = x$loss_trace)
    p <- ggplot2::ggplot(df, ggplot2::aes_string("Iteration", "Deviance")) +
      ggplot2::geom_line()
  } else if (type == "loadings") {
    df = data.frame(x$U)
    colnames(df) <- paste0("PC", 1:ncol(df))
    if (ncol(df) == 1) {
      df$PC2 = 0
      p <- ggplot2::ggplot(df, ggplot2::aes_string("PC1", "PC2")) + ggplot2::geom_point() +
        ggplot2::labs(y = NULL)
    } else {
      p <- ggplot2::ggplot(df, ggplot2::aes_string("PC1", "PC2")) + ggplot2::geom_point()
    }
  } else if (type == "scores") {
    df = data.frame(x$PCs)
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
print.gpca <- function(x, ...) {
  cat(nrow(x$PCs), "rows and ")
  cat(nrow(x$U), "columns of ")
  cat(x$family, "data\n")
  cat("Rank", ncol(x$U), "solution with M =", x$M, "\n")
  cat("\n")
  cat(round(x$prop_deviance_expl * 100, 1), "% of deviance explained\n", sep = "")
  cat(x$iters, "iterations to converge\n")

  invisible(x)
}

#' @title CV for generalized PCA
#'
#' @description
#' Run cross validation on dimension and \code{M} for generalized PCA
#'
#' @param x matrix of either binary, count, or continuous data
#' @param ks the different dimensions \code{k} to try
#' @param Ms the different approximations to the saturated model \code{M} to try
#' @param family exponential family distribution of data
#' @param weights an optional matrix of the same size as the \code{x} with data weights
#' @param folds if \code{folds} is a scalar, then it is the number of folds. If
#'  it is a vector, it should be the same length as the number of rows in \code{x}
#' @param quiet logical; whether the function should display progress
#' @param ... Additional arguments passed to \code{generalizedPCA}
#'
#' @return A matrix of the CV log likelihood with \code{k} in rows and
#'  \code{M} in columns
#' @examples
#' # construct a low rank matrix in the natural parameter space
#' rows = 100
#' cols = 10
#' set.seed(1)
#' mat_np = outer(rnorm(rows), rnorm(cols))
#'
#' # generate a count matrix
#' mat = matrix(rpois(rows * cols, c(exp(mat_np))), rows, cols)
#'
#' \dontrun{
#' loglikes = cv.gpca(mat, ks = 1:9, Ms = 3:6, family = "poisson", quiet = FALSE)
#' plot(loglikes)
#' }
#' @export
cv.gpca <- function(x, ks, Ms = seq(2, 10, by = 2), family = c("gaussian", "binomial", "poisson", "multinomial"),
                    weights, folds = 5, quiet = TRUE, ...) {
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
    cv = sample(1:folds, nrow(x), replace = TRUE)
  }

  if (missing(weights) || length(weights) == 1) {
    weights = matrix(1.0, nrow(x), ncol(x))
  }

  log_likes = matrix(0, length(ks), length(Ms),
                     dimnames = list(k = ks, M = Ms))
  for (k in ks) {
    for (M in Ms) {
      if (!quiet) {
        cat("k =", k, "M =", M, "")
      }
      for (c in unique(cv)) {
        if (!quiet) {
          cat(".")
        }
        gpca = generalizedPCA(x[c != cv, ], k = k, M = M, family = family, weights = weights[c != cv, ], ...)
        pred_theta = predict(gpca, newdat = x[c == cv, ], type = "link")
        log_likes[k == ks, M == Ms] = log_likes[k == ks, M == Ms] +
          exp_fam_log_like(x = x[c == cv, ], theta = pred_theta, family = family, weights = weights[c == cv, ])
      }
      if (!quiet) {
        cat("", log_likes[k == ks, M == Ms], "\n")
      }
    }
  }
  class(log_likes) <- c("matrix", "cv.gpca")
  which_min = which(log_likes == max(log_likes), arr.ind = TRUE)
  if (!quiet) {
    cat("Best: k =", ks[which_min[1]], "M =", Ms[which_min[2]], "\n")
  }

  return(log_likes)
}

#' @title Plot CV for generalized PCA
#'
#' @description
#' Plot cross validation results generalized PCA
#'
#' @param x a \code{cv.gpca} object
#' @param ... Additional arguments
#' @examples
#' # construct a low rank matrix in the natural parameter space
#' rows = 100
#' cols = 10
#' set.seed(1)
#' mat_np = outer(rnorm(rows), rnorm(cols))
#'
#' # generate a count matrix
#' mat = matrix(rpois(rows * cols, c(exp(mat_np))), rows, cols)
#'
#' \dontrun{
#' loglikes = cv.gpca(mat, ks = 1:9, Ms = 3:6, family = "poisson")
#' plot(loglikes)
#' }
#' @export
#' @importFrom utils type.convert
plot.cv.gpca <- function(x, ...) {
  # replaces reshape2::melt(-x, value.name = "NegLogLikelihood")
  Ms = utils::type.convert(colnames(x))
  ks = utils::type.convert(rownames(x))
  df = data.frame(k = rep(ks, times = length(Ms)),
                  m = rep(Ms, each = length(ks)),
                  NegLogLikelihood = as.vector(-x))

  if (ncol(x) == 1) {
    df$M = factor(df$M)
    p <- ggplot2::ggplot(df, ggplot2::aes_string("k", "NegLogLikelihood", colour = "m")) +
      ggplot2::geom_line()
  } else {
    df$k = factor(df$k)
    p <- ggplot2::ggplot(df, ggplot2::aes_string("m", "NegLogLikelihood", colour = "k")) +
      ggplot2::geom_line()
  }
  return(p)
}
