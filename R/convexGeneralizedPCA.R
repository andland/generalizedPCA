#' @title Convex Generalized Principal Component Analysis
#'
#' @description
#' Dimensionality reduction for exponential family data by extending Pearson's
#' PCA formulation to minimize deviance. The convex relaxation
#' to projection matrices, the Fantope, is used.
#'
#' @param x matrix of either binary, proportions, count, or continuous data
#' @param k number of principal components to return
#' @param M value to approximate the saturated model
#' @param family exponential family distribution of data
#' @param weights an optional matrix of the same size as the \code{x} with non-negative weights
#' @param quiet logical; whether the calculation should give feedback
#' @param use_irlba logical; if \code{TRUE}, the function uses the irlba package
#'   to more quickly calculate the eigen-decomposition
#' @param max_iters number of maximum iterations
#' @param conv_criteria convergence criteria. The difference between average deviance
#'   in successive iterations
#' @param random_start logical; whether to randomly inititalize the parameters. If \code{FALSE},
#'   function will use an eigen-decomposition as starting value
#' @param start_H starting value for the Fantope matrix
#' @param mu main effects vector. Only used if \code{main_effects = TRUE}
#' @param main_effects logical; whether to include main effects in the model
#' @param normalize logical; whether to weight the variables to they all have equal influence
#' @param ss_factor step size multiplier. Amount by which to multiply the step size. Quadratic
#'   convergence rate can be proven for \code{ss_factor = 1}, but I have found higher values
#'   sometimes work better. The default is \code{ss_factor = 4}.
#'   If it is not converging, try \code{ss_factor = 1}.
#'
#' @return An S3 object of class \code{cgpca} which is a list with the
#' following components:
#' \item{mu}{the main effects}
#' \item{H}{a rank \code{k} Fantope matrix}
#' \item{U}{a \code{ceiling(k)}-dimentional orthonormal matrix with the loadings}
#' \item{PCs}{the princial component scores}
#' \item{M}{the parameter inputed}
#' \item{iters}{number of iterations required for convergence}
#' \item{loss_trace}{the trace of the average deviance using the Fantope matrix}
#' \item{proj_loss_trace}{the trace of the average deviance using the projection matrix}
#' \item{prop_deviance_expl}{the proportion of deviance explained by this model.
#'    If \code{main_effects = TRUE}, the null model is just the main effects, otherwise
#'    the null model estimates 0 for all natural parameters.}
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
#' # run convex generalized PCA on it
#' cgpca = convexGeneralizedPCA(mat, k = 1, M = 4, family = "binomial")
#' @export
#' @importFrom stats var
convexGeneralizedPCA <- function(x, k = 2, M = 4, family = c("gaussian", "binomial", "poisson", "multinomial"),
                                 weights, quiet = TRUE, use_irlba = FALSE, max_iters = 1000,
                                 conv_criteria = 1e-6, random_start = FALSE, start_H, mu,
                                 main_effects = TRUE, normalize = FALSE, ss_factor = 1) {
  use_irlba = use_irlba && requireNamespace("irlba", quietly = TRUE)

  family = match.arg(family)
  check_family(x, family)

  x = as.matrix(x)
  missing_mat = is.na(x)
  n = nrow(x)
  d = ncol(x)
  ones = rep(1, n)

  if (missing(weights) || length(weights) == 1) {
    weights = 1.0
    sum_weights = sum(!is.na(x))
  } else {
    if (!all(dim(weights) == dim(x))) {
      stop("x and weights are not the same dimension")
    }
    weights[is.na(x)] <- 0
    if (any(is.na(weights))) {
      stop("Can't have NA in weights")
    }
    if (any(weights < 0)) {
      stop("weights must be non-negative")
    }
    sum_weights = sum(weights)
  }

  if (main_effects) {
    if (missing(mu)) {
      if (length(weights) == 1) {
        weighted_col_means = colMeans(x, na.rm = TRUE)
      } else {
        weighted_col_means = colSums(x * weights, na.rm = TRUE) / colSums(weights)
      }
      if (any(apply(x, 2, stats::var, na.rm = TRUE) == 0)) {
        stop("At least one variable has variance of 0")
      }
      if (family == "multinomial") {
        mu = as.numeric(saturated_natural_parameters(matrix(weighted_col_means, nrow = 1), family, M))
      } else {
        mu = as.numeric(saturated_natural_parameters(weighted_col_means, family, M))
      }
    }
  } else {
    mu = rep(0, d)
  }

  eta = saturated_natural_parameters(x, family, M = M) + missing_mat * outer(ones, mu)
  eta_centered = scale(eta, center = mu, scale = FALSE)
  mu_mat = outer(rep(1, n), mu)

  if (!missing(start_H)) {
    HU = project.Fantope(start_H, k)
    H = HU$H
  } else if (random_start) {
    U = matrix(rnorm(d * d), d, d)
    U = qr.Q(qr(U))
    HU = project.Fantope(U %*% t(U), k)
    H = HU$H
  } else {
    if (use_irlba) {
      udv = irlba::irlba(scale(eta, center = mu, scale = normalize), nu = k, nv = k)
    } else {
      udv = svd(scale(eta, center = mu, scale = normalize), nu = k, nv = k)
    }
    HU = project.Fantope(udv$v[, 1:k] %*% t(udv$v[, 1:k]), k)
    H = HU$H
  }

  eta_sat_nat = saturated_natural_parameters(x, family, M = Inf)
  sat_loglike = exp_fam_log_like(x, eta_sat_nat, family, weights)

  # when x is missing eta = mu. So eta_centered is 0
  eta_centered[missing_mat] <- 0

  theta = mu_mat + eta_centered %*% H

  # Initial step size
  if (family %in% c("binomial", "multinomial")) {
    # TODO: I think this works with proportions, but need to check
    init_ss = 2 / (sum(eta_centered^2) * max(weights))
  } else if (family == "gaussian") {
    init_ss = 0.5 / (sum(eta_centered^2) * max(weights))
  } else {
    second_dir = exp_fam_variance(theta, family, weights)
    init_ss = mean((t(second_dir) %*% (eta_centered^2))^(-1)) / 2
  }
  init_ss = init_ss * ss_factor

  # only sum over non-missing x. Equivalent to replacing missing x with 0
  x_zeros = x
  x_zeros[missing_mat] <- 0

  etatX = t(eta_centered) %*% (weights * x_zeros)

  loglike = exp_fam_log_like(x, theta, family, weights)
  min_loss = 2 * (sat_loglike - loglike) / sum_weights
  best_HU = HU
  best_loglike = loglike
  if (!quiet) {
    cat(0,"  ", min_loss, "\n")
  }

  loss_trace <- proj_loss_trace <- numeric(max_iters + 1)
  loss_trace[1] <- proj_loss_trace[1] <- min_loss

  H_lag = H
  for (m in 1:max_iters) {
    if (family == "poisson") {
      step = init_ss / m
    } else {
      step = init_ss
    }

    first_dir = exp_fam_mean(theta, family)

    first_dir[missing_mat] <- 0
    etat_dir = t(eta_centered) %*% (first_dir * weights)
    deriv = etatX - etat_dir
    deriv = deriv + t(deriv) - diag(diag(deriv))

    H = H + step * deriv
    HU = project.Fantope(H, k)
    H = HU$H

    theta = mu_mat + eta_centered %*% H
    loglike = exp_fam_log_like(x, theta, family, weights)
    loss_trace[m + 1] = 2 * (sat_loglike - loglike) / sum_weights

    proj_theta = mu_mat + eta_centered %*% tcrossprod(HU$U)
    proj_loglike = exp_fam_log_like(x, proj_theta, family, weights)
    proj_loss_trace[m + 1] = 2 * (sat_loglike - proj_loglike) / sum_weights

    if (!quiet) {
      cat(m, "  ", loss_trace[m + 1], "  ", proj_loss_trace[m + 1], "\n")
    }
    if (loss_trace[m + 1] < min_loss) {
      min_loss = loss_trace[m + 1]
      best_HU = HU
      best_loglike = loglike
    }
    if (abs(loss_trace[m+1]-loss_trace[m]) < 2 * conv_criteria | min_loss == 0) {
      break
    }
  }

  # calculate the null log likelihood for % deviance explained
  # assumes no missing data
  # if (main_effects) {
  #   null_proportions = x_bar
  # } else {
  #   null_proportions = rep(0.5, d)
  # }
  # null_loglikes <- null_proportions * log(null_proportions) +
  #   (1 - null_proportions) * log(1 - null_proportions)
  # null_loglike = sum((null_loglikes * colSums(q!=0))[!(null_proportions %in% c(0, 1))])
  null_loglike = exp_fam_log_like(x, mu_mat, family, weights)

  object = list(mu = mu,
                H = best_HU$H,
                U = best_HU$U,
                PCs = eta_centered %*% best_HU$U,
                M = M,
                family = family,
                iters = m,
                loss_trace = loss_trace[1:(m + 1)],
                proj_loss_trace = proj_loss_trace[1:(m + 1)],
                prop_deviance_expl = 1 - (best_loglike - sat_loglike) / (null_loglike - sat_loglike)
  )
  class(object) <- "cgpca"
  return(object)
}

#' @title Project onto the Fantope
#'
#' @description
#' Project a symmetric matrix onto the convex set of the rank k Fantope
#'
#' @param x a symmetric matrix
#' @param k the rank of the Fantope desired
#'
#' @return
#' \item{H}{a rank \code{k} Fantope matrix}
#' \item{U}{a \code{k}-dimentional orthonormal matrix with the first \code{k} eigenvectors of \code{H}}
#' @export
project.Fantope <- function(x, k) {
  eig = eigen(x, symmetric = TRUE)
  vals = eig$values
  lower = vals[length(vals)] - k / length(vals)
  upper = max(vals)
  while(TRUE) {
    theta = (lower+upper) / 2
    sum.eig.vals = sum(pmin(pmax(vals - theta, 0), 1))
    if (abs(sum.eig.vals-k) < 1e-10) {
      break
    } else if (sum.eig.vals>k) {
      lower = theta
    } else {
      upper = theta
    }
  }
  vals = pmin(pmax(vals - theta, 0), 1)
  return(list(H = eig$vectors %*% diag(vals) %*% t(eig$vectors),
              U = matrix(eig$vectors[, 1:ceiling(k)], nrow(x), ceiling(k))))
}

#' @title Predict Convex Generalized PCA scores or reconstruction on new data
#'
#' @description Predict Convex Generalized PCA scores or reconstruction on new data
#'
#' @param object convex generalized PCA object
#' @param newdata matrix with all binary entries. If missing, will use the
#'  data that \code{object} was fit on
#' @param type the type of fitting required. \code{type = "PCs"} gives the PC scores,
#'  \code{type = "link"} gives matrix on the logit scale and \code{type = "response"}
#'  gives matrix on the probability scale
#' @param ... Additional arguments
#' @examples
#' # construct a low rank matrices in the logit scale
#' rows = 100
#' cols = 10
#' set.seed(1)
#' loadings = rnorm(cols)
#' mat_logit = outer(rnorm(rows), loadings)
#' mat_logit_new = outer(rnorm(rows), loadings)
#'
#' # convert to a binary matrix
#' mat = (matrix(runif(rows * cols), rows, cols) <= inv.logit.mat(mat_logit)) * 1.0
#' mat_new = (matrix(runif(rows * cols), rows, cols) <= inv.logit.mat(mat_logit_new)) * 1.0
#'
#' # run generalized PCA on it
#' cgpca = convexGeneralizedPCA(mat, k = 1, M = 4, family = "binomial")
#'
#' PCs = predict(cgpca, mat_new)
#' @export
predict.cgpca <- function(object, newdata, type = c("PCs", "link", "response"), ...) {
  type = match.arg(type)

  if (type == "PCs") {
    if (missing(newdata)) {
      PCs = object$PCs
    } else {
      eta = ((as.matrix(newdata) * 2) - 1) * object$M
      eta_centered = scale(eta, center = object$mu, scale = FALSE)
      eta_centered[is.na(newdata)] <- 0
      PCs = eta_centered %*% object$U
    }
    return(PCs)
  } else {
    eta = ((as.matrix(newdata) * 2) - 1) * object$M
    eta_centered = scale(eta, center = object$mu, scale = FALSE)
    eta_centered[is.na(newdata)] <- 0
    theta = outer(rep(1, nrow(eta)), object$mu) + eta_centered %*% object$H
    if (type == "link") {
      return(theta)
    } else {
      return(inv.logit.mat(theta))
    }
  }
}

#' @title Plot convex generalized PCA
#'
#' @description
#' Plots the results of a convex generalized PCA
#'
#' @param x convex generalized PCA object
#' @param type the type of plot \code{type = "trace"} plots the algorithms progress by
#' iteration, \code{type = "loadings"} plots the first 2 PC loadings,
#' \code{type = "scores"} plots the first 2 PC scores
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
#' # run convex generalized PCA on it
#' cgpca = convexGeneralizedPCA(mat, k = 1, M = 4, family = "binomial")
#'
#' \dontrun{
#' plot(cgpca)
#' }
#' @export
plot.cgpca <- function(x, type = c("trace", "loadings", "scores"), ...) {
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
print.cgpca <- function(x, ...) {
  cat(nrow(x$PCs), "rows and ")
  cat(nrow(x$H), "columns\n")
  cat("Rank", ncol(x$U), "Fantope solution with M =", x$M, "\n")
  cat("\n")
  cat(round(x$prop_deviance_expl * 100, 1), "% of deviance explained\n", sep = "")
  cat(x$iters, "iterations to converge\n")

  invisible(x)
}

#' @title CV for convex generalized PCA
#'
#' @description
#' Run cross validation on dimension and \code{M} for convex generalized PCA
#'
#' @param x matrix with all binary entries
#' @param ks the different dimensions \code{k} to try
#' @param Ms the different approximations to the saturated model \code{M} to try
#' @param folds if \code{folds} is a scalar, then it is the number of folds. If
#'  it is a vector, it should be the same length as the number of rows in \code{x}
#' @param quiet logical; whether the function should display progress
#' @param ... Additional arguments passed to convexGeneralizedPCA
#'
#' @return A matrix of the CV log likelihood with \code{k} in rows and
#'  \code{M} in columns
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
#' loglikes = cv.cgpca(mat, ks = 1:9, Ms = 3:6)
#' plot(loglikes)
#' }
#' @export
cv.cgpca <- function(x, ks, Ms = seq(2, 10, by = 2), folds = 5, quiet = TRUE, ...) {
  # TODO: does not support weights
  q = 2 * as.matrix(x) - 1
  q[is.na(q)] <- 0

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
        clpca = convexGeneralizedPCA(x[c != cv, ], k = k, M = M, ...)
        pred_theta = predict(clpca, newdat = x[c == cv, ], type = "link")
        log_likes[k == ks, M == Ms] = log_likes[k == ks, M == Ms] +
          log_like_Bernoulli(q = q[c == cv, ], theta = pred_theta)
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
