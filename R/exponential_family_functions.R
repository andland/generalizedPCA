#' @title Inverse logit for matrices
#'
#' @description
#' Apply the inverse logit function to a matrix, element-wise. It
#' generalizes the \code{inv.logit} function from the \code{gtools}
#' library to matrices
#'
#' @param x matrix
#' @param min Lower end of logit interval
#' @param max Upper end of logit interval
#' @examples
#' (mat = matrix(rnorm(10 * 5), nrow = 10, ncol = 5))
#' inv.logit.mat(mat)
#' @export
inv.logit.mat <- function(x, min = 0, max = 1) {
  .Call(inv_logit_mat, x, min, max)
}

#' @title Bernoulli Log Likelihood
#'
#' @description
#' Calculate the Bernoulli log likelihood of matrix
#'
#' @param x matrix with all binary entries
#' @param theta estimated natural parameters with
#'  same dimensions as x
#' @param q instead of x, you can input matrix q which is
#'  -1 if \code{x = 0}, 1 if \code{x = 1}, and 0 if \code{is.na(x)}
#' @export
log_like_Bernoulli <- function(x, theta, q) {
  if (!missing(x)) {
    q = 2 * as.matrix(x) - 1
    q[is.na(q)] <- 0
  }
  .Call(compute_loglik, q, theta)
}


saturated_natural_parameters <- function(mat, family, M) {
  if (family == "gaussian") {
    eta = mat
  } else if (family == "binomial") {
    eta = abs(M) * (2 * mat - 1)
  } else if (family == "poisson") {
    eta = log(mat)
    eta[mat==0] = -abs(M)
  }
  return(eta)
}

exponential_family_mean <- function(theta, family) {
  if (family == "gaussian") {
    mean_mat = theta
  } else if (family == "binomial") {
    mean_mat = inv.logit.mat(theta)
    # (1 + exp(-theta))^(-1)
  } else if (family == "poisson") {
    mean_mat = exp(theta)
  }
  return(mean_mat)
}

exponential_family_variance <- function(theta, family, weights = 1.0) {
  if (family == "gaussian") {
    var_mat = matrix(1, nrow(theta), ncol(theta)) * weights
  } else if (family == "binomial") {
    mean_mat = inv.logit.mat(theta)
    var_mat = mean_mat * (1 - mean_mat) * weights
  } else if (family == "poisson") {
    var_mat = exp(theta) * weights
  }
  return(var_mat)
}

exponential_family_log_likelihood <- function(mat, theta, family, weights = 1.0) {
  if (family == "gaussian") {
    return(sum(weights * 0.5 * (mat - theta)^2, na.rm = TRUE))
  } else if (family == "binomial") {
    return(sum(weights * (mat * theta - log(1 + exp(theta))), na.rm = TRUE))
  } else if (family == "poisson") {
    return(sum(weights * (mat * theta - exp(theta)), na.rm = TRUE))
  }
}
