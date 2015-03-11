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

#' @export
check_family <- function(x, family) {
  distinct_vals = unique(c(x[!is.na(x)]))
  if (all(distinct_vals %in% c(0, 1))) {
    if (family != "binomial") {
      message("All entries are binary. Are you sure you didn't mean binomial?")
    }
  } else if (all(distinct_vals >= 0 & distinct_vals %% 1 == 0)) {
    if (family != "poisson") {
      message("All entries are counts. Are you sure you didn't mean poisson?")
    }
  }
}

#' @export
saturated_natural_parameters <- function(x, family, M) {
  if (family == "gaussian") {
    eta = x
  } else if (family == "binomial") {
    eta = abs(M) * (2 * x - 1)
  } else if (family == "poisson") {
    eta = log(x)
    eta[x==0] = -abs(M)
  }
  eta[is.na(x)] <- 0
  return(eta)
}

#' @export
exp_fam_mean <- function(theta, family) {
  if (family == "gaussian") {
    mean_mat = theta
  } else if (family == "binomial") {
    mean_mat = inv.logit.mat(theta)
  } else if (family == "poisson") {
    mean_mat = exp(theta)
  }
  return(mean_mat)
}

#' @export
exp_fam_variance <- function(theta, family, weights = 1.0) {
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

#' @export
exp_fam_log_like <- function(x, theta, family, weights = 1.0) {
  if (family == "gaussian") {
    # TODO: are these all loglikes and not negative?
    return(sum(weights * -0.5 * (x - theta)^2, na.rm = TRUE))
  } else if (family == "binomial") {
    return(sum(weights * (x * theta - log(1 + exp(theta))), na.rm = TRUE))
  } else if (family == "poisson") {
    return(sum(weights * (x * theta - exp(theta)), na.rm = TRUE))
  }
}
