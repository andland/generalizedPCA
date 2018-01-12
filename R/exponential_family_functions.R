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
  # .Call("inv_logit_mat", x, min, max, PACKAGE = "generalizedPCA")
  logisticPCA::inv.logit.mat(x, min, max)
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
# @export
log_like_Bernoulli <- function(x, theta, q) {
  logisticPCA::log_like_Bernoulli(x, theta, q)
}

check_family <- function(x, family) {
  distinct_vals = unique(c(x[!is.na(x)]))
  if (family %in% c("binomial", "multinomial") &&
      any(distinct_vals < 0 | distinct_vals > 1)) {
    stop(paste0(family, " family with data outside [0, 1]"))
  }
  if (any(family == "poisson" && distinct_vals < 0)) {
    stop("Negative data with poisson family")
  }
  if (family == "multinomial" && any(rowSums(x) > 1)) {
    stop(paste0(family, " family must have row sums <= 1"))
  }
  if (all(distinct_vals %in% c(0, 1))) {
    if (!(family %in% c("binomial", "multinomial"))) {
      message("All entries are binary. Are you sure you didn't mean binomial?")
    }
  } else if (all(distinct_vals >= 0 & distinct_vals %% 1 == 0)) {
    if (family != "poisson") {
      message("All entries are counts. Are you sure you didn't mean poisson?")
    }
  }
  if (any(distinct_vals %% 1 != 0) & family == "poisson") {
    message("Non-integer data with poisson data. Did you want family = guassian?")
  }
}

exp_fam_guess <- function(x) {
  distinct_vals = unique(c(x[!is.na(x)]))
  if (all(distinct_vals %in% c(0, 1))) {
    family = "binomial"
    message("Guessing data comes from a binomial distribution because it consists of all 0's and 1's")
  } else if (all(distinct_vals >= 0 & distinct_vals %% 1 == 0)) {
    family = "poisson"
    message("Guessing data comes from a poisson distribution because it consists of all non-negative integers")
  } else if (all(distinct_vals >=0 & distinct_vals <= 1) & all(rowSums(x, na.rm = TRUE) <= 1)) {
    family = "multinomial"
    message("Guessing data comes from a multinomial distribution because all the rows sum to less than 1")
  } else {
    family = "gaussian"
  }
  return(family)
}

saturated_natural_parameters <- function(x, family, M) {
  if (family == "gaussian") {
    eta = x
  } else if (family == "binomial") {
    eta = abs(M) * (2 * x - 1)
    non_binary = (x != 0 & x != 1 & !is.na(x))
    if (sum(non_binary) > 0) {
      logitvals = log(x) - log(1 - x)
      eta[non_binary] = logitvals[non_binary]
    }
  } else if (family == "poisson") {
    eta = log(x)
    eta[x==0] = -abs(M)
  } else if (family == "multinomial") {
    # TODO: check this! Correct if 0/1s.
    eta = abs(M) * (2 * x - 1)
    non_binary = (x != 0 & x != 1 & !is.na(x))
    if (sum(non_binary) > 0) {
      # TODO: give warning first time if any(rowSums(x) == 1)
      # TODO: doesn't deal with NAs,
      #   although typically a whole categorical variable would be NA

      # last_cat_prob = ifelse(rowSums(x) > 1 - as.numeric(inv.logit.mat(-as.numeric(M))),
      #                        as.numeric(inv.logit.mat(-as.numeric(M))), 1 - rowSums(x))
      # logitvals = sweep(log(x), 1, log(last_cat_prob), "-")

      # Below is in line with what is in tech report
      last_cat_prob = 1 - rowSums(x)
      last_cat_prob[last_cat_prob < exp(-as.numeric(M))] = exp(-as.numeric(M))
      logitvals = sweep(log(x), 1, log(last_cat_prob), "-")

      eta[non_binary] = logitvals[non_binary]
    }
  }
  eta[is.na(x)] <- 0
  return(eta)
}

# @export
exp_fam_mean <- function(theta, family) {
  if (family == "gaussian") {
    mean_mat = theta
  } else if (family == "binomial") {
    mean_mat = inv.logit.mat(theta)
  } else if (family == "poisson") {
    mean_mat = exp(theta)
  } else if (family == "multinomial") {
    exp_theta = exp(theta)
    mean_mat = sweep(exp_theta, 1, 1 + rowSums(exp_theta), "/")
  }
  return(mean_mat)
}

exp_fam_variance <- function(theta, family, weights = 1.0) {
  if (family == "gaussian") {
    var_mat = matrix(1, nrow(theta), ncol(theta)) * weights
  } else if (family %in% c("binomial", "multinomial")) {
    mean_mat = exp_fam_mean(theta, family)
    var_mat = mean_mat * (1 - mean_mat) * weights
  } else if (family == "poisson") {
    var_mat = exp(theta) * weights
  }
  return(var_mat)
}

# @export
exp_fam_log_like <- function(x, theta, family, weights = 1.0) {
  if (family == "gaussian") {
    return(-0.5 * sum(weights * (x - theta)^2, na.rm = TRUE))
  } else if (family == "binomial") {
    return(sum((weights * log(inv.logit.mat((2 * x - 1) * theta))), na.rm = TRUE))
    # the below does not work with x == 1 and theta > 750
    # return(sum(weights * (x * theta - log(1 + exp(theta))), na.rm = TRUE))
  } else if (family == "poisson") {
    return(sum(weights * (x * theta - exp(theta) - lfactorial(x)), na.rm = TRUE))
  } else if (family == "multinomial") {
    if (length(weights) > 1 && any(apply(weights, 1, stats::var) > 0)) {
      stop("weights should be the same for every variable withing each row")
    }
    return(sum(weights * (x * theta) - weights / ncol(x) *
                 outer(log(1 + rowSums(exp(theta))), rep(1, ncol(x))), na.rm = TRUE))
  }
}

#' Calculate exponential family deviance
#'
#' Calculates the deviance of data assuming it comes from an exponential family
#'
#' @param x a vector or matrix of data
#' @param theta natural parameters. Must be same dimensions as \code{x}
#' @param family exponential family distribution of data
#' @param weights weights of datapoints. Must be a scalar or the same dimensions as the data
#'
#' @export
exp_fam_deviance <- function(x, theta, family, weights = 1.0) {
  eta_sat_nat = saturated_natural_parameters(x, family, M = Inf)
  sat_loglike = exp_fam_log_like(x, eta_sat_nat, family, weights)

  -2 * (exp_fam_log_like(x, theta, family, weights) - sat_loglike)
}

# for solving for m
exp_fam_sat_ind_mat <- function(x, family) {
  if (family == "gaussian") {
    return(matrix(0, nrow(x), ncol(x)))
  } else if (family == "binomial") {
    # set 1's to 1, 0's to -1, and everything else to 0
    return((x == 1) - (x == 0))
  } else if (family == "poisson") {
    return((x == 0) * -1.0)
  } else if (family == "multinomial") {
    # first see if last cat == 0 and x between 0 and 1
    last_cat_prob_zero = rowSums(x) > (1 - 1e-3)
    q = outer(last_cat_prob_zero, rep(TRUE, ncol(x))) & x > 0 & x < 1
    # then add in exact 0's and 1's
    q = q + (x == 1) - (x == 0)
  }
}

# for exponential family harmoniums
exp_fam_sufficient_stat <- function(x, family) {
  if (family %in% c("gaussian", "binomial", "poisson")) {
    return(x)
  }
}

# for exponential family harmoniums
exp_fam_sample <- function(theta, family) {
  mean_mat = exp_fam_mean(theta = theta, family = family)

  if (family == "gaussian") {
    return(matrix(rnorm(prod(dim(theta)), mean = c(mean_mat), sd = 1), nrow(theta), ncol(theta)))
  } else if (family == "binomial") {
    # set 1's to 1, 0's to -1, and everything else to 0
    return(matrix(rbinom(prod(dim(theta)), size = 1, prob = c(mean_mat)), nrow(theta), ncol(theta)))
  } else if (family == "poisson") {
    return(matrix(rpois(prod(dim(theta)), lambda = c(mean_mat)), nrow(theta), ncol(theta)))
  }
}
