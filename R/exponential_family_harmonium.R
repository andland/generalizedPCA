#' Exponential Family Harmoniums
#'
#' @description Welling et al. (2005)'s Exponential Family Harmoniums
#'
#' @param x matrix of either binary, count, or continuous data
#' @param k dimension (number of hidden units)
#' @param family exponential family distribution of data
#' @param family_hidden exponential family distribution of hidden units
#' @param cd_iters number of iterations for contrastive divergence (CD) at each iteration.
#'   It should be a vector of two integers, which is the range of CD interations that the
#'   algorithm will perform from the beginning until the end, linearly interpolated
#' @param max_iters maximum number of iterations
#' @param learning_rate learning rate used for gradient descent
#' @param quiet logical; whether the calculation should give feedback
#' @param random_start whether to randomly initialize \code{W}
#' @param start_W initial value for \code{W}
#' @param mu specific value for \code{mu}, the mean (bias) vector of \code{x}
#' @param main_effects logical; whether to include main effects (bias terms) in the model
#'
#' @return An S3 object of class \code{efh} which is a list with the
#' following components:
#' \item{mu}{the main effects (bias terms) for dimensionality reduction}
#' \item{hidden_bias}{the bias for the hidden units (currently hard coded to 0)}
#' \item{W}{the \code{d}x\code{k}-dimentional matrix with the loadings}
#' \item{family}{the exponential family of the data}
#' \item{family_hidden}{the exponential family of the hidden units}
#' \item{iters}{number of iterations required for convergence}
#' \item{loss_trace}{the trace of the average deviance of the algorithm.
#'    Should be non-increasing}
#' \item{prop_deviance_expl}{the proportion of deviance explained by this model.
#'    If \code{main_effects = TRUE}, the null model is just the main effects, otherwise
#'    the null model estimates 0 for all natural parameters.}
#'
#' @export
#'
#' @references
#' Welling, Max, Michal Rosen-Zvi, and Geoffrey E. Hinton. "Exponential family
#' harmoniums with an application to information retrieval." Advances in neural
#' information processing systems. 2005.
#'
#' @examples
#' rows = 100
#' cols = 10
#' set.seed(1)
#' mat_np = outer(rnorm(rows), rnorm(cols))
#'
#' # generate a count matrix and binary response
#' mat = matrix(rpois(rows * cols, c(exp(mat_np))), rows, cols)
#' mat[1, 1] <- NA
#'
#' modp = exponential_family_harmonium(mat, k = 2, family = "poisson", quiet = FALSE, learning_rate = 0.001, rms_prop = F, max_iters = 250)
#' gmf = generalizedMF(mat, k = 1, family = "poisson", quiet = FALSE)
#'
exponential_family_harmonium <- function(x, k = 2,
                          family = c("gaussian", "binomial", "poisson"),
                          family_hidden = c("gaussian", "binomial", "poisson"),
                          cd_iters = 1, learning_rate = 0.001, max_iters = 100, rms_prop = FALSE,
                          quiet = TRUE, random_start = TRUE, start_W,  mu, main_effects = TRUE) {
  family = match.arg(family)
  family_hidden = match.arg(family_hidden)
  check_family(x, family)
  # x = mat; k = 1; family = "poisson"; family_hidden = "gaussian"; cd_iters = 1; learning_rate = 0.001; max_iters = 100; quiet = FALSE;

  stopifnot(cd_iters > 0.5, length(cd_iters) %in% c(1, 2, max_iters))
  if (length(cd_iters) == 1) {
    cd_iters = rep(round(cd_iters), max_iters)
  } else if (length(cd_iters) == 2) {
    cd_iters = round(cd_iters)
    cd_iters = pmax(min(cd_iters), pmin(max(cd_iters), round(seq(cd_iters[1] - 0.5, cd_iters[2] + 0.5, length.out = max_iters))))
  } else {
    cd_iters = round(cd_iters)
  }

  x = as.matrix(x)
  missing_mat = is.na(x)
  sum_weights = sum(!missing_mat)

  # missing values are ignored, which is equivalent to setting them to 0
  x_imputed = x
  x_imputed[missing_mat] <- 0
  # scaled to take the avg with different number of missing values
  x_imputed_scaled = scale(x_imputed, FALSE, colSums(!missing_mat))

  n = nrow(x)
  d = ncol(x)
  ones = rep(1, n)

  # calculate the null log likelihood for % deviance explained and normalization
  if (main_effects) {
    weighted_col_means = colMeans(x, na.rm = TRUE)
    null_theta = as.numeric(saturated_natural_parameters(matrix(weighted_col_means, 1), family, M = Inf))
  } else {
    null_theta = rep(0, d)
  }
  null_deviance = exp_fam_deviance(x, outer(ones, null_theta), family) / sum_weights

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

  # hard code to 0 for now
  hidden_bias = rep(0, k)

  if (!missing(start_W)) {
    stopifnot(dim(start_W) == c(d, k))
    W = as.matrix(start_W)
  } else if (random_start) {
    W = matrix(rnorm(d * k, 0, 0.1), d, k)
  } else {
    udv = svd(scale(saturated_natural_parameters(x_imputed, family, 4), TRUE, FALSE))
    W = udv$v[, 1:k, drop = FALSE]
  }

  loss_trace = numeric(max_iters + 1)
  theta = outer(ones, mu) + exp_fam_mean(outer(ones, hidden_bias) + x_imputed %*% W, family_hidden) %*% t(W)
  loss_trace[1] = exp_fam_deviance(x, theta, family) / sum_weights

  ptm <- proc.time()

  if (!quiet) {
    cat(0, "  ", loss_trace[1], "")
    cat("0 hours elapsed\n")
  }

  # RMSProp
  if (rms_prop) {
    W_grad_sq = matrix(1, nrow(W), ncol(W))
  } else {
    W_grad_sq = matrix(1, nrow(W), ncol(W))
  }

  for (ii in seq_len(max_iters)) {

    while (TRUE) {
      hidden_mean_0 = exp_fam_mean(outer(ones, hidden_bias) + x_imputed %*% W, family_hidden)
      visible_hidden_0 = t(x_imputed_scaled) %*% hidden_mean_0

      visible_cd = contrastive_divergence(x = x_imputed, W = W, mu = mu, hidden_bias = hidden_bias,
                                          family = family, family_hidden = family_hidden,
                                          num_iter = cd_iters[ii])
      hidden_mean_cd = exp_fam_mean(outer(ones, hidden_bias) + visible_cd %*% W, family_hidden)
      visible_hidden_cd = t(visible_cd) %*% hidden_mean_cd / n

      # the generated data can explode with gaussian hidden and poisson visible
      if (any(is.na(visible_hidden_cd))) {
        if (ii > 1) {
          W_grad_sq = W_grad_sq * 4
          W = W_lag + learning_rate * W_grad / sqrt(W_grad_sq)
        } else {
          W = W / 2
        }
        next
      }

      W_grad = visible_hidden_0 - visible_hidden_cd
      W = W + learning_rate * W_grad / sqrt(W_grad_sq)
      if (rms_prop) {
        if (ii == 1) {
          W_grad_sq = W_grad^2
        } else {
          W_grad_sq = 0.1 * W_grad^2 + 0.9 * W_grad_sq
        }
      }

      # Calc Deviance
      theta = outer(ones, mu) + exp_fam_mean(outer(ones, hidden_bias) + x_imputed %*% W, family_hidden) %*% t(W)
      loss_trace[ii + 1] <- exp_fam_deviance(x, theta, family) / sum_weights

      if (loss_trace[ii + 1] <= loss_trace[1]) {
        break
      } else {
        W_grad_sq = W_grad_sq * 4
      }
    }

    if (!quiet) {
      time_elapsed = as.numeric(proc.time() - ptm)[3]
      tot_time = max_iters / ii * time_elapsed
      time_remain = tot_time - time_elapsed
      cat(ii, "  ", loss_trace[ii + 1], "")
      cat(round(time_elapsed / 3600, 1), "hours elapsed. Max", round(time_remain / 3600, 1), "hours remain.\n")
    }

    if (ii > max_iters) {
      break
    } else {
      W_lag = W
    }
  }

  object <- list(
    mu = mu,
    W = W,
    hidden_bias = hidden_bias,
    family = family,
    family_hidden = family_hidden,
    iters = ii,
    cd_iters = cd_iters[1:ii],
    loss_trace = loss_trace[2:(ii + 1)],
    prop_deviance_expl = 1 - loss_trace[ii + 1] / null_deviance,
    W_grad_sq = W_grad_sq
  )
  class(object) <- "efh"
  object
}

contrastive_divergence <- function(x, W, mu, hidden_bias, family, family_hidden, num_iter) {
  visible = x
  ones = rep(1, nrow(x))
  for (cd_ii in seq_len(num_iter)) {
    hidden = exp_fam_sample(outer(ones, hidden_bias) + visible %*% W, family_hidden)
    visible = exp_fam_sample(outer(ones, mu) + hidden %*% t(W), family)
  }
  return(visible)
}

simulate_efh <- function(model, x, num_iter) {
  visible = x
  ones = rep(1, nrow(x))
  for (cd_ii in seq_len(num_iter)) {
    hidden = exp_fam_sample(outer(ones, model$hidden_bias) + visible %*% model$W, model$family_hidden)
    visible = exp_fam_sample(outer(ones, model$mu) + hidden %*% t(model$W), model$family)
  }
  return(visible)
}

#' contrastive_divergence(mat, )
