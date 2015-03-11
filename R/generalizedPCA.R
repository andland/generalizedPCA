#' @title Logistic Principal Component Analysis
#'
#' @description
#' Dimension reduction for exponential family data by extending Pearson's
#' PCA formulation
#'
#' @param x matrix of either binary, count, or continuous data
#' @param k number of principal components to return
#' @param M value to approximate the saturated model
#' @param family exponential family distribution of data
#' @param weights a matrix of the same size as the \code{x} with data weights
#' @param quiet logical; whether the calculation should give feedback
#' @param use_irlba logical; if \code{TRUE}, the function uses the irlba package
#'   to more quickly calculate the eigen-decomposition
#' @param max_iters number of maximum iterations
#' @param conv_criteria convergence criteria. The difference between average deviance
#'   in successive iterations
#' @param random_start logical; whether to randomly inititalize the parameters. If \code{FALSE},
#'   function will use an eigen-decomposition as starting value
#' @param start_U starting value for the orthogonal matrix
#' @param start_mu starting value for mu. Only used if \code{main_effects = TRUE}
#' @param main_effects logical; whether to include main effects in the model
#'
#' @return An S3 object of class \code{gpca} which is a list with the
#' following components:
#' \item{mu}{the main effects}
#' \item{U}{a \code{k}-dimentional orthonormal matrix with the loadings}
#' \item{PCs}{the princial components}
#' \item{M}{the parameter inputed}
#' \item{family}{the exponential family used}
#' \item{iters}{number of iterations required for convergence}
#' \item{loss_trace}{the trace of the average negative log likelihood of the algorithm.
#'    Should be non-increasing}
#' \item{prop_deviance_expl}{the proportion of deviance explained by this model.
#'    If \code{main_effects = TRUE}, the null model is just the main effects, otherwise
#'    the null model estimates 0 for all natural parameters.}
#'
#' @examples
#' @export
logisticPCA <- function(x, k = 2, M = 4, family = c("gaussian", "binomial", "poisson"),
                        majorizer = c("row", "all"), weights,
                        quiet = TRUE, use_irlba = FALSE,
                        max_iters = 1000, conv_criteria = 1e-5, random_start = FALSE,
                        start_U, start_mu, main_effects = TRUE) {

  object <- list(mu = mu,
                 U = U,
                 PCs = scale(eta, center = mu, scale = FALSE) %*% U,
                 M = M,
                 family = family,
                 iters = m,
                 loss_trace = loss_trace[1:(m + 1)],
                 prop_deviance_expl = 1 - loglike / null_loglike)
  class(object) <- "gpca"
  object
}


#' @title Predict generalized PCA scores or reconstruction on new data
#'
#' @param object generalized PCA object
#' @param newdata matrix of the same exponential family as in \code{object}.
#'  If missing, will use the data that \code{object} was fit on
#' @param type the type of fitting required. \code{type = "PCs"} gives the PC scores,
#'  \code{type = "link"} gives matrix on the natural parameter scale and
#'  \code{type = "response"} gives matrix on the response scale
#' @param ... Additional arguments
#' @examples
#' @export
predict.lpca <- function(object, newdata, type = c("PCs", "link", "response"), ...) {
  type = match.arg(type)


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
#' @export
fitted.lpca <- function(object, type = c("link", "response"), ...) {
  type = match.arg(type)
  n = nrow(object$PCs)

  theta = outer(rep(1, n), object$mu) + tcrossprod(object$PCs, object$U)

  if (type == "link") {
    return(theta)
  } else if (type == "response") {
    return(exponential_family_mean(theta))
  }
}

#' @title Plot generalized PCA
#'
#' @description
#' Plots the results of a generalized PCA
#'
#' @param object generalized PCA object
#' @param type the type of plot \code{type = "trace"} plots the algorithms progress by
#' iteration, \code{type = "loadings"} plots the first 2 principal component
#' loadings, \code{type = "scores"} plots the loadings first 2 principal component scores
#' @param ... Additional arguments
#' @examples
#' @export
plot.lpca <- function(object, type = c("trace", "loadings", "scores"), ...) {
  library("ggplot2")
  type = match.arg(type)

  if (type == "trace") {
    df = data.frame(Iteration = 0:object$iters,
                    NegativeLogLikelihood = object$loss_trace)
    p <- ggplot2::ggplot(df, aes(Iteration, NegativeLogLikelihood)) +
      geom_line()
  } else if (type == "loadings") {
    df = data.frame(object$U)
    colnames(df) <- paste0("PC", 1:ncol(df))
    if (ncol(df) == 1) {
      p <- ggplot2::qplot(PC1, 0, data = df, ylab = NULL)
    } else {
      p <- ggplot2::ggplot(df, aes(PC1, PC2)) + geom_point()
    }
  } else if (type == "scores") {
    df = data.frame(object$PCs)
    colnames(df) <- paste0("PC", 1:ncol(df))
    if (ncol(df) == 1) {
      p <- ggplot2::qplot(PC1, 0, data = df, ylab = NULL)
    } else {
      p <- ggplot2::ggplot(df, aes(PC1, PC2)) + geom_point()
    }
  }

  return(p)
}

#' @title Print generalized PCA object
#'
#' @param x generalized PCA object
#' @param ... Additional arguments
#'
#' @export
print.lpca <- function(x, ...) {
  cat(nrow(x$PCs), "rows and ")
  cat(nrow(x$U), "columns\n")
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
#' @param folds if \code{folds} is a scalar, then it is the number of folds. If
#'  it is a vector, it should be the same length as the number of rows in \code{x}
#' @param quiet logical; whether the function should display progress
#' @param ... Additional arguments passed to \code{generalizedPCA}
#'
#' @return A matrix of the CV log likelihood with \code{k} in rows and
#'  \code{M} in columns
#'
#' @examples
#' @export
cv.gpca <- function(x, ks, Ms = seq(2, 10, by = 2), family, folds = 5, quiet = TRUE, ...) {
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
        gpca = generalizedPCA(x[c != cv, ], k = k, M = M, family = family, ...)
        pred_theta = predict(gpca, newdat = x[c == cv, ], type = "link")
        log_likes[k == ks, M == Ms] = log_likes[k == ks, M == Ms] +
          exponential_family_log_likelihood(x = x[c == cv, ], theta = pred_theta, family = family)
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
#' @param object a \code{cv.gpca} object
#' @param ... Additional arguments
#'
#' @examples
#' @export
plot.cv.gpca <- function(object, ...) {
  library(ggplot2)
  library(reshape2)
  df = melt(-object, value.name = "NegLogLikelihood")
  if (ncol(object) == 1) {
    df$M = factor(df$M)
    p <- ggplot(df, aes(k, NegLogLikelihood, colour = M)) +
      geom_line()
  } else {
    df$k = factor(df$k)
    p <- ggplot(df, aes(M, NegLogLikelihood, colour = k)) +
      geom_line()
  }
  return(p)
}
