# Generalized PCA

`generalizedPCA` is an R package which extends principal component analysis to other types of data, much like [generalized linear models](http://en.wikipedia.org/wiki/Generalized_linear_model) extends linear regression. The package [logisticPCA](https://github.com/andland/logisticPCA) contains the extension to binary data, among other methods, and this package intends to generalize it to all exponential family distributions. Please note that it is still in the very early stages of development and the conventions will possibly change in the future.

## Installation

To install R, visit [r-project.org/](http://www.r-project.org/).

To install the package, first install `devtools` from CRAN. Then run the following commands.
```R
# install.packages("devtools")
library("devtools")
install_github("andland/generalizedPCA")
```

## Use
The main function is `generalizedPCA()`. Like in generalized linear models, you must specify the distribution of your data. `generalizedPCA()` currently supports `"gaussian"`, `"binomial"`, `"poisson"`, or `"multinomial"` data. Unlike standard PCA, it can incorporate weights and missing data.  If your data are proportions, you can use `family = "binomial"` with `weights` being a matrix of the number of opportunities. If your data is a multinomial variable with `d` levels, the input matrix should have `d - 1` columns.

The function returns `mu`, the variable main effects vector of length `d`, and `U`, the `d x k` loadings matrix.

## Details
`generalizedPCA()` estimates the natural parameters of an exponential family distribution in a lower dimensional space. This is done by projecting the natural parameters from the saturated model. A rank-`k` projection matrix, or equivalently a `d x k` orthogonal matrix, is solved for to minimize the deviance. 

For some distributions, the natural parameters from the saturated model are either negative or positive infinity, and an additional tuning parameter `M` is needed to approximate them. This occurs when `family = "binomial"` and your data include `0`'s or `1`'s or when `family = "poisson"` and your data include `0`'s. I usually use `cv.gpca()` to select `M`. Typical values are in the range `3` to `10`.

## Methods
The generalizedPCA class, `gpca`, has several methods to make data analysis easier.

* `print()`: Prints a summary of the fitted model.
* `fitted()`: Fits the low dimensional matrix of either natural parameters or response.
* `predict()`: Predicts the PCs on new data. Can also predict the low dimensional matrix of natural parameters or response on new data.
* `plot()`: Plots either the deviance trace by the number of iterations, the first two PC loadings, or the first two PC scores using the package `ggplot2`.

In addition, there are functions for performing cross-validation.

* `cv.gpca()`: Run cross-validation over the rows of the matrix to assess the fit of `M` and/or `k`.
* `plot.cv()`: Plots the results of the `cv.gpca()` function using the package `ggplot2`.
