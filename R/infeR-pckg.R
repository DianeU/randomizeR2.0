#' Randomization Based Inference in Clinical Trials
#'
#' \tabular{ll}{
#' Package: \tab infeR \cr
#' Type: \tab Package\cr
#' Version: \tab 1.0\cr
#' Date: \tab 2020-06-03\cr
#' License: \tab GPL (>= 3)\cr
#' LazyLoad: \tab yes \cr
#' }
#'
#' This package implements functions to conduct randomization tests in clinical trials.
#' Randomization tests are non-parametric tests that only rely on the randomization distribution
#' of the test statistic. The randomization distribution is induced by the randomization procedure
#' of the clinical trial. The infeR package relies heavily on the randomizeR R package. It
#' implements different test statistics found in the literature and provides a convenient interface
#' to conduct randomization tests, and simulation studies that use randomization tests.
#'
#'
#' @docType package
#' @name infeR-package
#' @aliases infeR
#' @title Randomization based inference for clinical trials
#' @author Diane Uschner \email{diane.uschner@@gmail.com}
#' @references Uschner, D: Randomization based inference in the presence of selection bias. Submitted to Statistics in Medicine, 2020.
#' @import randomizeR
NULL
