% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dr_scores}
\alias{dr_scores}
\title{Doubly-Robust Scores}
\usage{
dr_scores(Y, D, X, k = 5)
}
\arguments{
\item{Y}{Outcome vector.}

\item{D}{Treatment assignment vector.}

\item{X}{Covariate matrix (no intercept).}

\item{k}{Number of folds.}
}
\value{
A vector of scores.
}
\description{
Constructs doubly-robust scores via K-fold cross-fitting.
}
\details{
Honest regression forests are used to estimate the propensity score and the conditional mean function of the outcome.
}
\author{
Riccardo Di Francesco
}
