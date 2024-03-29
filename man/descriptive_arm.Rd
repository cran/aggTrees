% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/balance-measures.R
\name{descriptive_arm}
\alias{descriptive_arm}
\title{Descriptive Statistics by Treatment Arm}
\usage{
descriptive_arm(X, D)
}
\arguments{
\item{X}{Covariate matrix (no intercept).}

\item{D}{Treatment assignment vector.}
}
\value{
4xp array, storing the desired statistics.
}
\description{
Computes sample averages and standard deviations of the covariates across treatment arms.
}
\details{
Sample means and standard deviations across treatment arms are a first, useful insight to assess covariate balance.
}
\keyword{internal}
