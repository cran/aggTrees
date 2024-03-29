% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/balance-measures.R
\name{balance_measures}
\alias{balance_measures}
\title{Balance Measures}
\usage{
balance_measures(X, D)
}
\arguments{
\item{X}{Covariate matrix (no intercept).}

\item{D}{Treatment assignment vector.}
}
\value{
Prints LATEX code in the console.
}
\description{
Compute several balance measures to check whether the covariate distributions are balanced across
treatment arms.
}
\details{
For each covariate in \code{X}, \code{balance_measures} computes sample averages and standard deviations
for both treatment arms. Additionally, two balance measures are computed:
\describe{
  \item{\code{Norm. Diff.}}{Normalized differences, computed as the differences in the means of each covariate
  across treatment arms, normalized by the sum of the within-arm variances. They provide a measure of the
  discrepancy between locations of the covariate distributions across treatment arms.}
  \item{\code{Log S.D.}}{Log ratio of standard deviations are computed as the logarithm of the ratio of the
  within-arm standard deviations. They provide a measure of the
  discrepancy in the dispersion of the covariate distributions across treatment arms.}
  }

Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}.
}
\examples{
## Generate data.
set.seed(1986)

n <- 1000
k <- 3

X <- matrix(rnorm(n * k), ncol = k)
colnames(X) <- paste0("x", seq_len(k))
D <- rbinom(n, size = 1, prob = 0.5)
mu0 <- 0.5 * X[, 1]
mu1 <- 0.5 * X[, 1] + X[, 2]
y <- mu0 + D * (mu1 - mu0) + rnorm(n)

## Print table.
balance_measures(X, D)

}
\author{
Elena Dal Torrione, Riccardo Di Francesco
}
