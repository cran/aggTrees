## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(aggTrees)
library(grf)

## ----data-generation, eval = TRUE---------------------------------------------
## Generate data.
set.seed(1986)

n <- 500 # Small sample size due to compliance with CRAN notes.
k <- 3

X <- matrix(rnorm(n * k), ncol = k)
colnames(X) <- paste0("x", seq_len(k))
D <- rbinom(n, size = 1, prob = 0.5)
mu0 <- 0.5 * X[, 1]
mu1 <- 0.5 * X[, 1] + X[, 2]
Y <- mu0 + D * (mu1 - mu0) + rnorm(n)

## Sample split.
splits <- sample_split(length(Y), training_frac = 0.5)
training_idx <- splits$training_idx
honest_idx <- splits$honest_idx

Y_tr <- Y[training_idx]
D_tr <- D[training_idx]
X_tr <- X[training_idx, ]

Y_hon <- Y[honest_idx]
D_hon <- D[honest_idx]
X_hon <- X[honest_idx, ]

## ----estimate-cates, eval = TRUE----------------------------------------------
## Estimate the CATEs. Use only training sample.
forest <- causal_forest(X_tr, Y_tr, D_tr) 

cates_tr <- predict(forest, X_tr)$predictions
cates_hon <- predict(forest, X_hon)$predictions

## ----construct-sequence, eval = TRUE------------------------------------------
## Construct the sequence. Use doubly-robust scores (default option).
groupings <- build_aggtree(Y_tr, D_tr, X_tr, # Training sample. 
                           Y_hon, D_hon, X_hon, # Honest sample.
                           cates_tr = cates_tr, cates_hon = cates_hon) # Predicted CATEs.

## Print.
print(groupings)

## Plot.
plot(groupings) # Try also setting 'sequence = TRUE'.

## ----inference, eval = TRUE---------------------------------------------------
## Inference with 4 groups.
results <- inference_aggtree(groupings, n_groups = 4)

## LATEX.
print(results, table = "diff")

print(results, table = "avg_char")

