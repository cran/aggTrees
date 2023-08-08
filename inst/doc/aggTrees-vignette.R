## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(aggTrees)

## ----data-generation, eval = FALSE--------------------------------------------
#  ## Generate data.
#  set.seed(1986)
#  
#  n <- 1000
#  k <- 3
#  
#  X <- matrix(rnorm(n * k), ncol = k)
#  colnames(X) <- paste0("x", seq_len(k))
#  D <- rbinom(n, size = 1, prob = 0.5)
#  mu0 <- 0.5 * X[, 1]
#  mu1 <- 0.5 * X[, 1] + X[, 2]
#  y <- mu0 + D * (mu1 - mu0) + rnorm(n)

## ----estimate-cates, eval = FALSE---------------------------------------------
#  ## Sample split.
#  splits <- sample_split(length(y), training_frac = 0.5)
#  training_idx <- splits$training_idx
#  honest_idx <- splits$honest_idx
#  
#  y_tr <- y[training_idx]
#  D_tr <- D[training_idx]
#  X_tr <- X[training_idx, ]
#  
#  y_hon <- y[honest_idx]
#  D_hon <- D[honest_idx]
#  X_hon <- X[honest_idx, ]
#  
#  ## Estimate the CATEs. Use only training sample.
#  library(grf)
#  forest <- causal_forest(X_tr, y_tr, D_tr)
#  cates <- predict(forest, X)$predictions

## ----construct-sequence, eval = FALSE-----------------------------------------
#  ## Construct the sequence. Use doubly-robust scores.
#  groupings <- build_aggtree(y, D, X, method = "aipw",
#                             cates = cates, is_honest = 1:length(y) %in% honest_idx)
#  
#  ## Print.
#  print(groupings)
#  
#  ## Plot.
#  plot(groupings) # Try also setting 'sequence = TRUE'.

## ----inference, eval = FALSE--------------------------------------------------
#  ## Inference with 4 groups.
#  results <- inference_aggtree(groupings, n_groups = 4)
#  
#  ## LATEX.
#  print(results, table = "diff")
#  print(results, table = "avg_char")

