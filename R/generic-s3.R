#' Plot Method for aggTrees Objects
#'
#' Plots an \code{aggTrees} object.
#'
#' @param x An \code{aggTrees} object.
#' @param leaves Number of leaves of the desired tree. This can be used to plot subtrees.
#' @param sequence If \code{TRUE}, the whole sequence of optimal groupings is displayed in a short animation.
#' @param ... Further arguments from \code{\link[rpart.plot]{prp}}.
#'
#' @return
#' Plots an \code{aggTrees} object.
#'
#' @examples
#' \donttest{## Generate data.
#' set.seed(1986)
#'
#' n <- 1000
#' k <- 3
#'
#' X <- matrix(rnorm(n * k), ncol = k)
#' colnames(X) <- paste0("x", seq_len(k))
#' D <- rbinom(n, size = 1, prob = 0.5)
#' mu0 <- 0.5 * X[, 1]
#' mu1 <- 0.5 * X[, 1] + X[, 2]
#' Y <- mu0 + D * (mu1 - mu0) + rnorm(n)
#'
#' ## Training-honest sample split.
#' honest_frac <- 0.5
#' splits <- sample_split(length(Y), training_frac = (1 - honest_frac))
#' training_idx <- splits$training_idx
#' honest_idx <- splits$honest_idx
#'
#' Y_tr <- Y[training_idx]
#' D_tr <- D[training_idx]
#' X_tr <- X[training_idx, ]
#'
#' Y_hon <- Y[honest_idx]
#' D_hon <- D[honest_idx]
#' X_hon <- X[honest_idx, ]
#'
#' ## Construct sequence of groupings. CATEs estimated internally.
#' groupings <- build_aggtree(Y_tr, D_tr, X_tr,
#'                            Y_hon, D_hon, X_hon)
#'
#' ## Plot.
#' plot(groupings)
#' plot(groupings, leaves = 3)
#' plot(groupings, sequence = TRUE)}
#'
#' @details
#' Nodes are colored using a diverging palette. Nodes with predictions smaller than the ATE (i.e., the root
#' prediction) are colored in blue shades, and nodes with predictions larger than the ATE are colored in red
#' shades. Moreover, predictions that are more distant in absolute value from the ATE get darker shades.
#' This way, we have an immediate understanding of the groups with extreme GATEs.
#'
#' @import rpart.plot
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @seealso
#' \code{\link{build_aggtree}}, \code{\link{inference_aggtree}}
#'
#' @export
plot.aggTrees <- function(x, leaves = get_leaves(x$tree), sequence = FALSE, ...) {
  ## Handling inputs and checks.
  if (!(sequence %in% c(TRUE, FALSE))) stop("Invalid 'sequence'. This must be either TRUE or FALSE.", call. = FALSE)

  tree <- x$tree

  ## Plotting.
  if (sequence) {
    ## Node labels.
    prefix <- c("ATE = ", rep("GATE = ", times = (length(tree$frame$n) - 1)))
    sizes <- tree$frame$n
    percentages <- round(tree$frame$n / tree$frame$n[1] * 100, 0)
    suffix <- paste("\n", "n = ", sizes, "\n (", percentages, "%)", sep = "")

    ## Threshold values of cost-complexity parameter and node numbers.
    nodes <- as.numeric(rownames(tree$frame))
    alpha_values <- rev(tree$cptable[, "CP"])

    ## Settings for node colors.
    low_colors <- c("#002F70", "#6889D0", "#F6F6F6")
    high_colors <- c("#F6F6F6", "#CB6F70", "#5F1415")
    color_fun_low <- grDevices::colorRampPalette(low_colors)
    color_fun_high <- grDevices::colorRampPalette(high_colors)

    low_palette <- color_fun_low(sum(tree$frame$yval <= tree$frame$yval[1]))
    high_palette <- color_fun_high(sum(tree$frame$yval >= tree$frame$yval[1]))[-1]

    ## Plot sequence.
    for (alpha in alpha_values[1:length(alpha_values[-1])]) {
      ## Graying out pruned leaves.
      temp_nodes <- as.numeric(rownames(rpart::prune.rpart(tree, alpha)$frame))
      colors <- ifelse(nodes %in% temp_nodes, 1, "white")

      ## Get rid of pruned leaves, branches and split labels by using with colors. Also, order colors as in Details.
      box_colors <- c(low_palette, high_palette)
      ate <- tree$frame$yval[1]
      to_order_bottom <- sort(tree$frame$yval[-1][tree$frame$yval[-1] < ate])
      to_order_top <- sort(tree$frame$yval[-1][tree$frame$yval[-1] > ate])
      to_order_box_colors <- c(to_order_bottom, ate, to_order_top)
      box_colors <- box_colors[match(tree$frame$yval, to_order_box_colors)]
      box_colors[colors == "white"] <- "white"

      temp_parents <- temp_nodes %/% 2
      split_colors <- ifelse(nodes %in% temp_parents, "black", "white")

      ## Plot.
      grDevices::dev.hold()

      rpart.plot::prp(tree,
                      type = 2,
                      extra = 0,
                      under = FALSE,
                      fallen.leaves = TRUE,
                      round = 0,
                      leaf.round = 0,
                      prefix = prefix,
                      suffix = suffix,
                      box.col = box_colors,
                      branch = 0.3,
                      branch.col = colors,
                      split.col = split_colors,
                      col = colors,
                      ...)

      grDevices::dev.flush()
      Sys.sleep(1)
    }
  } else {
    ## Get subtree.
    subtree <- subtree(tree, leaves)

    ## Node labels.
    prefix <- c("ATE = ", rep("GATE = ", times = (length(subtree$frame$n) - 1)))
    sizes <- subtree$frame$n
    percentages <- round(subtree$frame$n / subtree$frame$n[1] * 100, 0)
    suffix <- paste("\n", "n = ", sizes, "\n (", percentages, "%)", sep = "")

    ## Settings for node colors.
    low_colors <- c("#002F70", "#6889D0", "#F6F6F6")
    high_colors <- c("#F6F6F6", "#CB6F70", "#5F1415")
    color_fun_low <- grDevices::colorRampPalette(low_colors)
    color_fun_high <- grDevices::colorRampPalette(high_colors)

    low_palette <- color_fun_low(sum(subtree$frame$yval <= subtree$frame$yval[1]))
    high_palette <- color_fun_high(sum(subtree$frame$yval >= subtree$frame$yval[1]))[-1]
    box_colors <- c(low_palette, high_palette)

    ate <- subtree$frame$yval[1]
    to_order_bottom <- sort(subtree$frame$yval[-1][subtree$frame$yval[-1] < ate])
    to_order_top <- sort(subtree$frame$yval[-1][subtree$frame$yval[-1] > ate])
    to_order_box_colors <- c(to_order_bottom, ate, to_order_top)
    box_colors <- box_colors[match(subtree$frame$yval, to_order_box_colors)]

    ## Plot.
    rpart.plot::prp(subtree,
                    type = 2,
                    extra = 0,
                    under = FALSE,
                    fallen.leaves = TRUE,
                    round = 0,
                    leaf.round = 0,
                    prefix = prefix,
                    suffix = suffix,
                    box.col = box_colors,
                    branch = 0.3,
                    ...)
  }
}


#' Summary Method for aggTrees Objects
#'
#' Summarizes an \code{aggTrees} object.
#'
#' @param object \code{aggTrees} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' Prints the summary of an \code{aggTrees} object.
#'
#' @seealso
#' \code{\link{build_aggtree}}, \code{\link{inference_aggtree}}
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @author Riccardo Di Francesco
#'
#' @export
summary.aggTrees <- function(object, ...) {
  if (is.null(object$honest_sample)) cat("Honest estimates:", FALSE, "\n") else cat("Honest estimates:", TRUE, "\n")
  summary(object$tree)
}


#' Print Method for aggTrees Objects
#'
#' Prints an \code{aggTrees} object.
#'
#' @param x \code{aggTrees} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' Prints an \code{aggTrees} object.
#'
#' @seealso
#' \code{\link{build_aggtree}}, \code{\link{inference_aggtree}}
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @author Riccardo Di Francesco
#'
#' @export
print.aggTrees <- function(x, ...) {
  if (is.null(x$honest_sample)) cat("Honest estimates:", FALSE, "\n") else cat("Honest estimates:", TRUE, "\n")
  print(x$tree)
}


#' Print Method for aggTrees.inference Objects
#'
#' Prints an \code{aggTrees.inference} object.
#'
#' @param x \code{aggTrees.inference} object.
#' @param table Either \code{"avg_char"} or \code{"diff"}, controls which table must be produced.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' Prints LATEX code.
#'
#' @examples
#' \donttest{## Generate data.
#' set.seed(1986)
#'
#' n <- 1000
#' k <- 3
#'
#' X <- matrix(rnorm(n * k), ncol = k)
#' colnames(X) <- paste0("x", seq_len(k))
#' D <- rbinom(n, size = 1, prob = 0.5)
#' mu0 <- 0.5 * X[, 1]
#' mu1 <- 0.5 * X[, 1] + X[, 2]
#' Y <- mu0 + D * (mu1 - mu0) + rnorm(n)
#'
#' ## Training-honest sample split.
#' honest_frac <- 0.5
#' splits <- sample_split(length(Y), training_frac = (1 - honest_frac))
#' training_idx <- splits$training_idx
#' honest_idx <- splits$honest_idx
#'
#' Y_tr <- Y[training_idx]
#' D_tr <- D[training_idx]
#' X_tr <- X[training_idx, ]
#'
#' Y_hon <- Y[honest_idx]
#' D_hon <- D[honest_idx]
#' X_hon <- X[honest_idx, ]
#'
#' ## Construct sequence of groupings. CATEs estimated internally.
#' groupings <- build_aggtree(Y_tr, D_tr, X_tr,
#'                            Y_hon, D_hon, X_hon)
#'
#' ## Analyze results with 4 groups.
#' results <- inference_aggtree(groupings, n_groups = 4)
#'
#' ## Print results.
#' print(results, table = "diff")
#' print(results, table = "avg_char")}
#'
#' @details
#' A description of each table is provided in its caption.\cr
#'
#' Some covariates may feature zero variation in some leaf. This generally happens to dummy variables used to split some
#' nodes. In this case, when \code{table == "avg_char"} a warning message is produced displaying the names of the covariates
#' with zero variation in one or more leaves. The user should correct the table by removing the associated standard errors.\cr
#'
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox},
#' \code{multirow}.
#'
#' @seealso
#' \code{\link{build_aggtree}}, \code{\link{inference_aggtree}}
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @author Riccardo Di Francesco
#'
#' @export
print.aggTrees.inference <- function(x, table = "avg_char", ...) {
  ## Check.
  if (!(table %in% c("avg_char", "diff"))) stop("Invalid 'table'. This must be either 'avg_char' or 'diff'.", call. = FALSE)

  ## Select appropriate sample (adaptive/honest) according to the output of build_aggtree.
  if (is.null(x$aggTree$honest_sample)) {
    X <- x$aggTree$training_sample[, -c(1:3)]
  } else {
    X <- x$aggTree$honest_sample[, -c(1:3)]
  }

  ## Extract information.
  leaves <- leaf_membership(x$groups, X)
  parms <- lapply(x$avg_characteristics, function(x) {stats::coef(summary(x))[, c("Estimate", "Std. Error")]})

  if (x$aggTree$method == "raw") {
    gates_idx <- which(sapply(names(x$model$coefficients), function(x) grepl(":D", x)))
  } else if (x$aggTree$method == "aipw") {
    gates_idx <- which(sapply(names(x$model$coefficients), function(x) grepl("leaf", x)))
  }

  gates_point <- format(round(x$model$coefficients[gates_idx], 3), nsmall = 3)
  gates_sd <- format(round(x$model$std.error[gates_idx], 3), nsmall = 3)

  gates_ci_lower <- format(round(as.numeric(gates_point) - 1.96 * as.numeric(gates_sd), 3), nsmall = 3)
  gates_ci_upper <- format(round(as.numeric(gates_point) + 1.96 * as.numeric(gates_sd), 3), nsmall = 3)

  if (length(x$boot_ci) != 0) gates_ci_lower_boot <- format(round(x$boot_ci$lower, 3), nsmall = 3) else gates_ci_lower_boot <- NA
  if (length(x$boot_ci) != 0) gates_ci_upper_boot <- format(round(x$boot_ci$upper, 3), nsmall = 3) else gates_ci_upper_boot <- NA

  ## Write table.
  if (table == "avg_char") {
    table_names <- rename_latex(colnames(X))

    cat("\\begingroup
  \\setlength{\\tabcolsep}{8pt}
  \\renewcommand{\\arraystretch}{1.1}
  \\begin{table}[b!]
    \\centering
    \\begin{adjustbox}{width = 1\\textwidth}
    \\begin{tabular}{@{\\extracolsep{5pt}}l ", rep("c ", (length(unique(leaves)) * 2)), "}
      \\\\[-1.8ex]\\hline
      \\hline \\\\[-1.8ex]
      & ", c(paste("\\multicolumn{2}{c}{\\textit{Leaf ", 1:(length(unique(leaves))-1), "}} & ", sep = ""), paste("\\multicolumn{2}{c}{\\textit{Leaf ", length(unique(leaves)), "}}", sep = "")) ," \\\\",  paste0("\\cmidrule{", seq(2, length(unique(leaves)) * 2, by = 2), "-", seq(3, (length(unique(leaves)) * 2)+1, by = 2), "} "), "
     ", rep(" & Mean & (S.E.)", length(unique(leaves))), " \\\\
      \\addlinespace[2pt]
      \\hline \\\\[-1.8ex] \n\n", sep = "")

    for (i in seq_len(length(table_names))) {
      temp_means <- format(round(parms[[i]][, 1], 3), nsmall = 3)
      temp_sds <- format(round(parms[[i]][, 2], 3), nsmall = 3)

      temp_line <- paste0("      \\texttt{", table_names[i], "} & " )

      for (j in seq_len(length(unique(leaves)))) {
        temp_line <- paste0(temp_line, temp_means[j], " & (", temp_sds[j], ") & ")
      }

      temp_line <- stringr::str_sub(temp_line, end = -4)

      cat(temp_line, " \\\\ \n", sep = "")
    }

    cat("\n      \\addlinespace[3pt]
      \\\\[-1.8ex]\\hline
      \\hline \\\\[-1.8ex]
    \\end{tabular}
    \\end{adjustbox}
    \\caption{Average characteristics of units in each leaf, obtained by regressing each covariate on a set of dummies denoting leaf membership", if (!is.null(x$aggTree$idx$honest_idx)) "using only the honest sample." else "." ,"Standard errors are estimated via the Eicker-Huber-White estimator. Leaves are sorted in increasing order of the GATEs.}
    \\label{table_average_characteristics_leaves}
    \\end{table}
\\endgroup \n\n")

    ## Warn the user for zero variation in leaves.
    no_variation_names <- names(unlist(lapply(parms, function(x) { if (any(x[, 2] == 0)) idx <- 1})))

    if (length(no_variation_names) > 1) {
      names_string <- paste0(no_variation_names, collapse = ", ")
      warning(paste0("Variables '", names_string, "' have no variation in one or more leaves. Please correct the table by removing the associated standard errors."))
    } else if (length(no_variation_names) == 1) {
      warning(paste0("Variable '", no_variation_names, "' has no variation in one or more leaves. Please correct the table by removing the associated standard errors."))
    }
  } else if (table == "diff") {
    cat("\\begingroup
  \\setlength{\\tabcolsep}{8pt}
  \\renewcommand{\\arraystretch}{1.2}
  \\begin{table}[b!]
    \\centering
    \\begin{adjustbox}{width = 1\\textwidth}
    \\begin{tabular}{@{\\extracolsep{5pt}}l", rep(" c", times = get_leaves(x$groups)), "}
      \\\\[-1.8ex]\\hline
      \\hline \\\\[-1.8ex] \n
      & ", paste0("\\textit{Leaf ", seq_len(get_leaves(x$groups)-1), "} & "), paste0("\\textit{Leaf ", get_leaves(x$groups), "}") ," \\\\
      \\addlinespace[2pt]
      \\hline \\\\[-1.8ex] \n\n", sep = "")

    cat("      \\multirow{3}{*}{GATEs} & ", paste(gates_point[1:(length(unique(leaves))-1)], " & ", sep = ""), gates_point[length(unique(leaves))], " \\\\
      & ", paste("[", gates_ci_lower[1:(length(unique(leaves))-1)], ", ", gates_ci_upper[1:(length(unique(leaves))-1)], "] & ", sep = ""), paste("[", gates_ci_lower[length(unique(leaves))], ", ", gates_ci_upper[length(unique(leaves))], "]", sep = ""), " \\\\
      & ", paste("\\{", gates_ci_lower_boot[1:(length(unique(leaves))-1)], ", ", gates_ci_upper_boot[1:(length(unique(leaves))-1)], "\\} & ", sep = ""), paste("\\{", gates_ci_lower_boot[length(unique(leaves))], ", ", gates_ci_upper_boot[length(unique(leaves))], "\\}", sep = ""), " \\\\ \n\n", sep = "")
    cat("      \\addlinespace[2pt]
      \\hline \\\\[-1.8ex] \n\n")

    for (i in seq_len(get_leaves(x$groups))) {
      cat(paste0("      \\textit{Leaf ", i, "}"), " & ", paste0(format(round(x$gates_diff_pairs$gates_diff[i, seq_len(get_leaves(x$groups)-1)], 2), nsmall = 2), " & "), format(round(x$gates_diff_pairs$gates_diff[i, get_leaves(x$groups)], 2), nsmall = 2), " \\\\
            & ", paste0("(", format(round(x$gates_diff_pairs$holm_pvalues[i, seq_len(get_leaves(x$groups)-1)], 3), nsmall = 3), ") & "), paste0("(", format(round(x$gates_diff_pairs$holm_pvalues[i, get_leaves(x$groups)], 3), nsmall = 3), ")"), " \\\\ \n", sep = "")
    }

    cat("\n      \\addlinespace[3pt]
      \\\\[-1.8ex]\\hline
      \\hline \\\\[-1.8ex]
    \\end{tabular}
    \\end{adjustbox}
    \\caption{Point estimates and $95\\%$ confidence intervals for the GATEs based on asymptotic normality (in square brackets) and on the percentiles of the bootstrap distribution (in curly braces). Leaves are sorted in increasing order of the GATEs. Additionally, the GATE differences across all pairs of leaves are displayed. $p$-values testing the null hypothesis that a single difference is zero are adjusted using Holm's procedure and reported in parenthesis under each point estimate.}
    \\label{table_differences_gates}
    \\end{table}
\\endgroup \n\n")
  }
}
