# Aggregation Trees <a href="https://riccardo-df.github.io/aggTrees/"><img src="man/figures/logo.svg" align="right" height="200" /></a>

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/aggTrees)](https://CRAN.R-project.org/package=aggTrees)
![CRAN Downloads overall](http://cranlogs.r-pkg.org/badges/grand-total/aggTrees)
[![R-CMD-check](https://github.com/riccardo-df/aggTrees/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/riccardo-df/aggTrees/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

R package to implement aggregation trees, a nonparametric approach to discovering heterogeneous subgroups in a selection-on-observables framework. 

`aggTrees` allows researchers to assess whether there exists relevant heterogeneity in treatment effects by generating a sequence of optimal groupings, one for each level of granularity. For each grouping, we obtain point estimation and inference about the group average treatment effects. Please reference the use as [Di Francesco (2022)](https://doi.org/10.2139/ssrn.4304256).

To get started, please check the online [short tutorial](https://riccardo-df.github.io/aggTrees/articles/aggTrees-vignette.html).

## Installation  
The package can be downloaded from CRAN:

```
install.packages("aggTrees")
```

Alternatively, the current development version of the package can be installed using the `devtools` package:

```
devtools::install_github("riccardo-df/aggTrees") # run install.packages("devtools") if needed.
```

## References

- Athey, S., & Imbens, G. W. (2016).
<b>Recursive Partitioning for Heterogeneous Causal Effects.</b>
<i>Proceedings of the National Academy of Sciences</i>, 113(27).
[<a href="https://doi.org/10.1073/pnas.1510489113">paper</a>]

- Athey, S., Tibshirani, J., & Wager, S. (2019).
<b>Generalized Random Forests.</b> 
<i>Annals of Statistics</i>, 47(2).
[<a href="https://doi.org/10.1214/18-AOS1709">paper</a>]

- Chernozhukov, V., Demirer, M., Duflo, E., & Fernandez-Val, I. (2017).
<b>Generic Machine Learning Inference on Heterogeneous Treatment Effects in Randomized Experiments.</b>
<i>arXiv preprint</i>.
[<a href="https://doi.org/10.48550/arXiv.1712.04802">paper</a>]

- Cotterman, R., & Peracchi, F. (1992).
<b>Classification and aggregation: An application to industrial classification in cps data.</b> 
<i>Journal of Applied Econometrics</i>, 7(1).
[<a href="https://doi.org/10.1002/jae.3950070105">paper</a>]

- Di Francesco, R. (2022).
<b>Aggregation Trees.</b> 
<i>CEIS Research Paper, 546.</i>
[<a href="https://doi.org/10.2139/ssrn.4304256">paper</a>]

- Holm, S. (1979).
<b>A Simple Sequentially Rejective Multiple Test Procedure.</b> 
<i>Scandinavian Journal of Statistics</i>, 6(2).

- Semenova, V., & Chernozhukov, V. (2021).
<b>Debiased Machine Learning of Conditional Average Treatment Effects and Other Causal Functions.</b>
<i>The Econometrics Journal</i>, 24(2).
[<a href="https://doi.org/10.1093/ectj/utaa027">paper</a>]
