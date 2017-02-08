# laggedLGCP
R package for fitting multivariate Poisson or Coalescent processes with a time lag.

## Installation
1. Install `phylodyn` as detailed at [mdkarcher/phylodyn](https://github.com/mdkarcher/phylodyn)
2. Install `rstan` -- http://mc-stan.org/interfaces/rstan.
3. Install `rstantools` via CRAN:
```{r}
install.packages("rstantools")
```
4. Install `laggedLGCP` via 
```{r}
devtools::install_github("AlexanderVR/laggedLGCP", args="--preclean", build_vignettes = FALSE)
```
