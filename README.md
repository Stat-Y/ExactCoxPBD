### ExactCoxPBD R Package

The *ExactCoxPBD* R package enables fitting the Cox model by optimizing the exact partial likelihood as proposed in the paper, *"An Accurate Computational Approach for Partial Likelihood Using Poisson-Binomial Distributions"* by Cho et al. This method is especially useful for data with heavy ties or high covariate variation, providing improved accuracy in Cox model estimation.

### R Version Requirements

Requires R version 4.3.2 or higher.

### Installing the Package

Ensure the following R packages are installed: `poibin`, `survival`, `numDeriv`, `KMsurv`, and `devtools`. Then, install the *ExactCoxPBD* package by running:

```r
devtools::install_github("Stat-Y/ExactCoxPBD")
```



