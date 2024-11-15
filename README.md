### ExactCoxPBD R Package

The *ExactCoxPBD* R package enables fitting the Cox model by optimizing the exact partial likelihood as proposed in the paper, *"An Accurate Computational Approach for Partial Likelihood Using Poisson-Binomial Distributions"* by Cho et al. This method is especially useful for data with heavy ties or high covariate variation, providing improved accuracy in Cox model estimation.

### R Version Requirements

Requires R version 4.3.2 or higher.

### Installing the Package

Ensure the following R packages are installed: `poibin`, `survival`, `numDeriv`, `KMsurv`, and `devtools`. Then, install the *ExactCoxPBD* package by running:

```r
devtools::install_github("Stat-Y/ExactCoxPBD")
```

### Help

To view the help file for the main function `coxph.PB`, use `?coxph.PB`.

### Usage

```r
coxph.PB(
  formula,
  data,
  baseline.data,
  bhf.initial = "efron",
  info.option = "pbd"
)
```

### Arguments

*formula*: A formula specifying the Cox model, following the syntax of `coxph` from the *survival* package.
  
*data*: The dataset for the Cox model, used in the same format as in `coxph` from the *survival* package.

*baseline.data*: Baseline data used to calculate the fitted baseline hazard of the initial Cox model. For factor covariates, specify the reference factor level; for continuous covariates, use zero.

*bhf.initial*: Tie correction method for the initial Cox model. The default is `"efron"`, with other options matching those available in the `method` argument of `coxph` from the *survival* package.

*info.option*: Method for computing the covariance estimator for the slope. The default, `"pbd"`, uses the inverse of the exact information matrix, computed numerically from the exact partial likelihood and evaluated at the PBD estimated slope. Alternatively, `"breslow"` uses the inverse of the Breslow information, also evaluated at the PBD estimated slope. This choice affects both the covariance estimation and the Wald statistic.

### Value

The function returns a list containing:

*coef*: Fitted slope based on optimizing the exact partial likelihood.
  
*vcov*: Covariance estimator for the slope.
  
*lambda0*: Fitted baseline hazard, optimized based on the exact likelihood.
  
*logPL*: Vector of exact log partial likelihood evaluated at the 0 vector and `coef`.
  
*lrt*: Likelihood ratio test statistic for the entire slope, including its degrees of freedom and p-value.
  
*wald.all*: Wald statistic for the entire slope, including degrees of freedom and p-value.
  
*wald.each*: Wald statistic and p-value for each slope.
  
*bhf.initial*: The specified `bhf.initial` method.
  
*info.option*: The specified `info.option` method.
  
*n*: Sample size of the dataset.
  
*nevent*: Total number of events in the dataset.
  
*dat*: Data object created by `coxph.pb.dat.setup` with the initial Cox model.
  
*formula*: The specified formula.
  
*call*: The function call.

### Example

The following example demonstrates fitting a Cox model based on the exact partial likelihood:

```r
library(poibin)
library(survival)
library(numDeriv)
library(KMsurv)

data("larynx")
fit_pbd <- coxph.PB(
  formula = Surv(time, delta) ~ age + diagyr + factor(stage),
  data = larynx,
  baseline.data = t(c(age = 0, diagyr = 0, stage = 1)),
  bhf.initial = "efron",
  info.option = "pbd"
)
print(fit_pbd)
```
