### ExactCoxPBD R Package

The `ExactCoxPBD` R package enables fitting the Cox model by optimizing the exact partial likelihood as proposed in the paper, *"An Accurate Computational Approach for Partial Likelihood Using Poisson-Binomial Distributions"* by Cho et al. (2024+). This method is especially useful for data with heavy ties or high covariate variation, providing improved accuracy in Cox model estimation.

### R Version Requirements

Requires R version 4.3.2 or higher.

### Installing the Package

Ensure the following R packages are installed: `poibin`, `survival`, `numDeriv`, `KMsurv`, and `devtools`. Then, install the `ExactCoxPBD` package by running:

```r
if (!require("poibin")) install.packages("poibin")
if (!require("survival")) install.packages("survival")
if (!require("numDeriv")) install.packages("numDeriv")
if (!require("KMsurv")) install.packages("KMsurv")
if (!require("devtools")) install.packages("devtools")
if (!require("ExactCoxPBD")) devtools::install_github("Stat-Y/ExactCoxPBD")
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

- `formula`: a formula for the Cox model, following the same syntax as `coxph` from the `survival` package.
  
- `data`: the dataset for the Cox model, used in the same way as in `coxph` from the `survival` package.

- `baseline.data`: baseline data used to calculate the fitted baseline hazard of the initial Cox model. For a factor covariate, use the reference factor level; for other covariates, use zero. See the example code below for a detailed usage example.

- `bhf.initial`: tie correction method for the initial Cox model. The default is `"efron"`; other options align with the `method` argument in `coxph` from the `survival` package.

- `info.option`: method for computing the covariance estimator for the slope. The default option, `"pbd"`, uses the inverse of the exact information matrix, computed numerically from the exact partial likelihood and evaluated at the PBD estimated slope. In contrast, the `"breslow"` option uses the inverse of the Breslow information, also evaluated at the PBD estimated slope. This choice influences both the covariance estimation and the Wald statistic.

### Value

The function returns a list containing:

- `coef`: fitted slope optimizing exact partial likelihood.

- `vcov`: covariance estimator for the slope.

- `lambda0`: fitted baseline hazard, optimized based on exact likelihood.

- `logPL`: vector of exact log partial likelihood evaluated at the 0 vector and `coef`.

- `lrt`: likelihood ratio test statistic for the entire slope, including its degrees of freedom and p-value.

- `wald.all`: wald statistic for the entire slope, including its degrees of freedom and p-value.

- `wald.each`: wald statistic and p-value for each slope.

- `bhf.initial`: the given `bhf.initial`.

- `info.option`: the given `info.option`.

- `n`: sample size of the given `data`.

- `nevent`: total number of events in the given `data`.

- `dat`: a data object created by `coxph.pb.dat.setup` with the initial Cox model.

- `formula`: the given `formula`.

- `call`: the call used to run the function.

### Example

The following example fits a Cox model based on the exact partial likelihood:

```r
library(poibin)
library(survival)
library(numDeriv)
library(KMsurv)
library(ExactCoxPBD)
data("larynx")
fit_pbd=coxph.PB(Surv(time, delta)~age+diagyr+factor(stage),
                data=larynx,
                baseline.data=t(c(age=0,diagyr=0,stage=1)),
                bhf.initial="efron",info.option="pbd")
print(fit_pbd)
```
