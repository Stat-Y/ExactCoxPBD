#' Fitting Cox Model with Exact Partial Likelihood.
#'
#' This function fits a Cox model by optimizing the exact partial likelihood for the slope and the exact likelihood for the baseline hazard. It returns estimates for both the slope and baseline hazard, along with covariance estimates for the slope. Additionally, the function provides Wald and Likelihood Ratio Test (LRT) statistics and their p-values for inference. Refer to the manuscript for detailed calculations of the estimators and test statistics.
#'
#' @param formula a formula for the Cox model, following the same syntax as `coxph` from the `survival` package.
#' @param data the dataset for the Cox model, used in the same way as in `coxph` from the `survival` package.
#' @param baseline.data baseline data used to calculate the fitted baseline hazard of the initial Cox model. For a factor covariate, use the reference factor level; for other covariates, use zero. See the example code below for a detailed usage example.
#' @param bhf.initial tie correction method for the initial Cox model. The default is "efron"; other options align with the `method` argument in `coxph` from the `survival` package.
#' @param info.option method for computing the covariance estimator for the slope. The default option, "pbd", uses the inverse of the exact information matrix, computed numerically from the exact partial likelihood and evaluated at the PBD estimated slope. In contrast, the "breslow" option uses the inverse of the Breslow information, also evaluated at the PBD estimated slope. This choice influences both the covariance estimation and the Wald statistic.
#'
#' @import poibin
#' @import survival
#' @import numDeriv
#' @import KMsurv
#'
#' @export
#'
#' @return The returned list contains the following components:
#'
#' @return `coef` fitted slope optimizing exact partial likelihood.
#' @return `vcov` covariance estimator for the slope.
#' @return `lambda0` fitted baseline hazard, optimized based on exact likelihood.
#' @return `logPL` vector of exact log partial likelihood evaluated at the 0 vector and `coef`.
#' @return `lrt` likelihood ratio test statistic for the entire slope, including its degrees of freedom and p-value.
#' @return `wald.all` wald statistic for the entire slope, including its degrees of freedom and p-value.
#' @return `wald.each` wald statistic and p-value for each slope.
#' @return `bhf.initial` the given `bhf.initial`.
#' @return `info.option` the given `info.option`.
#' @return `n` sample size of the given `data`.
#' @return `nevent` total number of events in the given `data`.
#' @return `dat` a data object created by `coxph.pb.dat.setup` with the initial Cox model.
#' @return `formula` the given `formula`.
#' @return `call` the call used to run the function.
#'
#' @examples
#' library(poibin)
#' library(survival)
#' library(numDeriv)
#' library(KMsurv)
#'data("larynx")
#'fit_pbd=coxph.PB(Surv(time, delta)~age+diagyr+factor(stage),
#'                 data=larynx,
#'                 baseline.data=t(c(age=0,diagyr=0,stage=1)),
#'                 bhf.initial="efron",info.option="pbd")
#'print(fit_pbd)
#'
coxph.PB=function(formula, data, baseline.data, bhf.initial="efron", info.option="pbd")
{
  # Estimate beta
  row.names(data)=NULL

  environment(formula) <- environment()
  sfit=coxph(formula=formula, data=data, ties=c(bhf.initial))
  coefname=names(sfit$coefficients)

  bfit=survfit(sfit, newdata=data.frame(baseline.data))

  coxph.pb.dat=coxph.pb.dat.setup(cox.fit=sfit, bl.fit=bfit)
  beta0=coxph.pb.dat$beta0

  sfit1=try(lifetime.mle(dat=coxph.pb.dat, minusloglik=minus.log.partial.PB.discretized, starts=beta0), silent=T)

  try.flag1=(attr(sfit1,"class")=="try-error")

  if(length(try.flag1)==0 & length(coxph.pb.dat$lambda0)>1)
  {
    sfit=sfit1
  }

  names(sfit$coef)=coefname

  # Estimate lambda
  coxph.pb.dat$beta0=as.vector(sfit$coef)
  lambda0=lambda0.update(dat=coxph.pb.dat)
  names(lambda0)=sfit$dat$time
  sfit$lambda0=lambda0
  coxph.pb.dat$lambda0=lambda0

  # Estimate covariance when using breslow information option
  if(info.option=="breslow"){
    sfit$vcov=solve(infob(dat=sfit$dat,pars=sfit$coef))
  }
  colnames(sfit$vcov)=coefname
  rownames(sfit$vcov)=coefname

  # LRT
  logPL=c(null=-minus.log.partial.PB.discretized(coxph.pb.dat,rep(0,length(sfit$coef))),
          full=-minus.log.partial.PB.discretized(coxph.pb.dat,sfit$coef))
  LRT=2*(logPL["full"]-logPL["null"])
  lrt=data.frame(
    t(c(LRT=LRT,
        df=length(sfit$coef),
        pvalue=1-pchisq(q=LRT,df=length(sfit$coef))
    )))

  # Wald
  Wald=t(sfit$coef)%*%solve(sfit$vcov)%*%sfit$coef
  wald.all=data.frame(
    t(c(Wald=Wald,
        df=length(sfit$coef),
        pvalue=1-pchisq(q=Wald,df=length(sfit$coef))
    )))
  se_w=sqrt(diag(sfit$vcov))
  z_w=sfit$coef/se_w
  pval.w=1-pchisq(z_w^2,df=1)

  wald.each=data.frame(coef=sfit$coef,
                       se=se_w,
                       z=z_w,
                       pvalue=pval.w)
  row.names(wald.each)=coefname

  # Output
  output=list(coef=sfit$coef,
              vcov=sfit$vcov,
              logPL=logPL,
              lrt=lrt,
              wald.all=wald.all,
              wald.each=wald.each,
              bhf.initial=bhf.initial,
              info.option=info.option,
              n=dim(sfit$dat$xmat)[1],
              nevent=sum(sfit$dat$nevent),
              dat=sfit$dat,
              lambda0=sfit$lambda0,
              formula=formula,
              call=sfit$call)

  return(output)
}
