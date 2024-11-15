#' Estimate the Baseline Hazard Minimizing the Negative Log Exact Likelihood.
#'
#' This function estimate the baseline hazard by minimizing the negative log of the exact likelihood for the baseline hazard. Since the true slope is unknown, it is substituted with the PBD estimated slope.
#'
#' @param dat a data object created by `coxph.pb.dat.setup`. In this object, the `beta0` represents the PBD estimated slope, rather than an initial Cox model slope. To substitute the initial slope with the PBD estimated slope, refer to the example code provided below.
#'
#' @import poibin
#' @import survival
#' @import numDeriv
#' @import KMsurv
#'
#' @export
#'
#' @return The returned value is the estimated baseline hazard.
#'
#' @examples
#' library(poibin)
#' library(survival)
#' library(numDeriv)
#' library(KMsurv)
#' data("larynx")
#' cox.fit=coxph(Surv(time, delta)~age+diagyr+factor(stage),data=larynx,method="efron")
#' bl.fit=survfit(cox.fit, newdata=data.frame(t(c(age=0,diagyr=0,stage=1))))
#' pb.data.setup=coxph.pb.dat.setup(cox.fit,bl.fit)
#' pb.slope.fit=lifetime.mle(dat=pb.data.setup, minusloglik=minus.log.partial.PB.discretized, starts=pb.data.setup$beta0)
#' pb.data.setup$beta0=as.vector(pb.slope.fit$coef)
#' pb.baseline.fit=lambda0.update(dat=pb.data.setup)
#' print(pb.baseline.fit)
#'
lambda0.update=function(dat)
{
  # Update baseline hazard

  failmat=dat$failmat
  riskmat=dat$riskmat
  xmat=dat$xmat
  nevent=dat$nevent
  lambda0=dat$lambda0
  beta0=dat$beta0

  risk.score=exp(xmat%*%beta0)
  risk.score=as.vector(risk.score)

  ndt=ncol(riskmat)

  res=double(ndt)

  for(i in 1:ndt)
  {
    idx1=(failmat[,i]==1)
    idx2=(riskmat[,i]==1)

    npp=risk.score[(!idx1)&idx2]
    if(length(npp)==0)
    {
      res[i]=1e6
    }else{
      tmp.dat=list(risk.score=risk.score, idx1=idx1, idx2=idx2, init=lambda0[i])
      res[i]=lambda0.i.update(dat=tmp.dat)
    }
  }

  return(res)

}

lambda0.i.update=function(dat)
{
  # Updates baseline hazard at each event time

  risk.score=dat$risk.score
  idx1=dat$idx1
  idx2=dat$idx2
  init=dat$init

  minus.obj.fun.lambda0.i.update(dat=dat, pars=log(init))

  ifit=lifetime.mle(dat=dat, minusloglik=minus.obj.fun.lambda0.i.update, starts=log(init))

  res=exp(ifit$coef)
  res=as.vector(res)

  return(res)
}

minus.obj.fun.lambda0.i.update=function(dat, pars)
{
  # Minus log exact likelihood for baseline hazard at each event time

  risk.score=dat$risk.score
  idx1=dat$idx1
  idx2=dat$idx2

  lambda0.i=exp(pars)
  risk.score=lambda0.i*risk.score

  tmp1=risk.score[idx1&idx2]
  tmp2=risk.score[(!idx1)&idx2]

  res=(-1)*sum(log(1-exp(-tmp1)))-sum(-tmp2)

  return(res)

}
