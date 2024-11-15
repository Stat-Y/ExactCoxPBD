#' Compute the Minimizer for a Given Loss Function.
#'
#' This function computes the minimizer for a given loss function, which can be the negative log of the exact partial likelihood for slope estimation or the negative log of the exact likelihood for baseline hazard estimation.
#'
#' @param dat the basic data setup result by `coxph.pb.dat.setup`.
#' @param minusloglik loss function to minimize, which depends on `dat`.
#' @param starts starting value for the optimization used in `optim` inside the function.
#' @param method the optimization method used in `optim` inside the function, with the default set to "BFGS".
#' @param hessian the Hessian option used in `optim` inside the function, with the default set to TRUE.
#' @param ... additional arguments for `optim` inside the function.
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
#' @return `call`: the call used in running the function.
#' @return `coef`: the fitted coefficient, which is the optimizer of the given loss function.
#' @return `vcov`: inverse of the hessian for the given loss function evaluated at `coef`.
#' @return `min`: the value of the given loss function evaluated at `coef`.
#' @return `dat`: the given `dat`.
#' @return `minusloglik`: the given `minusloglik`.
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
#' print(pb.slope.fit)
#'
lifetime.mle=function(dat, minusloglik, starts, method = "BFGS",hessian = TRUE,...)
{
  call=match.call()
  f = function(p) {
    minusloglik(dat,p)
  }
  oout = optim(starts, f, method = method, hessian = hessian,...)#,control=list(trace=T))
  coef = oout$par
  #browser()
  if(hessian)
  {
    vcov =solve(oout$hessian)
  }else{
    vcov=NULL
  }
  min = oout$value
  invisible(list(call = call, coef = coef,vcov = vcov, min = min,dat=dat,minusloglik=minusloglik))
}
