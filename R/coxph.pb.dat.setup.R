#' Basic Data Setup Prior to Computing the Minus Log Exact Partial Likelihood
#'
#' The exact partial likelihood involves a complex function of conditional probabilities based on the risk and event sets. This function provides essential data setup, such as distinct ordered event times, the count of individuals at risk and events at each time, risk and event sets, and the covariate matrix. These components serve as the groundwork for calculating the negative log of the exact partial likelihood. Since the exact partial likelihood depends on an unknown true baseline hazard, we substitute it with an initial Cox model, like Breslow or Efron, using its fitted baseline hazard. The output also includes the fitted slope and the baseline hazard from this initial Cox model.
#'
#' @param cox.fit fitted Cox model by `coxph` from the `survival` package.
#' @param bl.fit fitted baseline hazard of `cox.fit` by `survfit` from the `survival` package.
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
#' @return `time`: a vector of event times.
#' @return `nevent`: a vector representing the number of events at each event time.
#' @return `nrisk`: a vector representing the number of units at risk at each event time.
#' @return `failmat`: a matrix of event indicators, where rows represent subjects and columns represent event times.
#' @return `riskmat`: a matrix of at-risk indicators, where rows represent subjects and columns represent event times.
#' @return `xmat`: a covariate matrix, where rows represent subjects and columns represent covariates.
#' @return `beta0`: the fitted slope from the initial Cox model.
#' @return `lambda0`: the fitted baseline hazard from the initial Cox model at each event time.
#'
#' @examples
#' library(KMsurv)
#' data("larynx")
#' cox.fit=coxph(Surv(time, delta)~age+diagyr+factor(stage),data=larynx,method="efron")
#' bl.fit=survfit(cox.fit, newdata=data.frame(t(c(age=0,diagyr=0,stage=1))))
#' pb.data.setup=coxph.pb.dat.setup(cox.fit,bl.fit)
#' print(pb.data.setup)
#'
#'
coxph.pb.dat.setup=function(cox.fit, bl.fit)
{
  beta0=as.vector(cox.fit$coef)
  xmat=model.matrix(cox.fit)

  tmp=coxph.detail(object=cox.fit, riskmat=T, rorder=c("data", "time"))
  time=tmp$time
  nevent=tmp$nevent
  nrisk=tmp$nrisk
  riskmat=tmp$riskmat

  failmat=riskmat
  yy=tmp$y
  yy=yy[order(as.numeric(rownames(yy))),]

  for(i in 1:nrow(yy))
  {
    for(j in 1:ncol(failmat))
    {
      failmat[i,j]=0

      if(yy[i, "time"]==time[j] & yy[i, "status"]==1)
      {
        failmat[i,j]=1
      }

    }
  }

  lambda0=diff(c(0,bl.fit$cumhaz[bl.fit$n.event>0]))

  res=list(time=time,
           nevent=nevent,
           nrisk=nrisk,
           failmat=failmat,
           riskmat=riskmat,
           xmat=xmat,
           beta0=beta0,
           lambda0=lambda0)

  return(res)

}
