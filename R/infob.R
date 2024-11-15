#' Calculate the Breslow Information Matrix.
#'
#' This function calculates the Breslow information matrix for a given slope.
#'
#' @param dat the basic data setup result by `coxph.pb.dat.setup`, where `nevent`, `riskmat`, and `xmat` in it are used in Breslow information calculation.
#' @param pars given slope.
#'
#' @import poibin
#' @import survival
#' @import numDeriv
#' @import KMsurv
#'
#' @export
#'
#' @return The returned value is the Breslow information matrix at a given slope.
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
#' br.info=infob(pb.data.setup,rep(0,5))
#' print(br.info)
#'
infob=function(dat,pars){
  pp=length(pars)
  info=matrix(0,pp,pp)
  for(jj in 1:pp){
    for(kk in 1:pp){
      risk.score=dat$xmat%*%as.vector(pars)
      xtmpj=matrix(dat$xmat[,jj])%*%rep(1,length(dat$nevent))
      xtmpk=matrix(dat$xmat[,kk])%*%rep(1,length(dat$nevent))
      tmp=risk.score%*%rep(1,length(dat$nevent))
      den=apply(exp(tmp)*dat$riskmat,2,sum)
      nunj=apply(xtmpj*exp(tmp)*dat$riskmat,2,sum)
      nunk=apply(xtmpk*exp(tmp)*dat$riskmat,2,sum)
      nunjk=apply(xtmpj*xtmpk*exp(tmp)*dat$riskmat,2,sum)
      info[jj,kk]=sum((nunjk/den-(nunj/den)*(nunk/den))*dat$nevent)
    }
  }
  return(info)
}
