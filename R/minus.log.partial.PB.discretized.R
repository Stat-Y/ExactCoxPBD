#' Calculate the Minus Log Exact Partial Likelihood.
#'
#' This function calculates the negative log of the exact partial likelihood for a given slope. Since the true baseline hazard is unknown, it is substituted with a fitted baseline hazard from an initial Cox model.
#'
#' @param dat the basic data setup result by `coxph.pb.dat.setup`, where the `lambda0`: fitted baseline hazard from the initial Cox model in it, substitutes the unknown true baseline hazard in the minus log exact partial likelihood.
#' @param pars given slope.
#'
#' @import poibin
#' @import survival
#' @import numDeriv
#' @import KMsurv
#'
#' @export
#'
#' @return The returned value is minus log partial likelihood at a given slope.
#'
#' @examples
#' library(KMsurv)
#' data("larynx")
#' cox.fit=coxph(Surv(time, delta)~age+diagyr+factor(stage),data=larynx,method="efron")
#' bl.fit=survfit(cox.fit, newdata=data.frame(t(c(age=0,diagyr=0,stage=1))))
#' pb.data.setup=coxph.pb.dat.setup(cox.fit,bl.fit)
#' minus.log.PL=minus.log.partial.PB.discretized(pb.data.setup,rep(0,5))
#' print(minus.log.PL)
#'
minus.log.partial.PB.discretized=function(dat, pars)
{
  failmat=dat$failmat
  riskmat=dat$riskmat
  xmat=dat$xmat
  nevent=dat$nevent
  lambda0=dat$lambda0

  risk.score=exp(xmat%*%pars)
  risk.score=as.vector(risk.score)

  tmp=kronecker(risk.score, t(lambda0), "*")

  prob.mat=1-exp(-tmp)

  ##numeric correction
  prob.mat[prob.mat<=(1e-5)]=1e-5
  prob.mat[prob.mat>=(1-(1e-5))]=1-(1e-5)
  ###

  prob.mat=prob.mat*riskmat

  ndt=ncol(riskmat)

  res=double(ndt)


  for(i in 1:ndt)
  {
    #browser()

    idx1=(failmat[,i]==1)
    idx2=(riskmat[,i]==1)

    pp=prob.mat[idx2, i]

    fpp=prob.mat[idx1&idx2, i]
    npp=prob.mat[(!idx1)&idx2, i]
    di=nevent[i]

    tres1=sum(log(fpp))
    tres2=0
    if(length(npp)>=1)
    {
      tres2=sum(log(1-npp))
    }

    #print(pp)

    tmp2=dpoibin(kk=di, pp=pp)

    if(tmp2<=(1e-7))
    {
      tmp2=(1e-7)
    }

    tres3=log(tmp2)

    tres=tres1+tres2-tres3

    res[i]=(-1)*tres
  }


  val=sum(res)

  #browser()

  return(val)

}
