##' Approximates the models' posterior weights  by simple Monte Carlo integration
##'
##' if \eqn{J} is the number of models, the posterior weights are given by
##'\deqn{postW(i) = priorW(i)*lkl(i)/
##' (\sum_{j=1,\dots,J} priorW(j)*lkl(j)),}
##' where \eqn{lkl(i)} stands for the Monte-Carlo estimate of the
##' marginal likelihood of model \eqn{i} and \eqn{priorW(i)} is
##' the prior weight
##' defined in \code{priorweights[i]}. For more explanations, see the reference below 
##' The confidence intervals are obtained by adding/subtracting two times the estimated standard errors of the marginal likelihood estimates.
##' The latter are only estimates, which interpretation may be misleading:
##' See the note in section \code{\link{marginal.lkl}}
##' @title Posterior model weights
##' @inheritParams marginal.lkl
##' @inheritParams posteriorMCMC
##' @param HparList A list containing the hyper Parameter for the priors in each model. (list of lists).
##' @param lklList A list containing the likelihood functions of each model
##' @param priorList A list containing the prior definitions of each model.
##' @param priorweights A vector of positive weights, summing to one: the prior marginal weights of each model.
##' @param Nsim The maximum number of iterations to be performed.
##' @param displ Logical. Should convergence monitoring plots be issued ?
##' @references HOETING, J., MADIGAN, D., RAFTERY, A. and VOLINSKY, C. (1999). Bayesian model averaging: A tutorial. \emph{Statistical science 14, 382-401}.
##' @return A matrix of \eqn{6} columns and \code{length(priorweights)} rows. The columns contain respectively the posterior model weights (in the same order as in \code{priorweights}), the lower and the upper bound of the confidence interval (see \bold{Details}), the marginal model weights, the estimated standard error of the marginal likelihood estimators, and the number of simulations performed.
##'@examples
##'data(pb.Hpar)
##' data(nl.Hpar)
##' set.seed(5)
##' mixDat=rbind(rpairbeta(n=10,dimData=3, par=c(0.68,3.1,0.5,0.5)),
##'   rnestlog(n=10,par=c(0.68,0.78, 0.3,0.5)))
##' posteriorWeights (dat=mixDat,
##'                   HparList=list(get("pb.Hpar"),get("nl.Hpar")),
##'                   lklList=list(get("dpairbeta"), get("dnestlog")),
##'                   priorList=list(get("prior.pb"), get("prior.nl")),
##'                   priorweights=c(0.5,0.5),
##'                   Nsim=1e+3,
##'                   Nsim.min=5e+2, precision=0.1,
##'                   displ=FALSE)
##' \dontrun{posteriorWeights (dat=mixDat,
##'                   HparList=list(get("pb.Hpar"),get("nl.Hpar")),
##'                   lklList=list(get("dpairbeta"), get("dnestlog")),
##'                   priorList=list(get("prior.pb"), get("prior.nl")),
##'                   priorweights=c(0.5,0.5),
##'                   Nsim=20e+3,
##'                   Nsim.min=10e+3, precision=0.05,
##'                   displ=TRUE)}
##' @export
posteriorWeights <-
  function(dat,
           HparList=list(get("pb.Hpar"),get("nl.Hpar")),
           lklList=list(get("dpairbeta"), get("dnestlog")),
           priorList=list(get("prior.pb"), get("prior.nl")),
           priorweights=c(0.5,0.5),
           Nsim=20e+3,
           Nsim.min=10e+3, precision=0.05, 
           seed=1, kind = "Mersenne-Twister",
           show.progress=floor(seq(1,Nsim,length.out=10)),
           displ=FALSE)
  {
    if( (length(HparList) != length(lklList)) | (length(lklList) != length(priorweights) ) | (length(priorweights) != length(priorList) ) )
      {stop("arguments 'HparList', 'lklList', 'priorList'  and 'priorweights' should have  same length. ")}

    Nmodel=length(HparList)

    output=matrix(0,ncol=6,nrow=Nmodel,
      dimnames=list(NULL,
        c("postW", "lowConf", "upConf", "margLkl", "estErr", "nsim"))
      )

    
    
    for(i in 1:Nmodel)
      {
        set.seed(seed, kind =kind)
        
        lklMarg=marginal.lkl(dat=dat,
          likelihood=lklList[[i]],
          prior=priorList[[i]],
          Nsim=Nsim ,
          displ= displ,
          Hpar=HparList[[i]],
          Nsim.min=Nsim.min,
          show.progress=show.progress,
          precision=precision
          )
        if(displ)
          {
            mtext(paste("Model ", toString(i),
                        " : marginal likelihood",
                        sep=""), side=3, line=2)
          }
        
        output[i,4]=lklMarg$emp.mean
        output[i,5]=lklMarg$est.error
        output[i,6]=lklMarg$nsim
      }

    lklTot= sum(output[,4] * priorweights)
    
    lklTotSup=sum((output[,4]+ 2*output[,5]) * priorweights) 

    lklTotInf=sum((output[,4]- 2*output[,5]) * priorweights) 

    output[,1]=(output[,4]*priorweights)/lklTot

    output[,2] =priorweights* (output[,4]-2*output[,5] )/
      (lklTotSup -4*priorweights*output[,5] )

    output[,3] =priorweights* (output[,4]+2*output[,5] )/
      (lklTotInf + 4*priorweights*output[,5] )

    return(output)

  }
    
