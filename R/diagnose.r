##' @rdname diagnose.PBNLpostsample
##' @export
diagnose <- function(obj, ...)
  UseMethod("diagnose")



## @param Hpar The hyper parameter list. 
## @param dat The angular dataset
## @param model one of the character strings \code{"pairbeta"} or\code{"nestlog"}
## @param save Logical. Should the result be saved ?
## @param name.save The name under whic hthe result is to be saved. 
## @param save.directory The directory where the result is to be saved, without trailing slash.


##' The method  issues several convergence diagnostics, in the particular case when the PB or the NL model is used. The code may be easily modified for other angular models. 
##'
##' @title Diagnostics for the MCMC output in the PB and NL models.
##' @param obj an object of class \code{postsample}:  posterior sample, as produced by
##' \code{\link{posteriorMCMC.pb}} or \code{\link{posteriorMCMC.nl}}
##' @param true.par The true parameter. If \code{NULL}, it is considered as unknown.
##' @inheritParams discretize
##' @inheritParams add.frame
##' @inheritParams posterior.predictive3D
##' @param autocor.max The maximum accepted auto-correlation  for two successive parameters in the thinned sample.
##' @param default.thin The default thinning interval if the above condition cannot be satisfied.
##' @param predictive Logical. Should the predictive density be plotted ?
##' @param xlim.density The \code{xlim} interval for the density plots,
##' on the transformed scale.
##' @param ylim.density the \code{ylim} intervals for the density plots.
##' @param plot Logical. Should plots be issued ?
##' @param save Logical: should the result be saved ? Only used if the posterior sample has been saved itself (\emph{i.e.} if it contains \code{save=TRUE} in its arguments list)
##' @param ... Additional parameters to be passed to the functions
##' \code{\link{posterior.predictive.pb}} or \code{\link{posterior.predictive.nl}}.
##' @return A list made of \describe{
##' \item{predictive}{The posterior predictive, or \code{0} if \code{predictive=FALSE} }
##' \item{effective.size}{the effective sample size of each component}
##' \item{heidelTest}{The first part of the Heidelberger and Welch test (stationarity test). The first row indicates \dQuote{success} (1) or
##' rejection(0), the second line shows the number of iterations to be discarded, the third line is the p-value of the test statistic.}
##' \item{gewekeTest}{The test statistics from the Geweke stationarity test.}
##' \item{gewekeScore}{The p-values for the above test statistics}
##' \item{thin}{The thinning interval retained}
##' \item{correl.max.thin}{The maximum auto-correlation for a lag equal to \code{thin} }
##' \item{linked.est.mean}{The posterior mean of the transformed parameter (on the real line)}
##' \item{linked.est.sd}{The standard deviation of the transformed parameters}
##' \item{est.mean}{The posterior mean of the original parameters, as they appears in the expression of the likelihood}
##' \item{sample.sd}{the posterior standard deviation of the original parameters}
##' }
##' @export
##' @method diagnose PBNLpostsample
diagnose.PBNLpostsample <- function(obj,
                                    true.par=NULL,
                                    from =NULL, 
                                    to=NULL,    
                                    autocor.max = 0.2,
                                    default.thin=50,
                                    xlim.density=c(-4,4),
                                    ylim.density=NULL,
                                    plot=TRUE,
                                    predictive=FALSE,
                                    save=TRUE,
                                    ## npoints=60, eps=10^(-3),
                                    ## equi=TRUE,
                                    ## save=FALSE,
                                    ## name.save="pb.predictive",
                                    ## save.directory = "~",
                                    ...
                                    )
  {
    dat <- obj$arguments$dat
    Hpar <- obj$arguments$Hpar
    model <- obj$arguments$name.model

    likelihood <- obj$arguments$likelihood
    prior <- obj$arguments$prior

    
    save.directory <- obj$arguments$save.directory
    if(model=="pairbeta")
      {
        link <- function(x){log(x)}
        unlink <- function(x){exp(x)}
       
      }
    if(model=="nestlog")
      {
        link <- logit
        invlink <- invlogit
      }
    orig.vals <- obj$stored.vals
    vals <- link(obj$stored.vals)
    Nbin <- obj$arguments$Nbin
    Nsim <- obj$arguments$Nsim
    if(is.null(to))
      to <- Nsim
     
   if(is.null(from))
     from <- Nbin+1
      
 
    vals <- vals[(from-Nbin):(to-Nbin),]
    orig.vals <- orig.vals[(from-Nbin):(to-Nbin),]

    autocorr <- acf(vals, lag.max=default.thin, plot=plot) $ acf
    goodLags <- which( apply(abs(autocorr),1,max) < autocor.max)
  
    if(length(goodLags) == 0)
      thin <- default.thin
    
    else
      thin <- min(goodLags)
      
    correl.thin <- max(autocorr[thin,,])
    
    effective.size <- effectiveSize(vals)
    
    heidelTest <- apply(vals,2,heidel.diag, eps=0.1, pvalue= 0.05)[1:3,]
    gewekeTest <- geweke.diag(vals)
    gewekeScore <- pnorm(abs(gewekeTest[[1]]), lower.tail=FALSE)
   
   
    linked.cum.mean <- apply(vals,2,cumsum)/(1:(to-from+1))
    cum.mean <- apply(orig.vals,2,cumsum)/(1:(to-from+1))
    linked.est.mean <- linked.cum.mean[to-from+1,]
    est.mean <- cum.mean[to-from+1,]
    linked.est.sd <- apply(vals,2,sd)
    est.sd <- apply(orig.vals,2,sd)
    ## est.error.effsize = apply(orig.vals,2,function(X){
    ## summary(mcmc(X) )$statistics[4]})

    if(plot)
      {
        for(i in 1:ncol(vals))
          {
            dev.new()
            plot(from:to, linked.cum.mean[,i],type ="l",
                 main=paste("par", toString(i), sep=" "), ylab="") 
          }
        
    
        dprior.Talpha <-  function(x)
          {
            return( dnorm(x,
                          mean=Hpar$mean.alpha,
                          sd=Hpar$sd.alpha) )
          }
        dprior.Tbeta <-  function(x)
            {
              return( dnorm(x,
                            mean=Hpar$mean.beta,
                            sd=Hpar$sd.beta) )
            }
        XX <- seq(xlim.density[1], xlim.density[2], length.out=100)
        YY.alpha <- sapply(XX,dprior.Talpha)
        YY.beta <- sapply(XX,dprior.Tbeta)

        dev.new()
        par(mfrow = c(ceiling(sqrt(ncol(vals))),
              ceiling(ncol(vals) / ceiling(sqrt(ncol(vals))) )))
    
        for(i in 1:ncol(vals))
          {
    
            plot(density(vals[,(i)]),
                 main="",
                 xlab="",ylab="",
                 ylim=ylim.density, xlim=xlim.density ,# add=FALSE,
                 lty=1,lwd=2,col="black")

            if(i==1){
              lines(XX,YY.alpha,lty=3, lwd=1.5, col="black")
              title(main=switch(model, pairbeta="log(alpha)",
                      nestlog="logit(alpha)"))
            }
            else{
              lines(XX,YY.beta,lty=3, lwd=1.5, col="black")
              title(main=switch(model,
                      pairbeta= paste("log(beta[",i-1,"])", sep=""),
                      nestlog=paste("logit(beta[",i-1,"])", sep="")) )
            }
            if( ! is.null(true.par) )
              abline(v=link(true.par[i]),col="black",lwd=2, lty=2)
          }
      }
        time.idx = (from-Nbin):(to-Nbin)
        kept.idx = time.idx[ (time.idx %% thin == 0) ]
    
    if(predictive)
      {
        if(plot)
          dev.new()
        
        if(model == "pairbeta")
          predictiveDens <- posterior.predictive.pb(
            post.sample = obj,
            thin = thin,  displ=plot,
            ...
            )
        if(model=="nestlog")
          predictiveDens <- posterior.predictive.nl(
            post.sample = obj,
            thin = thin,  displ=plot,
            ...
            )
 
      }
    if(predictive)
      {
        pred.return <- predictiveDens
      }
    else
      {
        pred.return <- NA
      }
    result.list <- list(predictive = pred.return,
                        size.subsample.predictive=  length(kept.idx ),
                        effective.size=effective.size,
                        heidelTest = heidelTest,
                        gewekeTest=gewekeTest,
                        gewekeScore=gewekeScore,
                        thin = thin,
                        from=from,to=to,
                        thinnedCorrel=autocorr[thin,,],
                        correl.max.thin=correl.thin,
                        linked.est.mean=linked.est.mean,
                        linked.sample.sd=linked.est.sd,
                        est.mean=est.mean,
                        sample.sd=est.sd,
                        model=model
                        )

    class(result.list) <- "PBNLdiagnostic"

    if(save)
      save <- obj$arguments$save

   
    if(save)
      {
        name.save <- paste(obj$arguments$name.save,
                           ".diagnose", sep="")
        assign(name.save, result.list)
        save(result.list, list = name.save,
             file=paste( obj$arguments$save.directory,
               "/", name.save, ".rda",
               sep = ""))

        loglist <- c(list(true.par=true.par,
          autocor.max =autocor.max,
          default.thin=default.thin), list(...) )
    
        name.log=paste(name.save, ".log", sep="")
        assign(name.log, loglist  )
        save(loglist,list=name.log,
             file=paste( obj$arguments$save.directory,
               "/", name.log, ".rda",
               sep = ""))
      }
    
    return(result.list)
  }

##' @export print PBNLdiagnostic 
print.PBNLdiagnostic <- function(x,...)
  {
    cat("\n Model:", x$model, "\n", sep="\t")
    cat("\nPredictive angular density:",
        ifelse(is.matrix(x$predictive), "not shown", "not computed"),
        "\n", sep=" ")

    cat("Thinning:", x$thin, "\t from:",x$from, "\tto:", x$to,
        "\n", sep=" ")
    cat("Corresponding autocorrelations \n",
        diag(x$thinnedCorrel),"\n",  sep="\t")
    cat("Size of thinned sample used for the predictive:",
          x$size.subsample.predictive, "\n",  sep=" ")
    cat("Effective sizes (marginally):\n",
        x$effective.size,"\n",  sep="\t")
    cat("Heidelberger and Welches (stationarity part) p-values:\n",
        x$heidelTest[3,], "\n",  sep="\t")
    cat("Geweke p-values:\n",
        x$gewekeScore,"\n",  sep="\t")
    cat("posterior mean (transformed sample, defined on the real line)\n",
        x$linked.est.mean,   sep="\t")
    cat("\nposterior standard deviation  (idem)\n",
        x$linked.sample.sd,"\n",  sep="\t")
    
     cat("posterior mean (original parametrization)\n",
        x$est.mean,"\n",  sep="\t")
    cat("posterior standard deviation  (idem)\n",
        x$sample.sd,"\n",  sep="\t")
        
  }




## posterior.predictive.compare <-
##   function( pb.postsample = pb.postsample.Leeds,
##            nl.postsample=nl.postsample.Leeds,
##            dat=Leeds,
##            true.par=NULL,
##            npoints=60, eps=10^(-3),
##            equi=TRUE,
##            from =NULL,## 10e+3,
##            to=NULL,##20e+3,
##            autocor.max = 0.2,
##            default.thin=10,
##            ...
##              )
##   {
##     dev.new()
##     par(mai=c(0.7,0.1,0.7,0.1), mar=c(3.5,3,2,3))
##     if(equi)
##       {
##         Points = (apply((dat[, 1:2]), 1, transf.to.equi))

##       }
##     else
##       {
##         Points=t(dat)
##       }
    
##     pb.predictive <-
##       postsample.diagnose(model= "pairbeta", ##"nestlog"
##                           postsample = pb.postsample,
##                           dat=dat,
##                           true.par=true.par,
##                           npoints=npoints, eps=10^(-3),
##                           equi=equi,
##                           from =from, ## 10e+3,
##                           to=to,    ##20e+3,
##                           autocor.max = autocor.max,
##                           default.thin=default.thin,
##                           predictive=TRUE,
##                           plot=FALSE,
##                           save=FALSE, ...
##                           )
##     points(Points[1, ], Points[2, ], col=gray(0.5),
##            cex=0.8, pch=21, bg="white")

##     dev.new()
##     par(mai=c(0.7,0.1,0.7,0.1), mar=c(3.5,3,2,3))

##     nl.predictive <-
##       postsample.diagnose(model= "nestlog",
##                           postsample = nl.postsample,
##                           dat=dat,
##                           true.par=true.par,
##                           npoints=npoints, eps=10^(-3),
##                           equi=equi,
##                           from =from, ## 10e+3,
##                           to=to,    ##20e+3,
##                           autocor.max = autocor.max,
##                           default.thin=default.thin,
##                           predictive=TRUE,
##                           plot=FALSE,
##                           save=FALSE,...
##                           )

##      points(Points[1, ], Points[2, ],
##             col=gray(0.5), cex=0.8, pch=21, bg="white")
##   }
    
