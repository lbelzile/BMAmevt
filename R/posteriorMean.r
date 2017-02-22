##' Computes an approximation of the posterior mean of a parameter functional, based on a  posterior parameters sample.
##'
##' Only a sub-sample is used: one out of \code{thin} parameters is used
##' (thinning). Further, only the parameters produced between time
##' \code{from} and time \code{to} (included) are kept. 
##' @title Posterior predictive density on the simplex, for three-dimensional extreme value  models.
##' @param post.sample A posterior sample as returned by \code{\link{posteriorMCMC}}
##' @param FUN a parameter functional returning a vector.
##' @param ... Additional parameters to be passed to \code{FUN}.
##' @param from Integer or \code{NULL}. If \code{NULL}, the default value is used. Otherwise,  should be greater than \code{post.sample$Nbin}. Indicates the  index where the averaging process should start. Default to \code{post.sample$Nbin +1}
##' @param to Integer or \code{NULL}. If \code{NULL}, the default
##' value is used. Otherwise, must be lower than \code{Nsim+1}.
##' Indicates  where the averaging process should stop.
##' Default to \code{post.sample$Nsim}.
##' @param thin Thinning interval.
##' @param displ logical. Should a plot be produced ?
##' @return A list made of \describe{
##' \item{values}{A matrix : each column is the result of \code{FUN} applied to a parameter from the posterior sample.}
##' \item{est.mean}{The posterior mean}
##' \item{est.sd}{The posterior standard deviation }
##' }
##' @seealso \code{\link{posteriorMCMC}}.
##' @export 
posteriorMean <-
  function(post.sample , 
           FUN=function(par,...){par},
           from = NULL, to = NULL, thin=50,
           displ=TRUE,
           ...
           )

  {
    if(is.null(from))
      from <- post.sample$Nbin +1
    if(is.null(to))
      {
        to <- post.sample$Nsim
      }
    if((from<= post.sample$Nbin) | (to>post.sample$Nsim) )
      {stop(" argument \"from\" or \"to\" out of range")}
    
    time.idx <- (from-post.sample$Nbin):(to-post.sample$Nbin)
    kept.idx <- time.idx[ (time.idx %% thin == 0) ]
#    mat.postSample <- post.sample $ stored.vals[kept.idx,] 

    values <- sapply(kept.idx,
                     function(i){as.vector(FUN(par=post.sample $
                                     stored.vals[i,], ...) )},
                     simplify=TRUE
                     )
#    mean.res=matrix(0,ncol=npoints, nrow=npoints)

    if(is.vector(values))
      values <- matrix(values,nrow=1)
    
    est.mean <- apply(values, 1, mean)
    est.sd <- apply(values,1,sd)
    if(displ)
      {
        N <- length(kept.idx)
        shuffle <- sample(N,N,replace=F)
        for(i in 1:nrow(values))
          {
            
            dev.new()
            ylim <- range(values[i,])
            cummean <- cumsum(values[i,shuffle])/(1:N)
            plot(kept.idx, cummean, type="l" )
            esterr <- sqrt(cumsum( (values[i,shuffle]- cummean)^2 ) )/1:N
            lines( cummean+2*esterr, col=gray(0.5) )
            lines( cummean-2*esterr, col=gray(0.5) )
            
          }

      }

    return(list(values=values, emp.mean=est.mean, emp.sd=est.sd))
  }
