##' Density and generating function of the prior distribution.
##'
##' The parameters components are independent, log-normal.
##' @title Prior parameter distribution for the Pairwise Beta model
##' @param type One of the character strings \code{"r"}, \code{"d"}  
##' @param n The number of parameters to be generated. Only used
##' if \code{type == "r"}. 
##' @param par A vector with positive components:
##' The parameter where the density is to be taken.
##' Only used if \code{type=="d"}. In the Pairwise Beta model,
##' \code{par} is of length
##' \code{choose(p,2)+1}. The first element is the global dependence
##' parameter, the subsequent ones are the pairwise dependence
##' parameters, in lexicographic order (\emph{e.g.}
##' \eqn{\beta_{1,2}, \beta_{1,3}, \beta_{2,3}}.
##' @param Hpar list of Hyper-parameters : see \code{\link{pb.Hpar}} for a template.
##' @param log logical. Should the density be returned on the log scale ?
##' Only used if \code{type=="d"}
##' @param dimData The dimension of the sample space. (one more than the dimension of the simplex)
##' @return Either a matrix with \code{n} rows containing a random parameter sample generated under the prior (if type == "d"), or the (log)-density of the parameter \code{par}.
##' @examples \dontrun{prior.pb(type="r", n=5 ,Hpar=get("pb.Hpar"), dimData=3 ) }
##' \dontrun{prior.pb(type="d", par=rep(1,choose(4,2), Hpar=get("pb.Hpar"), dimData=4 ) }
##' @author Anne Sabourin
##' @export
prior.pb <-
  function(type=c("r","d"), n ,par, Hpar, log, dimData )
  {

    if(type =="r")
      {
        p <- dimData

        lengthPar <- choose(p,2)+1
       ## res <- matrix(0,ncol=lengthPar, nrow = n)
        alpha <- exp(rnorm(n, mean=Hpar$mean.alpha, sd=Hpar$sd.alpha))
        beta <-exp( matrix(rnorm(
                             (lengthPar-1)*n,
                             mean=Hpar$mean.beta, sd=Hpar$sd.beta),
                           ncol=3))
        res <- cbind(alpha,beta )
        return(res)
      }
    
    if(type =="d")
      {
        lpar <- log(par)
        ld1 <- dnorm(lpar[1], mean=Hpar$mean.alpha,
                     sd=Hpar$sd.alpha, log=TRUE)
        ld2 <- dnorm(lpar[-1], mean=Hpar$mean.beta,
                     sd=Hpar$sd.beta, log=TRUE)
##        browser()
        if(log)
          return(ld1 +sum(ld2))
        else
          return( exp(ld1+sum(ld2)))

        
        
       ## return(ifelse(log, sum(logdprior), exp(sum(logdprior))) )
  
      }
    stop("wrong 'type' argument")
  }
