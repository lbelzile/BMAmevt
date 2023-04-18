##' Bijective Transformation from  \eqn{(0,1)} to the real line, defined
##' by
##' \eqn{logit(p) = log( p / (1-p) )}.
##'
##' @title Logit transformation
##' @return  A real number
##' @param p A real number in \eqn{[0,1]}
##' @export
logit <- function (p) 
{
    if ( any(p < 0)  ||     any(p > 1) ) 
        stop("p must be in [0,1].")
    return( log(p)-log( 1 - p) )
}


##' Inverse transformation of the \code{logit} function.
##'
##' @title Inverse logit transformation
##' @return A number between \eqn{0} and \eqn{1}.  
##' @param x A real number
##' @export
invlogit <-  function (x) 
{
  return( 1/(1 + exp(-x)) )
}





##' Density and generating function of the prior distribution.
##'
##' The four  parameters are independent, the logit-transformed parameters follow a normal distribution.
##' @title Prior parameter distribution for the NL  model
##' @param type One of the character strings \code{"r"}, \code{"d"}  
##' @param n The number of parameters to be generated. Only used
##' if \code{type == "r"}. 
##' @param par A vector of length four, with component comprised
##' between \eqn{0} and \eqn{1} (both end points excluded for the first
##' element and \eqn{1} included for the others):
##' The parameter where
##' the density is to be taken.
##' Only used if \code{type=="d"}.
##'
##' In the NL model,
##' \code{par} is of length \eqn{4}.
##' The first element is the global dependence
##' parameter, the others  are  partial dependence parameter between pairs (12), (13), (23) respectively.
##'
##' In the NL model,
##' \code{par} is of length \eqn{4}.
##' The first element has the same interpretation as in the NL model, the subsequent ones are dependence parameters between 
##' @param Hpar list of Hyper-parameters : see \code{\link{nl.Hpar}} for a template.
##' @param log logical. Should the density be returned on the log scale ?
##' Only used if \code{type=="d"}
##' @param dimData The dimension of the sample space, equal to \eqn{3}.
##' Only for compatibility with \emph{e.g.} \code{\link{posteriorMCMC}}.
##' @return Either a matrix with \code{n} rows containing a random parameter sample generated under the prior (if type == "d"), or the (log)-density of the parameter \code{par}.
##' @examples \dontrun{prior.nl(type="r", n=5 ,Hpar=get("nl.Hpar")) }
##'
##' \dontrun{prior.trinl(type="r", n=5 ,Hpar=get("nl.Hpar")) }
##' \dontrun{prior.pb(type="d", par=rep(0.5,2), Hpar=get("nl.Hpar")) }
##' @author Anne Sabourin
##' @export
 prior.nl <- function(type=c("r","d"), n ,par, Hpar, log, dimData=3)
  {
    if(type =="r")
      {
        stopifnot(dimData == 3)
        alpha <- invlogit(rnorm(n, mean=Hpar$mean.alpha, sd=Hpar$sd.alpha) )
        beta <- invlogit(matrix(rnorm(3*n, mean=Hpar$mean.beta,
                                      sd=Hpar$sd.beta), ncol=3))
        res <- c(alpha,beta )
        names(res) <-  c("alpha","beta12","beta13", "beta23")
        return(res)
      }

    if(type=="d")
      {
        lpar <- logit(par)
        
        ld1 <- dnorm(lpar[1], mean=Hpar$mean.alpha,
                     sd=Hpar$sd.alpha, log=TRUE)
        ld2 <- dnorm(lpar[-1], mean=Hpar$mean.beta,
                     sd=Hpar$sd.beta, log=TRUE)
        
        if(log)
          return(ld1 +sum(ld2))
        else
          return( exp(ld1+sum(ld2)))

      }
  }
