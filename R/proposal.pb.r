##' Density of the proposal distribution \code{q(cur.par,prop.par)} and random generator for MCMC algorithm in the PB model.
##'
##' The components \code{prop.par[i]} of the proposal parameter are generated independently, from the lognormal distribution:
##'
##' \code{prop.par = rlnorm(length(cur.par), meanlog=log(cur.par),
##' sdlog=rep(MCpar$sdlog,length(cur.par)))}
##' @title PB model: proposal distribution 
##' @param type One of the character strings \code{"r"}, \code{"d"} 
##' @param cur.par Current state of the chain
##' @param prop.par Candidate parameter
##' @param MCpar A list made of a single element: MCMC parameter for the standard deviation of the log-normal proposition, on the log scale. See \code{\link{pb.MCpar}} for the default value
##' @param log Logical. Only used when \code{type =="d"}. Should the result be returned on the log-scale ?
##' @return Either the (log-)density of the proposal \code{prop.par}, given \code{cur.par} (if \code{type == "d"}), or a proposal  parameter (a vector), if \code{type =="r"}.
##' @examples \dontrun{ proposal.pb(type = "r",
##' cur.par = rep(1,4), MCpar=get("pb.MCpar"))
##' }
##' \dontrun{ proposal.pb(type = "d", cur.par = rep(1,4),
##' prop.par=rep(1.5,4), MCpar=get("pb.MCpar"))
##' }
##' @export
proposal.pb <-
  function(type = c("r","d"),
           cur.par,
           prop.par,
           MCpar, log=TRUE)
  {
    sd <- rep(MCpar$sd,length(cur.par))


    transfo <- function(x){log(x)}
    invtransfo <- function(x){exp(x)}

    mean <- transfo(cur.par)
    
    if(type =="r")
      {
        return(invtransfo(rnorm(length(cur.par), mean=mean,sd=sd)))
      }
    if(type == "d")
      {
        vect.res=sapply(1:length(prop.par),
          function(j){dnorm(transfo(prop.par[j]), mean=mean[j],
                             sd=sd[j],
                             log=TRUE)})

        
          return(ifelse(log,sum(vect.res),exp(sum(vect.res)) ))
       
      }
    
    
    stop("wrong type specification")
  }
    ## new.lAlpha=rnorm(1,mean=log(par[1]),
   ##   sd=sqrt(MCpar$sig_jumpA))   

   ## new.lBeta=rnorm(
   ##   length(par)-1,
   ##   mean=log(par[-1]),
   ##   sd=sqrt(MCpar$sig_jumpB)
   ##   )

   ## return(exp(c(new.lALpha, new.lBeta)))
   
   
 ## } 

