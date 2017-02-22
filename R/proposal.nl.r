##' Density of the proposal distribution \code{q(cur.par,prop.par)} and random generator for MC MC algorithm in the NL3 model.
##'
##' The two components of proposal parameter
##' \code{(alpha*, beta12*, beta13*, beta23*)}  are generated independently, under a beta distribution with mode at the current parameter's value. 
##' 
##' Let  \eqn{\epsilon =} \code{MCpar$eps.recentre}. To generate \code{alpha*}, given the current state \code{alpha(t)},
##' let \eqn{m(t) = \epsilon /2 + (1-\epsilon) * \alpha(t)} be the mean
##' of the Beta proposal distribution  and \eqn{\lambda = 2/\epsilon}
##' (a scaling constant). Then
##' \deqn{
##' \alpha^* \sim  \textrm{Beta}(\lambda m(t), (1-\lambda) m(t))}{%
##' \alpha *  ~ Beta(\lambda m(t), \lambda(1-m(t))).}
##' The  \code{betaij*}'s are generated similarly.
##' @title NL3  model: proposal distribution.
##' @param type One of the character strings \code{"r"} or \code{"d"}.
##' @param cur.par Current state of the chain.
##' @param prop.par Candidate parameter.
##' @param MCpar A list made of a single element: MC MC parameter. Re-centering parameter for the proposal distribution.
##' @param log Logical. Only used when \code{type =="d"}. Should the result be returned on the log-scale ?
##' @return Either the (log-)density of the proposal parameter \code{prop.par}, given \code{cur.par} (if \code{type == "d"}), or a proposal  parameter (a vector), if \code{type =="r"}.
##' @export
proposal.nl <-
function(type = c("r", "d"),
         cur.par,
         prop.par,
         MCpar=get("nl.MCpar"), log=TRUE )
  {
    sd <- rep(MCpar$sd,length(cur.par))


    transfo <- function(x){logit(x)}
    invtransfo <- function(x){invlogit(x)}

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

                
    
##     mu.a=MCpar $ eps.recentre* rep(1/2,4) +
##       (1-MCpar $ eps.recentre)* cur.par
##     mu.b=1-mu.a
    
##     Nu.ab=2/ MCpar $eps.recentre
       
##     vect.a=mu.a*Nu.ab
##     vect.b=mu.b*Nu.ab
    
##     if(type=="r")
##       {
##         proposal=rbeta(4,shape1=vect.a,shape2=vect.b)
##         return(proposal)
##       }
##     if(type =="d")
##       {
##         vect.res = sapply(1:4,
##           function(j){dbeta(prop.par[j],
##                             shape1 = vect.a[j], shape2=vect.b[j],
##                             log=TRUE)} )
## #browser()
##         if(log)
##           return(sum(vect.res))

##         else
##           return(exp(sum(vect.res)))
##       }
    
  }                                     #nestlog.proposal()


