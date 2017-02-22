
## rdirichlet.dm <- function (n, alpha) 
## {
##     l <- length(alpha)
##     x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
##     sm <- apply(x, 1, sum)
##     return(x/ sm )
## }


##' @param n The number of angular points to be generated
##' @rdname ddirimix
##' @export
rdirimix <-
function(n=10,
         par=get("dm.expar.D3k3"),
         wei=par$wei,
         Mu=par$Mu ,
         lnu=par$lnu
         )
                 
 {
 if(is.vector(Mu))
    {
     Mu=matrix(Mu,ncol=1)
    }
 p=nrow(Mu)
 k=length(wei)
## dat=matrix(0,ncol=p, nrow=n)
 ## for(i in 1:n)
 ##    {
 ##     u=runif(1)
 ##     cum_wei=cumsum(wei)
 ##     lowers=cum_wei[cum_wei<u]
 ##     m=length(lowers)+1
 ##     dat[i,]=rdirichlet(1,exp(lnu[m])*Mu[,m] )
 ##    }
 u <- runif(n)
 cum_wei <- cumsum(wei)
 ms <- sapply(1:n, function(i){
   length(which(cum_wei<u[i])) + 1})
 if(n>1){
#   matpars <- matrix(Mu[,ms], ncol=length(ms))%*%diag(exp(lnu[ms]))
      matpars <-Mu[,ms]%*%diag(exp(lnu[ms]))
 }
 else matpars <- matrix(Mu[,ms]*exp(lnu[ms]))
 dat <- rdirichlet(n=n, alpha=matpars)
   return(dat )
 
 } 


