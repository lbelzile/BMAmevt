##' @title Dirichlet distribution: random generator
##' @param n Number of draws
##' @param alpha Dirichlet parameter: a vector of positive number
##' @return A matrix with \code{n} rows and \code{length(alpha)} columns
##' @export
rdirichlet <- function (n=1, alpha) 
{
  if(is.vector(alpha))
    l <- length(alpha)
  else
    l <-  nrow(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- apply(x, 1, sum)
  return(x/ sm )
}



##' @rdname dpairbeta
##' @param n The number of points on the simplex to be generated.
##' @param dimData the dimension of the sample space, which is \eqn{1 + } the dimension of the simplex.
##' @keywords datagen distribution
##' @export
rpairbeta <-
function(n=1, dimData=3, par= c(1, rep(1,3)) )
  {
    p <- dimData
    if((choose(p,2)+1) != length(par) )
      { stop("wrong parameter length")}
    
    alpha <- par[1]
    beta <- par[-1]
    res=matrix(0,ncol=p,nrow=n)         #to store the result
    for(ll in 1:n)
      {
        ##choose uniformly a pair i,j:
        i_0=floor(runif(1,1,(p+1)))
        jStar=floor(runif(1,1,p))
        if(jStar<i_0){j_0=jStar}else{j_0=jStar+1}
        i=min(i_0,j_0);j=max(i_0,j_0)
        place=(i-1)*(p-(i)/2)+j-i ##place=\sum_{k=0}^{i_1}(p-k)  + (j-i)

        ##simulate r=wi+wj and theta=wi/r
        r_ij=rbeta(1,2*alpha+1,(p-2)*alpha) #2*alpha+1
        theta=rbeta(1,beta[place],beta[place])
    
        wi=theta*r_ij
        wj=(1-theta)*r_ij
    
                                        #simulate  the remaining components (p-2), uniformly on the simplex : "sum(w_i)=1-r"
        remains=(1-r_ij)*rdirichlet(1,rep(1,(p-2)))
    
                                        #put together r,theta and the remainings
        W0=append(as.vector(remains),wi,after=(i-1))
        W=append(W0,wj,after=(j-1))
    
                                        #keep the result...
        res[ll,]= W
      }
    return(res)
  
  
  }

## @importFrom  MCMCpack rdirichlet

