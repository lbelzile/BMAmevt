##' Random variable generator 
##'
##' An alpha-stable random variable \eqn{S} with index \eqn{\alpha} is
##' defined by its  Laplace transform
##' \eqn{E(exp(tS))) = exp(-t^\alpha)}. The algorithm used here is directly derived from Stephenson (2003).
##' @title Positive alpha-stable distribution.
##' @param alpha The parameter of the alpha-stable random variable 
##' @return A realization of the alpha-stable random variable.
##' @references STEPHENSON, A. (2003). Simulating multivariate extreme value distributions of logistic type. \emph{Extremes 6, 49-59}.
##' @export
rstable.posit <-
function(alpha=0.5)
  {
   if(alpha==1)
    {
    return(1)
    }
   else
    {
     U=runif(1,0,pi)
     W=rexp(1)
     
     A1=sin((1-alpha)*U)/W
     A2=sin(alpha*U)
     A3=(sin(U))^(1/alpha)
     return(A1^((1-alpha)/alpha) *A2 / A3 )
    } 
  }

test.rstable.posit <-
  function()
  {
    X=c()
    for(i in 1:5000)
      {
        X=c(X,rstable.posit(alpha=2/3))
      }
    hist(X,probability=TRUE) ##,xlim=c(0,10), break = 100)
  }


##test.rstable.posit()
