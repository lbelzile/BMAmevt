##' The integral is approximated by a rectangular method, using the values stored in matrix \code{density}. 
##' 
##' Integration is made with respect to the Lebesgue measure on the projection of the simplex onto the plane \eqn{(x,y): x > 0, y > 0, x+y < 1}.
##' It is assumed that \code{density} has been constructed on a
##' grid obtained \emph{via} function  \code{\link{discretize}},
##' with  argument \code{equi} set to \code{FALSE} and \code{npoints}
##' and \code{eps} equal to those passed to \code{rect.integrate}.
##' @title Density integration on the two-dimensional simplex 
##' @inheritParams dgridplot
##' @inheritParams discretize
##' @return  The value of the estimated integral of \code{density}.
##' @examples
##' wrapper <- function(x, y, my.fun,...)
##'   {
##'     sapply(seq_along(x), FUN = function(i) my.fun(x[i], y[i],...))
##'   }
##' 
##' grid <- discretize(npoints=40,eps=1e-3,equi=FALSE)
##' 
##' Density <- outer(grid$X,grid$Y,FUN=wrapper,
##'                                  my.fun=function(x,y){10*((x/2)^2+y^2)})
##' 
##' rect.integrate(Density,npoints=40,eps=1e-3)
##'
##' @export
rect.integrate <-
  function(density,npoints,eps)
  {
    
    wrapper <- function(x, y, my.fun,...)
      {
        sapply(seq_along(x), FUN = function(i) my.fun(x[i], y[i],...))
      }
    
    grid <- discretize(npoints=npoints,eps=eps,equi=FALSE)
    mask <- outer(grid$X,grid$Y,FUN=wrapper,
                  my.fun=function(x,y){(x+y)<1})

    maskDensity=mask*density
    
    sum(maskDensity[2:npoints,2:npoints]) * ((1-3*eps)/(npoints-1))^2 +
      sum(maskDensity[1,(2:npoints)])* eps * (1-3*eps)/(npoints-1) +
        sum(maskDensity[(2:npoints),1]) * eps * (1-3*eps)/(npoints-1) +
          maskDensity[1,1]* eps^2  
  } 
          

