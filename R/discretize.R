##' Builds a discretization grid covering the two-dimensional unit simplex, with specified number of points and minimal distance from the boundary.
##'
##' The \code{npoints*npoints} grid  covers either
##' the equilateral representation of
##' the simplex, or the  right angled one.
##' In any case, the grid is 
##' \emph{rectangular}: some nodes lie outside the triangle.
##' Density computations on such a grid should handle the case when
##' the point passed as argument is outside the simplex (typically,
##' the function should return zero in such a case). 
##' @title Discretization grid builder.
##' @param npoints The number of grid nodes on the squared grid containing the desired triangle. 
##' @param eps Positive number:  minimum
##' distance from any node inside the simplex to  the simplex boundary
##' @param equi logical. Is the simplex represented as an equilateral triangle (if \code{TRUE}) or a right triangle (if \code{FALSE}) ?
##' @note In case \code{equi==TRUE}, \code{epsilon} is the  minimum
##' distance from any node inside the simplex to  the simplex boundary,
##' \emph{after transformation} to the right-angled representation.
##' @return A list containing two elements: \code{X} and \code{Y}, vectors of size \code{npoints}, the Cartesian coordinates of the grid nodes.
##' @export
discretize <- function(npoints=40,eps=1e-3,equi=FALSE)
  {
   if(!equi)
        {
        X=seq(eps,1-2*eps,length.out=npoints)
        Y=X
        return(list(X=X,Y=Y))
        }
   else
        {
         ## X=sqrt(2)*seq(eps,1-eps,length.out=npoints)
        ###  Y=sqrt(3/2)*seq(eps/3,1-2*eps/3,length.out=npoints)
          ##Y=seq(sqrt(2/3)*eps, sqrt(3/2) -sqrt(2)*eps, length.out=npoints)
          X=seq(eps*3/sqrt(2), sqrt(2)-eps*3/sqrt(2), length.out=npoints)
          Y=seq(sqrt(3/2)*eps,sqrt(3/2)-eps*3/sqrt(2)  ,length.out=npoints)
          return(list(X=X,Y=Y))
        } 
     
  }

