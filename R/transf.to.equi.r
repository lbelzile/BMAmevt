##' Switching coordinates system  between equilateral and right-angled
##' representation of the two-dimensional simplex.
##'
##' If  \code{transf.to.rect}, is called, \code{vect} must belongs to  the triangle \eqn{[(0,0), (\sqrt(2), 0), (\sqrt(2)/2,\sqrt(3/2) ) ]}{[(0,0), (\sqrt(2), 0), (\sqrt(2)/2,\sqrt(3/2) ) ]}
##' and the result lies in \eqn{([(0,0), (1,0), (0,1)]}{([(0,0), (1,0), (0,1)]}. \code{transf.to.equi} is the reciprocal.
##' @title Linear coordinate transformations
##' @param vect a bi-variate vector, giving the first two coordinates of the angular point to be transformed. 
##' @return The vector obtained by linear transformation. 
##' @author Anne Sabourin
##' @aliases transf.to.rect
##' @examples \dontrun{ transf.to.equi(c(sqrt(2)/2, sqrt(3/8) ) )}
##' @export transf.to.equi
transf.to.equi <-
  function(vect)
    {
      M=rbind(
        c(sqrt(2), sqrt(2)/2) ,
        c(0 ,      sqrt(3/2)) )
      return(as.vector(  M%*%t(t(vect)) ) )   
    }
