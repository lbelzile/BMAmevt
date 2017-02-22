##' Plots a univariate Dirichlet mixture (in other words, a Beta mixture) angular density for extreme bi-variate data.
##'
##' @title Univariate projection or marginalization of a  Dirichlet mixture density on  on \code{[0,1]}
##' @inheritParams ddirimix
##' @param coord A vector of size 2:
##' the indices of the coordinates upon which the marginalization or projection is to be done if the dimension of the sample space is greater than two.
##' @param marginal logical. If \code{TRUE}, the angular density corresponds to the marginal intensity measure of the extreme Poisson process, over coordinates \code{coord}. Otherwise, it is only the projection of the full dimensional angular measure (hence the moments constraints is not satisfied anymore). 
##' @param npoints number of points on the 1D  discretization grid.
##' @param eps the minimum value ( = 1- the maximum value) of the grid points.
##' @param invisible Logical: should the result be returned as invisible ?
##' @param displ Logical: should a plot be issued ?
##' @param add Logical: should the density be added to the currently active plot ?
##' @param ... Additional arguments to be passed to \code{plot}
##' @return The discretized density on \code{[eps, 1-eps]} (included in [0,1])
##' @export
ddirimix.grid1D <-
  function(par= get("dm.expar.D2k4"),
           wei=par$wei,
           Mu= par$Mu,
           lnu=par$lnu,
           npoints=30,
           eps=10^(-3),
           coord=c(1,2), marginal=TRUE,
           invisible=TRUE,
           displ=TRUE, add=FALSE,
           ...
           )
  {
    dim <- nrow(Mu)
    if(dim>2)
      {
        if(length(coord) !=2)
          {
            warning("coord is not of length 2")
            coord <- coord[1:2]
          }
        MMu <- matrix(Mu[coord,], nrow=length(coord))
        if(marginal)
          {
            multconst <- apply(MMu,2,sum)
            if(length(multconst)>1)
              Mu <- MMu %*% diag(1/multconst)
            else
              Mu <- MMu/multconst
            
            lnu <- lnu + log(multconst)
            wei <- dim/2*wei*multconst
          }
        else
          {
            Mu <- MMu
          }
      }

    if(eps<0 | eps>=1)
      warning("eps should be in [0,1)")
    X.grid <- seq(eps,1-eps, length.out=npoints)
    nu <- exp(lnu)
    k <- length(lnu)
    Density <- double(npoints)

    C.out <- .C("ddirimix_grid1D", as.double(X.grid), 
                as.integer(npoints), as.double(Mu),
                as.integer(k), as.double(wei), as.double(nu),
                result=Density) 

                                 
    if(displ){
      if(!add)
        plot(X.grid, C.out$result, type="l", ...)
    
    else
      lines(X.grid, C.out$result, type="l", ...)
    }

    if(!invisible)
      return( C.out$result)
      
    else
      invisible( C.out$result)
  }

 
