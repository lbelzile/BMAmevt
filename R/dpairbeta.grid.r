##' The two functions compute respectively  the NL and PB spectral
##' densities,  in the three-dimensional case,  on a discretization grid.
##' A plot is issued (optional).
##'
##' @title PB and NL spectral  densities on the two-dimensional simplex
##' @inheritParams discretize
##' @inheritParams dpairbeta
##' @param displ logical. Should a plot be produced ?
##' @param invisible logical. If \code{TRUE}, the result is returned as \code{invisible}.
##' @param ... Additional arguments to be passed to \code{\link{dgridplot}}
##' @return A \code{npoints*npoints} matrix containing the
##' considered density's values on the grid.
##' The row (resp. column) indices increase
##' with  the first (resp. second) coordinate on the simplex.
##' @note If \code{equi==TRUE}, the density is relative to the Hausdorff
##' measure on the simplex itself: the values obtained with
##' \code{equi = FALSE} are thus divided by
##' \eqn{\sqrt 3}.
##' @export
##' @examples
##' 
##' dpairbeta.grid(par=c( 0.8, 8, 5, 2),
##' npoints=70, eps = 1e-3, equi = TRUE, displ = TRUE, invisible=TRUE)
##' 
##' ##  or ...
##' 
##' Dens <- dpairbeta.grid(par=c(0.8, 8, 5, 2),
##' npoints=70, eps = 1e-3, equi = TRUE, displ = FALSE)
##' Grid=discretize(npoints=70,eps=1e-3,equi=TRUE)
##' dev.new()
##' image(Grid$X, Grid$Y, Dens)
##' contour(Grid$X, Grid$Y, Dens, add=TRUE)
##' add.frame(equi=TRUE, npoints=70, axes=FALSE)
##' 
##' 
dpairbeta.grid <- function(par,
  npoints=50, eps = 1e-3, equi = TRUE, displ = TRUE, invisible=TRUE,
   ...) 
  {
    discr <- discretize(npoints=npoints,eps=eps,equi=equi)

    X.grid <- discr $X
    Y.grid <- discr $ Y

    if(length(par) != 4 | any(par<0))
      stop("misspecified parameter")
    
    C.out <- .C("d_pairbeta_grid", as.double(X.grid),
      as.double(Y.grid),
      as.integer(length(X.grid)),
      as.double(par[1]),
      as.double(par[-1]), as.integer(equi), 
     result =  as.double(rep(0,npoints*npoints)) )

    density <- matrix(C.out $ result, ncol = npoints, byrow=F)

    mult.dens <- 1/sqrt(3)*equi + 1*(!equi)
    density <- mult.dens*density
    if(displ)
      {
        dgridplot(density=density,
##                   npoints=npoints,
                   eps=eps,
                   equi=equi,
                  ...)

      }
    if(invisible)
      {invisible(density)}
    else
      return(density)

  }


