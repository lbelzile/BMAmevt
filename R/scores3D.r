##' Computes the Kullback-Leibler divergence and the \eqn{L^2} distance between the "true" density (\code{true.dens}) and an estimated density (\code{est.dens}).
##'
##' The integration is made \emph{via} \code{\link{rect.integrate}}: The discretization grid corresponding to the two matrices must be constructed
##' with \code{discretize(npoints, eps, equi=FALSE)}. 
##' @title Logarithmic score and \eqn{L^2} distance  between two densities on the simplex (trivariate case).
##'  @param true.dens A \code{npoints*npoints} matrix: The reference density, typically the distribution from which data was simulated. Must be a valid \code{density} argument to be passed to \code{dgridplot}, with \code{equi=FALSE}. 
##' @param est.dens The estimated density: of the same type as \code{true.dens}. 
##' @param npoints Number of grid points  used to construct the density matrices (see \code{\link{discretize}}).
##' @param eps Minimum distance from a grid point to the simplex boundary (see \code{\link{discretize}}).
##' @return A list made of
##' \itemize{
##' \item \code{check.true}: The result of the rectangular integration of
##' \code{true.dens}. It should be equal to one. If not, re size the grid.
##' \item \code{check.true}:
##' Idem, replacing \code{true.dens} with \code{est.dens}.
##' \item \code{L2score}: The estimated \eqn{L^2} distance.
##' \item \code{KLscore}: The estimated Kullback-Leibler divergence between the two re-normalized densities, using \code{check.true} and \code{check.est} as normalizing constants (this ensures that the divergence is always positive). 
##' }
##' @examples
##' dens1=dpairbeta.grid(par=c(0.8,2,5,8),npoints=150,eps=1e-3,
##'                      equi=FALSE)
##' dens2=dnestlog.grid(par=c(0.5,0.8,0.4,0.6),npoints=150,eps=1e-3, equi=FALSE)
##' 
##' scores3D(true.dens=dens1,
##'   est.dens=dens2,
##'   npoints=150, eps=1e-4)
##' 
##' @export
scores3D=function(true.dens,
  est.dens,
  npoints, eps)
{
  check.true=rect.integrate(density=true.dens,
    npoints=npoints,eps=eps)
  check.est=rect.integrate(density=est.dens,
    npoints=npoints,eps=eps)

  L2.score= ( rect.integrate(density= (est.dens - true.dens)^2,
      npoints=npoints,eps=eps) )^(1/2)

  KL.score= rect.integrate(
    density= (log(true.dens+(true.dens==0))-log(check.true)-
    log(est.dens+ (est.dens==0))+
    log(check.est)  ) * true.dens/check.true,
    npoints=npoints,eps=eps)

  return(list(check.true=check.true,
              check.est=check.est,
              L2score=L2.score,
              KLscore=KL.score))
                            
}

## dens1=dpairbeta.grid(par=c(0.8,2,5,8),npoints=150,eps=1e-3, equi=FALSE)
## dens2=dnestlog.grid(par=c(0.5,0.8),npoints=150,eps=1e-3, equi=FALSE)

## scores3D(true.dens=dens1,
##   est.dens=dens2,
##   npoints=150, eps=1e-4)



## graphics.off()
## dens1=dpairbeta.grid(par=c(0.8,2,5,8),npoints=100,eps=1e-3, equi=T,displ=T)
## dev.new()
## dens2=dnestlog.grid(par=c(0.5,0.8),npoints=100,eps=1e-3, equi=T,displ=T)



## dens2=dpairbeta.grid(par=c(0.8,5,10,2),npoints=100,eps=1e-3, equi=T,displ=T)








