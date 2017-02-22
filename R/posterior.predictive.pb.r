##'  Wrappers for \code{\link{posterior.predictive3D}} in the PB and NL models.
##'
##' The posterior predictive density is approximated by averaging the densities corresponding to the parameters stored in \code{post.sample}. See
##' \code{\link{posterior.predictive3D}} for details.
##' @title Posterior predictive densities in the three dimensional
##' PB, NL and NL3  models
##' @inheritParams posterior.predictive3D
##' @return A \code{npoints*npoints} matrix: the posterior predictive density.
##' @seealso \code{\link{posterior.predictive3D}},  \code{\link{posteriorMCMC.pb}}.
##' @export
posterior.predictive.pb= function(
  post.sample , 
  from = post.sample$Nbin+1, to = post.sample$Nsim, thin = 50,  
  npoints=40,eps=10^(-3), equi = T, displ=T,
  ...
 )

  {

    predictive <-
      posterior.predictive3D(post.sample = post.sample , 
                             densityGrid =dpairbeta.grid,
                             from = from, to = to, thin = thin,  
                             npoints=npoints, eps=eps, equi = equi,
                             displ=displ,
                             ...
                             )

    return(predictive)
  }

    
