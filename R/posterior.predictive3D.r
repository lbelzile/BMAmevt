##' Computes an approximation of the predictive density based on a  posterior parameters sample. Only allowed in the three-dimensional case.
##'
##' The posterior predictive density is approximated by averaging the
##' densities produced by the function
##' \code{densityGrid(par, npoints, eps, equi, displ,invisible, ...)} for
##' \code{par} in a  subset of  the parameters sample stored in
##' \code{post.sample}. The arguments of \code{densityGrid} must be
##' \itemize{
##' \item \code{par}:  A vector containing the parameters.
##' \item \code{npoints, eps, equi}: Discretization parameters
##' to be passed to  \code{\link{dgridplot}}.
##' \item \code{displ}:  logical. Should a plot be produced ?
##' \item \code{invisible}: logical. Should the result be returned as \code{invisible} ?
##' \item \code{...} additional arguments to be passed to
##' \code{\link{dgridplot}}
##' }
##' Only a sub-sample is used: one out of \code{thin} parameters is used
##' (thinning). Further, only the parameters produced between time
##' \code{from} and time \code{to} (included) are kept. 
##' @title Posterior predictive density on the simplex, for three-dimensional extreme value  models.
##' @param post.sample A posterior sample as returned by \code{\link{posteriorMCMC}}
##' @param densityGrid A function returning a \code{npoints*npoints}
##' matrix, representing a discretized version of the spectral density
##' on the two dimensional simplex.
##' The function should be compatible with \code{\link{dgridplot}}.
##' In particular, it must  use \code{\link{discretize}} to produce
##' the discretization grid. It must be of type  \cr
##' \code{function(par, npoints, eps, equi, displ,invisible,
##' ... )}.
##' See \bold{Details} below.
##' @param from Integer or \code{NULL}. If \code{NULL}, the default value is used. Otherwise,  should be greater than \code{post.sample$Nbin}. Indicates the  index where the averaging process should start. Default to \code{post.sample$Nbin +1}
##' @param to Integer or \code{NULL}. If \code{NULL}, the default
##' value is used. Otherwise, must be lower than \code{Nsim+1}.
##' Indicates  where the averaging process should stop.
##' Default to \code{post.sample$Nsim}.
##' @param thin Thinning interval.
##' @inheritParams dgridplot
##' @inheritParams discretize
##' @param displ logical. Should a plot be produced ?
##' @note The computational burden may be high: it is proportional to
##' \code{npoints^2}. Therefore, the function assigned to
##' \code{densityGridplot} should be
##' optimized, typically by calling \code{.C} with an internal,
##' user defined \code{C} function. 
##' @return A \code{npoints*npoints} matrix: the posterior predictive density.
##' @seealso \code{\link{dgridplot}},  \code{\link{posteriorMCMC}}.
##' @author Anne Sabourin
##' @include dgridplot.R
##' @export 
posterior.predictive3D <-
  function(post.sample , 
           densityGrid,
           from = post.sample$Nbin+1, to = post.sample$Nsim, thin = 40,  
           npoints=40,eps=10^(-3), equi = T, displ=T,
           ...
           )

  {
    
    Nsim = post.sample$Nsim
    Nbin = post.sample $ Nbin
    if(is.null(from))
      {
        from=Nbin+1
      }
    if(is.null(to))
      {
        to=Nsim
      }
    if((from<=Nbin) | (to>Nsim) )
      {stop(" argument \"from\" or \"to\" out of range")} 
    time.idx = (from-Nbin):(to-Nbin)
    kept.idx = time.idx[ (time.idx %% thin == 0) ]
    mat.postSample = post.sample $ stored.vals[kept.idx,] 
    mean.res=matrix(0,ncol=npoints, nrow=npoints)
  
    for(i in 1:length(kept.idx))
      {
        mean.res=mean.res+
          densityGrid(par=mat.postSample[i,],
                      npoints=npoints, eps=eps,equi=equi,
                      displ=FALSE,invisible=FALSE)
      }
  
    mean.res=mean.res/length(kept.idx)

    if(displ)
      {dgridplot(density = mean.res, ##pred.density,
#                 npoints = npoints, 
                 eps = eps,
                 equi = equi,
                 ... )
     }

    return(mean.res)
  }
