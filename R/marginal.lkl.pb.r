##'  Wrappers for \code{\link{marginal.lkl}}, in the specific cases of the PB and NL models,
##' with parameter \code{likelihood} set to \code{dpairbeta} or
##' \code{dnestlog}, and \code{prior} set to \code{prior.pb} or
##' \code{prior.nl}. See  \code{\link{MCpriorIntFun}} for more  details.
##' 
##' @title Marginal likelihoods of the PB and NL models.
##' @inheritParams marginal.lkl
##' @inheritParams MCpriorIntFun
##' @inheritParams posteriorMCMC
##' @return The list returned by
##' \code{\link{marginal.lkl}}, \emph{i.e.}, the one returned by \code{\link{MCpriorIntFun}} 
##' @export
##' @seealso \code{\link{marginal.lkl}}, \code{\link{MCpriorIntFun}} .
##' @examples
##' \dontrun{
##' 
##' marginal.lkl.pb(dat=Leeds ,
##'          Nsim=20e+3 ,
##'          displ=TRUE, Hpar = get("pb.Hpar") ,
##'           )
##' 
##' marginal.lkl.nl(dat=Leeds ,
##'          Nsim=10e+3 ,
##'          displ=TRUE, Hpar = get("nl.Hpar") ,
##'           )
##' }
marginal.lkl.pb <-
function(dat ,
         Nsim=10e+3 ,
         displ=TRUE, Hpar = get("pb.Hpar") ,
         Nsim.min=Nsim,
         precision=0,
         show.progress = floor(seq(1, Nsim, length.out = 20 ) )
         )
  {
    marginal.lkl(dat=dat,
                 likelihood=dpairbeta,
                 prior=prior.pb,
                 Nsim=Nsim,
                 displ=displ,
                 Hpar=Hpar,
                 Nsim.min=Nsim.min,
                 precision=precision,
                 show.progress=show.progress
                 )
  }    

