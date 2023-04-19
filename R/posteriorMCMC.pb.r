##' The functions generate  parameters samples approximating the posterior distribution in the PB model or the NL model.
##'
##' The two functions are wrappers simplifying the use of
##' \code{\link{posteriorMCMC}} for the two models implemented in this package.
##' @title MCMC posterior samplers for the pairwise beta and the negative logistic models.
##' @inheritParams posteriorMCMC
##' @param ...  Additional arguments to be passed to
##' \code{\link{posteriorMCMC}} instead of their
##' default values (must not contain any of
##' \code{"prior",
##' "likelihood", "proposal",
##' "name.model"} or \code{"class"}). 
##' @return an object with class attributes \code{"postsample"} and
##' \code{"PBNLpostsample"}: The posterior sample  and some statistics
##' as returned by function \code{\link{posteriorMCMC}}
##' @seealso \code{\link{posteriorMCMC}}
##' @note For the Leeds data set, and for simulated data sets with
##' similar features, setting \code{Nsim=50e+3} and \code{Nbin=15e+3}
##' is enough (possibly too much),
##' with respect to the Heidelberger and Welch tests implemented in
##' \code{\link[coda]{heidel.diag}}.
##' @examples
##' \dontrun{
##' data(Leeds)
##' data(pb.Hpar)
##' data(pb.MCpar)
##' data(nl.Hpar)
##' data(nl.MCpar)
##' pPB <- posteriorMCMC.pb(Nsim=5e+3, dat=Leeds, Hpar=pb.Hpar,
##' MCpar=pb.MCpar)
##' 
##' dim(pPB[1])
##' pPB[-(1:3)]
##'
##' pNL <- posteriorMCMC.nl(Nsim=5e+3, dat=Leeds, Hpar=nl.Hpar,
##' MCpar=nl.MCpar)
##'
##' dim(pNL[1])
##' pNL[-(1:3)]
##' }
##' @export
posteriorMCMC.pb <-
  function(Nsim,
           dat,
           Hpar,
           MCpar,
           ...
           )
  {
    
    postsample <-
      posteriorMCMC(Nsim=Nsim,
                    dat=dat,
                    prior=prior.pb,
                    proposal=proposal.pb,
                    likelihood=dpairbeta,
                    Hpar=Hpar,
                    MCpar=MCpar,
                    name.model="pairbeta",
                    class="PBNLpostsample",
                    ...)

    return(postsample)
  }


