##' Logarithm of the acceptance probability
##'
##' \code{lAccept.ratio} is a functional: \code{likelihood,proposal,prior} are user defined functions. Should not be called directly, but through the MCMC sampler \code{\link{posteriorMCMC}} generating the posterior.
##' @title Acceptance probability in the MCMC algorithm.
##' @param cur.par The current parameter in the Markov chain
##' @param prop.par The candidate parameter
##' @return The log-acceptance probability.
##' @inheritParams posteriorMCMC
##' @export
##' @keywords internal
lAccept.ratio <-
  function(cur.par,
           prop.par,
           llh.cur,
           lprior.cur,
           dat,
           likelihood,
           proposal,
           prior,
           Hpar, MCpar
           )
  {
    p <- ncol(dat)
    
    ## old.ll=likelihood(x=dat, par=cur.par,
    ##   log=TRUE, vectorial=FALSE)


    new.ll <- likelihood(x=dat, par=prop.par,
      log=TRUE, vectorial=FALSE)


    ## old.lprior = prior(type = "d", par = cur.par,
    ##    Hpar = Hpar, log = TRUE, dimData=p)

     new.lprior <- prior(type="d", par = prop.par,
       Hpar = Hpar, log = TRUE, dimData=p)

    proposal.oldToProp <- proposal (type = "d", 
           cur.par = cur.par,
           prop.par=prop.par,
           MCpar=MCpar, log=TRUE)


    proposal.propToOld <- proposal (type = "d", 
           cur.par = prop.par,
           prop.par=cur.par,
           MCpar=MCpar, log=TRUE)

    
    return(list(lrho=new.ll  - llh.cur +
                new.lprior-lprior.cur +
                proposal.propToOld - proposal.oldToProp,
                llh= new.ll,
                lprior=new.lprior) )
  }        

