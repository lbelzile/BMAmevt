##' @export
##' @rdname marginal.lkl.pb
marginal.lkl.nl <-
function(dat,
         Nsim=10e+3,
         displ=TRUE, Hpar = get("nl.Hpar"),
         Nsim.min=Nsim,
         precision=0,
         show.progress = floor(seq(1, Nsim, length.out = 20 ) )
         )
  {
    marginal.lkl(dat=dat,
                 likelihood=dnestlog,
                 prior=prior.nl,
                 Nsim=Nsim,
                 displ=displ,
                 Hpar=Hpar,
                 Nsim.min=Nsim.min,
                 precision=precision,
                 show.progress=show.progress
                 )

  }

