##' @export
##' @rdname posteriorMCMC.pb
posteriorMCMC.nl <-
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
                    prior=prior.nl,
                    proposal=proposal.nl,
                    likelihood=dnestlog,
                    Hpar=Hpar,
                    MCpar=MCpar,
                    name.model="nestlog",
                    class="PBNLpostsample",
                    ...)
    
    ## oldclasses <- class(postsample)
    ## class(postsample) <- c("PBNLpostsample", oldclasses)
    return(postsample)
  }


## posteriorMCMC.nl <-
##   function(Nsim=500,Nbin=200,
##            dat=get("Leeds"),
##            par.start=NULL, 
##            Hpar=get("nl.Hpar"),
##            MCpar = get("nl.MCpar"),
##            show.progress = floor(seq(1, Nsim, length.out = 20 ) ),
##            seed=1, kind ="Mersenne-Twister",
##            name.save=NULL,
##            save.directory="~", name.dat=NULL
##            )
##   {
##     post <-
##       posteriorMCMC(Nsim=Nsim,Nbin=Nbin,
##                     dat=dat,
##                     par.start=par.start,
##                     prior=prior.nl,
##                     proposal=proposal.nl,
##                     likelihood=dnestlog,
##                     Hpar=Hpar,
##                     MCpar=MCpar,
##                     show.progress = show.progress,
##                     seed=seed, kind=kind, name.save=name.save,
##                     save.directory=save.directory, name.dat=name.dat,
##                     name.model="nestlog")
    
##     return(post)
##   }
