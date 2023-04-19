##' Generates a posterior parameters sample, and computes the posterior  mean and component-wise variance on-line.
##'
##' @title MCMC sampler for parametric spectral measures
##' @param Nsim Total number of iterations to perform.
##' @param Nbin Length of the burn-in period.
##' @param par.start Starting point for the MCMC sampler.
##' @param dat An angular data set, \emph{e.g.}, constructed by
##' \code{\link{cons.angular.dat}}: A matrix which rows are the Cartesian coordinates of points on the unit simplex (summing to one). 
##' @param likelihood The likelihood function.
##' Should be of type\cr
##' \code{function(x, par, log, vectorial)}, where \code{log} and
##' \code{vectorial} are logical flags indicating respectively if
##' the result is to be  returned on the log-scale, and if the
##' value is   a vector of length \code{nrow(x)} or a single number
##' (the likelihood, or the log-likelihood, for the data set \code{x}).
##' See \code{\link{dpairbeta}} or \code{\link{dnestlog}}
##' for templates.
##' @param proposal The proposal function: of type \cr
##' \code{function(type = c("r","d"),  
##' cur.par, prop.par, MCpar, log)
##' }.
##' Should
##' return the (logarithm of) the proposal density for the move
##' \code{cur.par --> prop.par} if \code{type == "d"}. If
##' \code{type =="r"}, \code{proposal} must return a candidate
##' parameter, depending on \code{cur.par}, as a vector.
##' See \code{\link{proposal.pb}} or \code{\link{proposal.nl}}
##' for templates.
##' @param prior The prior distribution: of type \cr
##' \code{function(type=c("r","d"), 
##' n ,par, Hpar, log, dimData
##' )},
##' where \code{dimData} is the dimension of the sample
##' space (\emph{e.g.}, for
##' the two-dimensional simplex (triangle), \code{dimData=3}.
##' Should return either a matrix with \code{n} rows containing a
##' random parameter sample generated under the prior
##' (if \code{type == "d"}), or the density of the
##' parameter \code{par} (the logarithm of the density if
##' \code{log==TRUE}.
##' See \code{\link{prior.pb}} and \code{\link{prior.nl}} for templates.
##' @param Hpar A list containing  Hyper-parameters to be passed to
##' \code{prior}.
##' @param MCpar A list containing  MCMC tuning parameters to be
##' passed to \code{proposal}.
##' @param show.progress An vector of integers containing the times
##' (iteration numbers) at  which a message showing progression
##'  will be printed on the standard output.
##' @param seed The seed to be set \emph{via}
##' \code{\link[base]{set.seed}}. 
##' @param kind The kind of random numbers generator. Default to
##' "Mersenne-Twister". See \code{\link[base]{set.seed}} for details.
##' @param save Logical. Should the result be saved ?
##' @param class Optional character string: additional class attribute to be assigned to the result. A predefined class \code{"PBNLpostsample"} exists, for which a method performing  convergence diagnostics is defined (see \code{\link{diagnose}} )
##' @param save.directory A character string giving the directory where the result is to be saved (without trailing slash). 
##' @param name.save A character string giving the name under which
##' the result is to be saved. If \code{NULL} (default),
##' nothing is saved. Otherwise, the result is saved in file \cr
##' \code{paste(save.directory,"/",
##' name.save,".rda",sep="")}.
##' A "log" list  is also saved,  named \cr
##' \code{paste(name.save, ".log", sep="")},  in file \cr
##' \code{paste(save.directory,"/", name.log,".rda",sep="")}.
##' @param name.dat A character string naming  the data set used for inference. Default to \code{""}. 
##' @param name.model A character string naming the model. Default to \code{""}.
##' @return A list made of
##' \itemize{
##' \item \code{stored.vals}: A \code{(Nsim-Nbin)*d} matrix, where
##' \code{d}
##' is the dimension of the parameter space.
##' \item \code{llh} A vector of size \code{(Nsim-Nbin)} containing the loglikelihoods evaluated at each parameter of the posterior sample.
##' \item \code{lprior} A vector of size \code{(Nsim-Nbin)} containing the logarithm of the prior densities evaluated at each parameter of the posterior sample.
##' \item \code{elapsed}: The time elapsed, as given by
##' \code{proc.time} between the start and the end of the run.
##' \item \code{Nsim}: The same as the passed argument
##' \item \code{Nbin}: idem.
##' \item\code{n.accept}: The total number of accepted proposals.
##' \item \code{n.accept.kept}: The number of accepted proposals after the burn-in period.
##' \item \code{emp.mean} The estimated posterior parameters mean
##' \item \code{emp.sd} The empirical posterior sample  standard deviation.}
##' @export
##' @seealso \code{\link{posteriorMCMC.pb}},
##' \code{\link{posteriorMCMC.pb}} for specific uses
##' in the PB and the NL models.
##' @examples
##' data(Leeds)
##' data(pb.Hpar)
##' data(pb.MCpar)
##' postsample1 <- posteriorMCMC(Nsim=1e+3,Nbin=500,
##'          dat= Leeds,
##'          prior = prior.pb,
##'          proposal = proposal.pb,
##'          likelihood = dpairbeta,
##'          Hpar=pb.Hpar,
##'          MCpar=pb.MCpar)
##' 
##' dim(postsample1[[1]])
##' postsample1[-1]
##'
##' \dontrun{
##' ## a more realistic one:
##' 
##' postsample2 <- posteriorMCMC(Nsim=50e+3,Nbin=15e+3,
##'          dat= Leeds,
##'          prior = prior.pb,
##'          proposal = proposal.pb,
##'          likelihood = dpairbeta,
##'          Hpar=pb.Hpar,
##'          MCpar=pb.MCpar)
##' dim(postsample2[[1]])
##' postsample2[-1]
##' }
##' 
##' @keywords htest multivariate
posteriorMCMC <-
    function(prior = function(type=c("r","d"), n , 
                 par, Hpar,
                 log, dimData){
        NULL},
             
             proposal = function(type = c("r","d"),
                 cur.par,
                 prop.par,
                 MCpar, log){
                 NULL} ,
         
             likelihood = function(x, par,
                 log, vectorial){
                 NULL},
             Nsim,
             dat,
             Hpar,
             MCpar,
             Nbin=0,
             par.start=NULL,
             show.progress = floor(seq(1, Nsim,
                 length.out = 20 ) ),
             seed=NULL, kind ="Mersenne-Twister",
             save=FALSE, class=NULL,
             name.save=NULL,
             save.directory = "~", name.dat = "",
             name.model= ""
             )
  {
    ## keep track of arguments
    argnames <-ls()
    arglist <- list()
    for(i in 1:length(argnames))
      {
        arglist[[i]] <- get(argnames[i])
      }
    names(arglist) <- argnames

    
####### initialize ####
    if(!is.null(seed))
      set.seed(seed,kind=kind)
    
    start.time <- proc.time()
    p <- ncol(dat)
    if(is.null(par.start))
      {
        repeat{
          par.start <- prior (type = "r", n=1, Hpar = Hpar, dimData=p)
          condit <-
            likelihood(x=dat, par=par.start,
                       log=TRUE, vectorial=FALSE) +
                         prior(type = "d", par = par.start,
                               Hpar = Hpar, log = TRUE, dimData=p)
          if(is.finite(condit)) break
        }
      }
    cur.par <- par.start

    llh.cur <- likelihood(x=dat, par=cur.par,
      log=TRUE, vectorial=FALSE)
    
    lprior.cur <- prior(type = "d", par = cur.par,
       Hpar = Hpar, log = TRUE, dimData=p)

    nsim=1
    n.accept=0
    n.accept.kept = 0
   
    leng=length(cur.par)
    mean.res=rep(0,leng)
    emp.variance=rep(0,leng)
    emp.variance.unNorm=rep(0,leng)

    stored.vals <- matrix(0,nrow=Nsim-Nbin,ncol=leng)
    if(!is.null(names(par.start)) & length(names(par.start)) == leng){
     colnames(stored.vals) <- names(par.start)
    }
    llh <- double(Nsim-Nbin)
    lprior <- double(Nsim-Nbin)
    stored.vals[1,]=cur.par
    lprior[1] <- lprior.cur
    llh[1] <- llh.cur
############## start MCMC ############# 
    
    while(nsim<=Nsim)
      {
## show progression 
        if(any(nsim==show.progress) )
 
          {
            cat(paste("iter",nsim,
                      ": n.accepted=", n.accept, "\n",
                      sep = " " ))
          }
######### propose move         
        
        prop.par <- proposal(type="r", cur.par = cur.par, prop.par=NULL,
          MCpar=MCpar)
          ###### acceptance probability
        ratio.list <-
          lAccept.ratio(
            cur.par = cur.par,
            prop.par=prop.par,
            llh.cur=llh.cur,
            lprior.cur=lprior.cur,
            dat=dat,
            likelihood=likelihood, proposal=proposal, prior=prior,
            Hpar=Hpar, MCpar=MCpar
            )
      
        if( (is.finite(ratio.list$lrho) )&&( (ratio.list$lrho>0) || (log(runif(1))<=ratio.list$lrho) ) )
          ##then: proposal accepted: move, ie. update 
          {
            n.accept <- n.accept+1

            if(nsim > Nbin) n.accept.kept = n.accept.kept +1

            cur.par <- prop.par
            llh.cur <- ratio.list$llh
            lprior.cur <- ratio.list$lprior
          }
         
       
        
###### store results after burn-in
        
        if(nsim>Nbin) 
          {
            n <- nsim-Nbin
            new.mean.res=mean.res+1/n*(cur.par-mean.res)

            emp.variance.unNorm <- emp.variance.unNorm +
              (cur.par-new.mean.res)* (cur.par- mean.res)
            
            mean.res <- new.mean.res
            if(nsim == Nsim)
              {                  
                emp.variance <- emp.variance.unNorm / (n-1)
              }
        
            stored.vals[n,] <- cur.par
            llh[n] <- llh.cur
            lprior[n] <- lprior.cur
                       
          }                        

        nsim <- nsim+1
      }
########### end MCMC ###############
    
    end.time <- proc.time()
    print(end.time-start.time)

    res <- list(stored.vals=stored.vals,
                llh=llh,
                lprior=lprior,
                 arguments=arglist,
                elapsed =end.time-start.time, 
                Nsim=Nsim,
                Nbin=Nbin,
                n.accept=n.accept,
                n.accept.kept = n.accept.kept,
                emp.mean=mean.res,
                emp.sd=sqrt(emp.variance)
                ##est.error=sqrt(emp.variance/n)
                )
    class(res) <- c(class,"postsample")
    
    if( save && !is.null(name.save))
      {
        loglist <- list(model=name.model,
          Nsim=Nsim,Nbin=Nbin,dat=name.dat,
          Hpar=Hpar,MCpar=MCpar, seed=seed, kind=kind)
        
        name.log <- paste(name.save, ".log", sep="")
    
        assign(name.save,res)
        assign(name.log,loglist)
        save(list=name.save,
             file=paste(save.directory,"/",name.save,".rda",sep=""))
        save(list=name.log,
             file=paste(save.directory,"/", name.log,".rda",sep="")) 
        
      }
    
    return(res)
  }             

##' @export
print.postsample <- function(x, ...)
  {
    cat("\n A posterior sample of class (S3) \"postsample\" \n ")
    cat("Acceptance ratio after burn-in:",
        x$n.accept.kept/(x$Nsim- x$Nbin-1), "\n",  sep=" ")
    cat("time elapsed:\n")
    print(x$elapsed)
    cat("List names: ")
    cat(names(x), sep=" , ")
    cat("\n")
   }
