##' @rdname MCpriorIntFun.pb
##' @inheritParams MCpriorIntFun
##' @export
MCpriorIntFun.nl <-
function(Nsim=200,
         FUN=function(par,...){par},
         store=TRUE,
         Hpar = get("nl.Hpar"),
         show.progress = floor(seq(1, Nsim, length.out = 20 ) ) ,
         Nsim.min=Nsim,
         precision=0,
         ...)
  {
    mcres <- MCpriorIntFun(Nsim=Nsim,
                           prior=prior.nl,
                           Hpar=Hpar,
                           dimData=3,
                           FUN=FUN,
                           show.progress=show.progress,
                           Nsim.min=Nsim.min,
                           precision=precision,
                           ...)
    return(mcres)

  }
## ############ intialize  ############

##     start.time=proc.time()
##     not.finite=0
##     param = rparameter.prior.nl(n=1, Hpar=Hpar)
##     temp.res=FUN(param,...)
##     dim.res=dim(temp.res)
    
##     if(is.null(dim.res) || (sum(dim.res!=1) ==1)  )
##       {
##         emp.mean=rep(0,length(temp.res))
##       }
##     else
##       {
##         store=FALSE
##         emp.mean=array(0,dim=dim.res)
##       }
##     emp.variance= emp.mean
##     emp.variance.unNorm=emp.variance

##     if(store)
##       {
##         stored.vals=matrix(0,nrow=Nsim,ncol=length(emp.mean))
##       }
      
## ################ start MC ########
##     n=1
##     while((n<=Nsim) &&
##           ( (n<=Nsim.min) ||
##            (max( sqrt(emp.variance/(n-1))/emp.mean) > precision) )
##           )
##       {
##         ## show progression      
##         if(any(n==be.patient) || n==50 || n==100 || n==200)
##           {
##             cat(paste((n-1), "iterations done", "\n", sep = " " ))
##           }


##         flag=TRUE
##         count = 0
##         while(flag & (count<=10))
##           {
##             param = rparameter.prior.nl(n=1, Hpar=Hpar)
##             temp.res=FUN(param,...)
##             flag = (any(sapply(as.vector(temp.res),
##                                function(x){ ! is.finite(x) } ) ) )
##             if(flag)
##               {
##                 not.finite = not.finite+1
##               }
##             count = count+1
##           }
##         if(flag)
##           stop("more than 10 non finite values produced")
        
##         cur.res=temp.res ##FUN(alpha=param[1],beta=param[2:(ll+1)],...)
        
##         new.emp.mean=emp.mean+1/n*(cur.res-emp.mean)
        
##         emp.variance.unNorm=emp.variance.unNorm +
##           (cur.res-new.emp.mean)* (cur.res- emp.mean)
##         emp.variance = emp.variance.unNorm/(n-1)
##         emp.mean = new.emp.mean
        
##         if(store)
##           {
##             stored.vals[n,]= as.vector(cur.res)
##           }
##         n=n+1
##       }
## ######### end MC #########
##     end.time = proc.time()
##     elapsed=end.time-start.time
##     print(elapsed)
##     if(store)
##       {
##         return(list( stored.vals=stored.vals[1:(n-1), ],
##                     elapsed=elapsed,
##                     nsim = n-1,
##                     emp.mean=emp.mean,
##                     emp.stdev=sqrt(emp.variance),
##                     est.error=sqrt(emp.variance/(n-1)),
##                     not.finite = not.finite))
##       }
##     else
##       {
##         return(list(  nsim=n-1,
##                     elapsed=elapsed,
##                     emp.mean=emp.mean,
##                     emp.stdev=sqrt(emp.variance),
##                     est.error=sqrt(emp.variance/(n-1)),
##                     not.finite = not.finite ))
                  
##       }
##   }
