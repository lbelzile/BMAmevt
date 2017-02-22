
##'@inheritParams posterior.predictive3D
##' @export
##' @rdname posterior.predictive.pb
posterior.predictive.nl <-
  function(
    post.sample , ## from pb.post.sample  
    from = post.sample$Nbin+1, to = post.sample$Nsim,
    thin = 50,  
    npoints=40,eps=1e-3, equi = T, displ=T,
    ...
    )

  {

    predictive <-
      posterior.predictive3D(post.sample = post.sample , 
                             densityGrid =dnestlog.grid,
                             from = from, to = to, thin = thin,  
                             npoints=npoints, eps=eps, equi = equi,
                             displ=displ,
                             ...
                             )

    return(predictive)


  }
##   discr=discretize(npoints=npoints,eps=eps,equi=equi) 
##   X.grid = discr $X
##   Y.grid = discr $ Y

##   Nsim = post.sample$Nsim
##   Nbin = post.sample $ Nbin
  
##   time.idx = (from-Nbin):(to-Nbin)
##   kept.idx = time.idx[ (time.idx %% thin == 0) ]
##   mat.postSample = post.sample $ stored.vals[kept.idx,] 

  
##   density.grid.fun=function(v)
##     {
##       C.out = .C("d_nestlog_grid", as.double(X.grid),
##         as.double(Y.grid),
##         as.integer(npoints),
##         as.double(v[1]),
##         as.double(v[-1]),
##         as.integer(equi), 
##         result =  as.double(rep(0,npoints*npoints)) )

##       return(C.out $ result)
  
##     }

##   ## bindvect.res = apply(mat.postSample,1,density.grid.fun)

##   ## mean.vect.res = apply(bindvect.res, 1,mean)
##   mean.vect.res=rep(0, npoints*npoints)
  
##   for(i in 1:length(kept.idx))
##     {
##       mean.vect.res=mean.vect.res+
##         density.grid.fun(mat.postSample[i,])
##     }
  
##   mean.vect.res=mean.vect.res/length(kept.idx)


##   pred.density = matrix(mean.vect.res, ncol =npoints, byrow=T)
##   if(displ)
##     {plot.dens(Density = pred.density,
##             npoints = npoints, 
##             eps = eps,
##             equi = equi,
##              ... )
##    }
##   return(pred.density)
## }

