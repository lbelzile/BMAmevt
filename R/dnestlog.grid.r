##' @inheritParams dnestlog
##' @rdname dpairbeta.grid
##' @export
dnestlog.grid <- function(par,
                             npoints=50,eps=1e-3, equi = TRUE,
                             displ=TRUE, invisible=TRUE,
                             ... )
  {
  discr <- discretize(npoints=npoints,eps=eps,equi=equi)

  X.grid <- discr $X
  Y.grid <- discr $ Y
  if((length(par) !=4) || any(par>1) ||  any(par<0) )
    stop("misspecified parameters")
       
  C.out <- .C("d_trinestlog_grid", as.double(X.grid),
    as.double(Y.grid),
    as.integer(length(X.grid)),
    as.double(par[1]),
    as.double(par[2:4]), as.integer(equi), 
    result =  as.double(rep(0,npoints*npoints)) )

  density <- matrix(C.out $ result, ncol = npoints, byrow=F)

  
  mult.dens <- 1/sqrt(3)*equi + 1*(!equi)
  density <- mult.dens*density
  
  if(displ)
    {
      dgridplot(density=density,
##                npoints=npoints,
                eps=eps,
                equi=equi,
                ...)

    }

  if(invisible)
    {invisible(density)}
  else
    return(density)

  }

## dnestlog.grid <-
##   function(par=c(0.5,0.5),
##            npoints=30,eps=10^(-3), equi = FALSE,
##            displ=FALSE, invisible=FALSE,
##            ... ) 
## {
##   discr <- discretize(npoints=npoints,eps=eps,equi=equi)

##   X.grid <- discr $X
##   Y.grid <- discr $ Y

       
##   C.out <- .C("d_nestlog_grid", as.double(X.grid),
##     as.double(Y.grid),
##     as.integer(length(X.grid)),
##     as.double(par[1]),
##     as.double(par[2]), as.integer(equi), 
##     result =  as.double(rep(0,npoints*npoints)) )

##   density <- matrix(C.out $ result, ncol = npoints, byrow=F)

  
##   mult.dens <- 1/sqrt(3)*equi + 1*(!equi)
##   density <- mult.dens*density
  
##   if(displ)
##     {
##       dgridplot(density=density,
##                 npoints=npoints,
##                 eps=eps,
##                 equi=equi,
##                 ...)

##     }

##   if(invisible)
##     {invisible(density)}
##   else
##     return(density)

## }
