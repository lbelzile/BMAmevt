##' Only valid in the tri-variate case
##'
##' @title Plots the Dirichlet mixture density on a discretization grid 
##' @inheritParams ddirimix
##' @inheritParams discretize
##' @param marginal logical. If \code{TRUE}, the angular density corresponds to the marginal intensity measure, over coordinates \code{coord}. Otherwise, it is only the projection of the full dimensional angular measure (hence the moments constraints is not satisfied anymore). 
##' @param coord A vector of size 3: 
##' the indices of the coordinates upon which the marginalization is to be done.
##' @param invisible Logical: should the result be returned as invisible ?
##' @param displ Logical: should a plot be issued ?
##' @param ... Additional arguments to be passed to
##' \code{\link{dgridplot}}.
##' @return The discretized density
##' @export
ddirimix.grid <-
  function(par= get("dm.expar.D3k3"),
           wei=par$wei,
           Mu= par$Mu,
           lnu=par$lnu,
           npoints=30,
           eps=10^(-3),
           equi=TRUE,
           marginal=TRUE, coord=c(1,2,3),
           invisible=TRUE,
           displ=TRUE,
           ...
           )
  {

    if(length(coord) != 3)
      {warning("coord is not  of length 3;  density projection/marginalisation will be made on the two-dimensional face corresponding to the first two indices.")
       if(length(coord)<3)
         stop("at least 3 coordinates should be specified")
       
       coord <- coord[1:3]
     }
    
    if(!marginal & (nrow(Mu)>3))
      {
        equi <- FALSE
        mu3 <- 1 - Mu[coord[1],] - Mu[coord[2],]

        Mu <- rbind(matrix(Mu[coord[1:2],], nrow=2),  matrix(mu3, nrow = 1))
      
      }
      
    else
      {
        MMu <- matrix(Mu[coord,], nrow=3)
        p <- nrow(Mu)
        multconst <- apply(MMu,2,sum)
        if(length(multconst)>1)
          Mu <- MMu %*% diag(1/multconst)
        else
          Mu <- MMu/multconst
        
        lnu <- lnu + log(multconst)
        wei <- p/3*wei*multconst
      }
      
    
    discr <- discretize(npoints=npoints,eps=eps,equi=equi)

    X.grid <- discr $X
    Y.grid <- discr $ Y
    nu <- exp(lnu)
    k <- length(lnu)
    Density <- double(npoints*npoints)

    C.out <- .C("ddirimix_grid", as.double(X.grid), as.double(Y.grid),
                as.integer(npoints), as.double(Mu),
                as.integer(k), as.double(wei), as.double(nu),
                as.integer(equi), result=Density) 

    Density <- matrix(C.out $ result, ncol = npoints, byrow=F)
    
    mult.dens <- 1/sqrt(3) * equi + 1 * (!equi)
                             
    if(displ)
      {
        dgridplot(density=mult.dens*Density,
##                  npoints=npoints,
                  eps=eps,
                  equi=equi,
                  ...)
      }

    if(!invisible)
      return(mult.dens*Density)
      
    else
      invisible(mult.dens*Density)
  }

 
      ## proj.fun=function(u,v)
      ##   {
      ##     if(!equi)
      ##       {
      ##         if(u<0 || v<0 || u+v>1)
      ##           return ( 0 )
              
      ##         else
      ##           return(ddirimix(x=c(u,v,1-(u+v)),
      ##                           wei=wei,
      ##                         Mu=Mu,
      ##                           lnu=lnu,
      ##                           log=FALSE,
      ##                           vectorial=FALSE))
      ##       }
          
      ##     else
      ##       {
      ##         W=transf.to.rect(c(u,v)) 
      ##         if(W[1]<0 || W[2]<0 || sum(W)>1)
      ##           return(0)

      ##         else 
      ##         return (ddirimix(x=c(W[1],W[2],1-sum(W)),
      ##                          wei=wei,
      ##                          Mu=Mu,
      ##                          lnu=lnu,
      ##                          log=FALSE,
      ##                          vectorial=FALSE))
      ##       }
      ##   }
        

      ## discr <- discretize(npoints=npoints,eps=eps,equi=equi) 
      ## X_grid <- discr$X
      ## Y_grid <- discr$Y
      ## wrapper <-
      ##   function(x, y, my.fun,...)    #internal
      ##     {
      ##       sapply(seq_along(x), FUN = function(i) my.fun(x[i], y[i],...))
      ##     }

      
      ## Density= outer(X_grid,Y_grid,
      ##   FUN=wrapper,
      ##   my.fun=proj.fun)

## @param project Logical: if \code{TRUE} (default) the density is the projection of a higher dimensional dirichlet mixture. Otherwise, it is the angular density on S3 associated to the  exponent measure( i.e. intensity of the limiting point process) obtained by marginalisation of the initial exponent measure. In such a case, integration is done over all directions but the one indicated by \code{coord}. 
