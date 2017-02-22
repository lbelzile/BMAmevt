
##' @export
##' @rdname dpairbeta
dnestlog <- function(x=rbind(c(0.1,0.3,0.6),c(0.3,0.3,0.4)) ,
        par=c(0.5,0.5,0.2,0.3),
        log=FALSE, vectorial = T)
  {
  xvect = as.double(as.vector(t(x) ))
  
  if(is.vector(x))
    { dim = as.integer( length(x) ) 
      n = as.integer(1)
    }
  
  else
    {
      dim = as.integer(ncol(x) )
      n = as.integer(nrow(x) )
    }

    if((length(par) !=4) || any(par>1) ||  any(par<0) )
    stop("misspecified parameters")

  alpha <- as.double(par[1])
  beta <- as.double(par[2:4])

  if(vectorial)
    {
      result=double(n)
    }
  else
    {
      result=double(1)
    }

  C.out = .C("d_trinestlog", x=xvect, pnx = n, 
    alpha = alpha, beta = beta,
    take_logs = as.integer(log),
    return_vector = as.integer(vectorial),
    result = result)
##browser()
  return(C.out$result)
  }

## dtrinestlog(x=rbind(c(0.1,0.3,0.6),c(0.3,0.4,0.3)) ,
##         par=c(0.5,0.1,0.4,0.8),
##         log=FALSE, vectorial = T)
