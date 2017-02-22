
##' @rdname transf.to.equi
##' @export transf.to.rect
transf.to.rect=function(vect)
    {
      
     M=1/sqrt(3)*rbind(c(sqrt(3/2),-sqrt(2)/2),
                        c(0,sqrt(2))
                        )
      return(as.vector(
                         M%*%t(t(vect))
                       )
            )   
    
    }
