## @param n The number of points on the simplex to be generated.

##' @rdname dpairbeta
##' @param threshold The radial threshold \eqn{r} above which the simulated points should be kept to build the angular dataset. Should be set to a high value,  for the asymptotic approximation  
##' \deqn{P(W \in B |\; ||X|| >r)\simeq H(B)}{P(W \in B |  ||X|| >r) ~ H(B)} to hold.
##' @param return.points logical: should the censored vectorial dataset corresponding to the angular one be returned ?
##' @export 
rnestlog <- function(n=5,par=c(0.2,0.3,0.4,0.5),
                        threshold=1000,
                        return.points=FALSE)
  {
    Angles <- matrix(0,ncol=3,nrow=n)
    Points <- Angles
    alpha <- par[1]
    beta <- par[2:4]
        for(i in 1:n)
      {
        repeat
          {
            s <- rstable.posit(alpha=alpha)
            s12 <- rstable.posit(alpha=beta[1])
            s13 <- rstable.posit(alpha=beta[2])
            s23 <- rstable.posit(alpha=beta[3])

            X1.12= (s12*(s/2)^(1/beta[1])/rexp(1) )^(alpha*beta[1]) 
            X1.13= (s13*(s/2)^(1/beta[2])/rexp(1) )^(alpha*beta[2])

            X1 = max(X1.12,X1.13)

            X2.12= (s12*(s/2)^(1/beta[1])/rexp(1) )^(alpha*beta[1]) 
            X2.23= (s23*(s/2)^(1/beta[3])/rexp(1) )^(alpha*beta[3])

            X2 = max(X2.12,X2.23)

            X3.13= (s13*(s/2)^(1/beta[2])/rexp(1) )^(alpha*beta[2]) 
            X3.23= (s23*(s/2)^(1/beta[3])/rexp(1) )^(alpha*beta[3])

            X3 = max(X3.13,X3.23)
            
            X=c(X1,X2,X3)

            if(sum(X)>threshold)
              {
                break
              }
          }
     
        Angle=X/sum(X)
        Angles[i,]=Angle
        Points[i,]=X
      }
    
      if(return.points)
        {
          return(list(Angles=Angles,Points=Points))
        }  
      else
        {
          return(Angles)
        }



  }






## rnestlog <-
##   function(n=3, par=c(0.5,0.5),
##            threshold=1000, return.points=FALSE)
##   {
##     Angles <- matrix(0,ncol=3,nrow=n)
##     Points <- Angles
##     alpha0 <- par[1]
##     alpha12 <- par[2]
##     for(i in 1:n)
##       {
##         repeat
##           {
##             s0=rstable.posit(alpha=alpha0)
##             s12=rstable.posit(alpha=alpha12)
          
##             E1=rexp(1)
##             E2=rexp(1)
##             E3=rexp(1)
##             X1=(s0^(1/alpha12 )*s12 / E1)^(alpha0*alpha12) #^(1/alpha12 )
##             X2=(s0^(1/alpha12)*s12 / E2)^(alpha0*alpha12) #s0^(1/alpha12)
##             X3= (s0/E3)^(alpha0)
##             X=c(X1,X2,X3)
##             if(sum(X)>threshold)
##               {
##                 break
##               }
##           }
     
##         Angle=X/sum(X)
##         Angles[i,]=Angle
##         Points[i,]=X
##       }
    
##       if(return.points)
##         {
##           return(list(Angles=Angles,Points=Points))
##         }  
##       else
##         {
##           return(Angles)
##         }

##   }
                     
