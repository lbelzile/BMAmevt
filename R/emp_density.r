emp_density<-
  function(dat, L=NULL)  ##density with respect to Lebesgue
  ##on the projected    unit  simplex
  {ndat=nrow(dat)
   if(is.null(L))
     {
       L=max(5,floor(sqrt(ndat/100)))  #200
     }
   Mdat=L*dat[,c(1,2)]
   Dens=matrix(0,nrow=L,ncol=L)
   for(k in (1:ndat))
     {
       i0=ceiling(Mdat[k,1])
       j0=ceiling(Mdat[k,2])
       Dens[i0,j0]= Dens[i0,j0]+1#*((i0+j0)<=L)
     }
   return(list(l_grid=L,
               emp_dens=Dens/(ndat)*L^2 )) 
 }
