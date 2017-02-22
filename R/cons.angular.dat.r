##' Builds an angular data set, retaining the points with largest radial component.
##' 
##' The data set \code{frechetDat} is assumed to be marginally unit  Frechet distributed.
##' @title Angular data set generation from unit Frechet data.
##' @param coordinates Index vector of the columns in \code{frechetDat} to be retained to construct the angular data set.
##' @param frechetDat The data set. A matrix: each row is a multivariate record. May contain \var{NA}'s. 
##' @param n The number of desired observations in the final angular data set. Should be less than \code{nrow(frechetDat)}
##' @param displ logical. Should the angular data set be plotted ? 
##' @param invisible logical. Should the result be returned as invisible ?
##' @param add logical. Only used when \code{displ==TRUE}. Should the points be added to the current plot ?
##' @param ... Additional graphical parameters and  arguments to be passed to function \code{\link[graphics]{plot.window}} and \code{\link[graphics]{points}}.
##' @return The angular data set: A \code{n*length(coordinates)} matrix, containing values between zero and one, which rows sum to one: Each row is thus a point on the unit simplex of dimension \code{length(coordinates)-1}. Returned as  invisible if \code{invisible==TRUE}.
##' @import stats
##' @export
##' @include add.frame.r
##' @include transf.to.equi.r
##' @inheritParams add.frame
##' @examples  \dontrun{cons.angular.dat()}
##' @keywords datagen manip multivariate
cons.angular.dat <-
  function(coordinates=c(1,2,3), 
           frechetDat=get("frechetdat"), 
           n=100,
           displ=TRUE,
           invisible=TRUE,
           add=FALSE,
           lab1="w1", lab2="w2",lab3="w3",npoints=60,
           col.polygon="white",
           ...)
  {
    
    ff <- na.omit(frechetDat[, coordinates])

    rr <- apply(ff, 1, sum)
  
    ww <- ff/rr
    
    ssl <- sort.list(rr, decreasing = T)
    
    sortW <- ww[ssl,]

    dat <- sortW[1:n,]
  
    if(displ)      ###&(length(coordinates)==3))
      {

        Points=(apply((dat[,1:2]),1,transf.to.equi))
        mult.x=sqrt(2)
        mult.y=sqrt(3/2)

        if(!add)
          {
            plot.new()
            plot.window( xlim=c(0,mult.x),ylim=c(0,mult.y), asp=1,
                        bty ="n",adj=1, ...)
          }

        points(Points[1,],Points[2,],... )
           
        if(!add)
          {
            ## add.frame(equi=TRUE,
            ##       lab1=lab1, lab2=lab2,lab3=lab3,npoints=npoints,
            ##       col.polygon=col.polygon)
            ## if(!equi)
            ##   {
            ##     segments(0,0,0,1)
            ##     segments(0,0,1,0)
            ##     segments(0,1,1,0)
            ##     axis(1, at = c(0, 0.25,0.5, 0.75, 1),tcl=NA,cex.lab=0.5,
            ##          line=0.1,mgp=c(2,0.5,0),cex.axis=0.8)
            ##     axis(2, at = c(0,0.25,0.5, 0.75,  1),tcl=NA,cex.lab=0.5,
            ##          line=0.1,mgp=c(2,0.5,0), cex.axis=0.8)
            ##     mtext(lab1, side = 1, line = 2, outer = FALSE,adj=0.8)
            ##     mtext(lab2, side = 2, line = 2, outer = FALSE,adj=0.8)

            ##   }
            ## else
##              {
                segments(0,0,sqrt(2),0)
                segments(0,0,sqrt(2)/2,sqrt(3/2))
                segments(sqrt(2)/2,sqrt(3/2),sqrt(2),0)
                axis(1,at=round(mult.x*c(0, 0.25, 0.5, 0.75, 1),2),tcl=NA,
                     padj=1,cex.axis=0.8,mgp=c(1.5, 0.2, 0),
                     line=0
                     )
          mtext(lab3,cex=1.2, side = 1, line = 2.5, outer = FALSE,adj=0)
          mtext(lab1,cex=1.2, side = 1, line = 2.5, outer = FALSE,adj=1)
          mtext(lab2,cex=1.2, side = 3, line = 0, outer = FALSE,adj=0.5)


  ##            }
          }
         
      }
    
    if(invisible)
      {
       return( invisible(dat))
      }
    
    return(dat)
    
  }
