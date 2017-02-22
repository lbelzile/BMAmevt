##' Adds graphical elements to the current plot (on the two-dimensional simplex).
##'
##' Generic graphical tool for obtaining nice plots of the two-dimensional simplex
##' @title Adds graphical elements to a plot of the two dimensional simplex.
##' @inheritParams discretize
##' @param lab1 Character string: label for first component.
##' @param lab2 Character string: label for second component.
##' @param lab3 Character string: label for third component.
##' @param col.polygon The background color outside the simplex.
##' @param axes logical. Should axes be added ? 
##' @examples
##' plot.new()
##' add.frame()
##' plot.new()
##' mult.x=sqrt(2); mult.y=sqrt(3/2)
##' plot.window( xlim=c(0,mult.x),ylim=c(0,mult.y), asp=1,bty ="n")
##' add.frame(equi=TRUE)
##' @keywords aplot
##' @export add.frame
add.frame <-
  function(equi=FALSE,
           lab1="w1",lab2="w2",lab3="w3",npoints=60,
           col.polygon="black", axes=TRUE)
    {
      mult.x <- sqrt(2)*equi+1*(!equi)
      mult.y <- sqrt(3/2)*equi+1*(!equi)
      mult.dens <- 1/sqrt(3)*equi+1*(!equi) 
      if(!equi)
        {
         
          ## polygon(x=c(0,1,1,0),y=c(1,0,1,1),
          ##         col=col.polygon, border=gray(0.5)) 
          polygon(x=c(-1/npoints, 1+1/npoints,1+1/npoints, -1/npoints),
                  y=c(1+1/npoints,-1/npoints,1+1/npoints,1+1/npoints),
                  col=col.polygon, border="black") 
          polygon(x=c(-1/npoints, -1/npoints,0, 0, -1/npoints),
                  y=c(1+1/npoints,-1/npoints,-1/npoints,1+1/npoints,1+1/npoints),
                  col=col.polygon, border="black") 

          polygon(x=c(-1/npoints, 1+1/npoints,1+1/npoints, -1/npoints, -1/npoints),
                  y=c(-1/npoints,-1/npoints,0,0,-1/npoints),
                  col=col.polygon, border="black") 
         
          if(axes)
            {
              axis(1, at = c(0, 0.25,0.5, 0.75, 1),tcl=NA,cex.lab=0.5,
                   line=0.1,mgp=c(2,0.5,0),cex.axis=0.8)
              axis(2, at = c(0,0.25,0.5, 0.75,  1),tcl=NA,cex.lab=0.5,
                   line=0.1,mgp=c(2,0.5,0), cex.axis=0.8)
            }
          mtext(lab1, side = 1, line = 2, outer = FALSE,adj=0.8)
          mtext(lab2, side = 2, line = 2, outer = FALSE,adj=0.8)
        }   
        
      else
        {

          polygon(x=c(-1/npoints,0,sqrt(2)/2,sqrt(2),sqrt(2)+1/npoints,
                    sqrt(2)+1/npoints,sqrt(2)/2,-1/npoints,-1/npoints),
                  y=c(0,0,sqrt(3/2),0 ,0,
                    sqrt(3/2)+1/npoints,sqrt(3/2)+1/npoints,
                    sqrt(3/2)+1/npoints,0),
                  col=col.polygon, border="black", lwd=0.3)

          polygon(x=c(-1/npoints,sqrt(2)+1/npoints,sqrt(2)+1/npoints ,
                    -1/npoints,-1/npoints),
                  y=c(-1/npoints,-1/npoints,0 ,
                    0,-1/npoints),
                  col=col.polygon, border="black", lwd=0.3)
          if(axes)
            {
              axis(1,at=round(mult.x*c(0, 0.25, 0.5, 0.75, 1),2),tcl=NA,
                   padj=1,cex.axis=0.8,mgp=c(1.5, 0.2, 0),
                   line=0
                   )
            }
            

          mtext(lab3,cex=1.2, side = 1, line = 2.5, outer = FALSE,adj=0)
          mtext(lab1,cex=1.2, side = 1, line = 2.5, outer = FALSE,adj=1)
          mtext(lab2,cex=1.2, side = 3, line = 0, outer = FALSE,adj=0.5)

        }
    }  
