##' Plots contours or gray-scale level sets  of a spectral density on the two-dimensional simplex.
##'
##' The function  interprets the \code{density} matrix as 
##' \code{\link[graphics]{contour}} does, \emph{i.e.} as a table of
##' \code{f(X[i], Y[j])} values, with column 1 at the bottom,
##' where \code{X} and \code{Y} are
##' returned by \code{\link{discretize}} and \code{f} is the
##' density function.
##' @title Image and/or Contour plots of  spectral densities in  trivariate  extreme value  models
##' @inheritParams discretize
##' @inheritParams add.frame
##' @param density A \code{npoints*npoints} matrix containing the
##' density's values  scattered on the discretization grid defined by
##' \code{npoints, equi, eps} (see \code{\link{discretize}}).
##' @param add Logical. Should the contours be added to a currently active plot ?
##' @param breaks Set of breakpoints for the gray scale colors.
##' See \code{\link[graphics]{image}}
##' @param levels Levels to which plot the contour lines. See \code{\link[graphics]{contour}}
##'  @param labcex \code{cex} for contour labeling.
##' See  \code{\link[graphics]{contour}}.
##' @param col.lines The color to be used for the contour lines.
##' @param background Logical. Should a the background be filled
##' inside the simplex \emph{via} a call to
##' \code{\link[graphics]{image}} ? 
##' @param ... Additional graphical parameters and arguments to be passed
##' to \code{\link[graphics]{contour}}  and \code{\link[graphics]{image}}.
##' @examples 
##' wrapper <- function(x, y, my.fun,...) 
##'       {
##'        sapply(seq_along(x), FUN = function(i) my.fun(x[i], y[i],...))
##'       }
##' 
##' grid <- discretize(npoints=40,eps=1e-3,equi=FALSE)
##' 
##' Density <- outer(grid$X,grid$Y,FUN=wrapper,
##'                 my.fun=function(x,y){10*((x/2)^2+y^2)*((x+y)<1)})
##' 
##' dgridplot(density= Density,npoints=40, equi=FALSE)
##' @import graphics
##' @export
dgridplot <-
  function(density= matrix(5*sin(1/73*(1:(40*40)))^2,
               ncol=40, nrow=40),
#           npoints=40,
           eps=10^(-3),
           equi=TRUE,
           add=FALSE,
           breaks=seq(-0.01,5.1,length.out=1000), 
           levels=seq(0,6, length.out=13),
           col.lines="black", labcex= 0.8,
           background=FALSE,
           col.polygon=gray(0.5),
           lab1="w1", lab2="w2", lab3="w3",
           ...)
  {
    npoints <- dim(density)[1]
    discr=discretize(npoints=npoints,eps=eps,equi=equi) 
    X_grid=discr$X
    Y_grid=discr$Y
    mult.dens=1##  1/sqrt(3)*equi+1*(!equi)
    mult.x=sqrt(2)*equi+1*(!equi)
    mult.y=sqrt(3/2)*equi+1*(!equi)
    
    if(add)
      {
        contour(x = X_grid,
                y = Y_grid,
                density *mult.dens,
                levels=levels,
                labels = NULL,col=col.lines,labcex=labcex,
                add=TRUE,
                ...)
         
        add.frame (equi=equi,
                   lab1="",lab2="",lab3="",npoints=npoints,
                   col.polygon=col.polygon)

      }
    else
      {
        plot.new()
        plot.window( xlim=c(0,mult.x),ylim=c(0,mult.y), asp=1,
                    bty ="n",adj=1, mgp=c(1.5, 0.4, 0) )

        
        if(background)
          {
            image(X_grid,Y_grid,density*mult.dens,breaks=breaks, #
                  col=gray((1:(length(breaks)-1))/(length(breaks)-1)), 
                  xlab="",
                  ylab="",
                  cex.axis=0.8,
                  xaxt ="n", yaxt = "n",
                  bty ="n",adj=1,asp=1,
                  ...)
              
            contour(x = X_grid,
                    y = Y_grid,
                    density*mult.dens ,
                    levels=levels,
                    labels = NULL,
                    labcex=labcex,
                    col=col.lines,
                    add=TRUE,...)

          }
        else
          {
            contour(x = X_grid,
                    y = Y_grid,
                    density*mult.dens ,
                    xlab="",
                    ylab="",
                    xaxt ="n", yaxt = "n",
                    bty ="n",  adj=1,asp=1,
                    levels=levels,
                    labels = NULL,
                    labcex=labcex,
                    col=col.lines,
                    add=FALSE,...)
          }
        add.frame(equi=equi,
                  lab1=lab1, lab2=lab2,lab3=lab3,npoints=npoints,
                  col.polygon=col.polygon)
           

      }
  }  
