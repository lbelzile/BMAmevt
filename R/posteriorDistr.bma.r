##' Builds an empirical distribution defined as a sum of weighted Dirac masses from posterior samples in individual models.
##'
##' @title Posterior distribution in the average model
##' @param postweights a vector of positive real numbers, summing to one: the posterior weights (in the same order as the elements of \code{post.distrs}) of the  individual models.
##' @param post.distrs A list of same length as \code{postweights}.
##' Each element must be a vector which will be used as a posterior sample.
##' @return A matrix with two rows and as many columns as the sum of the lengths of the elements of \code{post.distrs}. The second line contains the weighted  posterior sample in the BMA;  the first line contains the weights to be assigned to each corresponding value on the second one.
##' @export
posteriorDistr.bma <- function(postweights=c(0.5,0.5),
                              post.distrs=list()
                              )
{
  M <- length(post.distrs)
  sizes <- double(M)
  sortedlist <- list()
  for (i in 1:M)
    {
      sizes[i] <- length(post.distrs[[i]])
      sortedlist[[i]] <- sort(post.distrs[[i]], decreasing=FALSE)
    }
  Nunion <- sum(sizes)
  sortedsample <- double(Nunion)
  sortedweights <- double(Nunion)
  for(i in 1:Nunion)
    {
      candidates <- sapply(1:M, function(m){sortedlist[[m]][1]}
                           )
      minind <- which.min(candidates)
      sortedsample[i] <- candidates[minind]
      sortedweights[i] <- postweights[minind]/sizes[minind]
      sortedlist[[minind]] <- sortedlist[[minind]][-1]
    }
  res <- rbind(sortedsample,sortedweights)
  rownames(res) <- c("values","weights")
  return(res)
}




## posteriorMean.bma <- function(postweights=c(0.5,0.5),
##                               post.samples=list(),
##                               FUN=function(par,model,...){NULL},
##                               models=c("pairbeta", "trinestlog"),
##                               N=1e+3,
##                               froms=c(NULL,NULL),
##                               tos=c(NULL,NULL),
##                               ...)
##   {
##     if( (length(postweights) != length(models) )| (length(models) != length(froms) ) | (length(froms) != length(tos) ) )
##        {stop("arguments 'postweights', 'models', 'froms'  and 'tos' should have  same length. ")
##       }
##     Nmodel <- length(postweights)
    
##     for( m in 1:Nmodel)
##       {
##         if(is.null(froms[m]))
##           froms[m] <- post.samples[[m]]$Nbin +1
##         if(is.null(tos[m]))
##           {
##             tos[m] <- post.samples[[m]]$Nsim
##           }
##         if((froms[m]<= post.samples[[m]]$Nbin) |
##            (tos[m]>post.samples[[m]]$Nsim) )
##           {stop(" argument \"from\" or \"to\" out of range")}
##         }

    

##     idx <- 1:N
##     cumweights <- cumsum(postweights)
##     intern.fun <- function(i)
##       {
##         u <- runif(1)
##         m <- min((1:Nmodel)[u<cumweights])
##         ind.i <- sample(froms[m]:tos[m],1)
##         par.i <- (post.samples[[m]])$stored.vals[ind.i,]
##         res.i <- as.vector(FUN(par=par.i, model=models[m], ...) )
##         return(res.i)
##       }

##     values <- sapply(idx,
##                      intern.fun,
##                      simplify=TRUE
##                      )
##     if(is.vector(values))
##       values <- matrix(values,nrow=1)
    
##     est.mean <- apply(values, 1, mean)
##     est.sd <- apply(values,1,sd)
##     if(displ)
##       {
##         for(i in 1:nrow(values))
##           {
##             dev.new()
##             ylim <- range(values[i,])
##             cummean <- cumsum(values[i,])/(1:N)
##             plot(kept.idx, cummean, type="l" )
##             esterr <- sqrt(cumsum( (values[i,]- cummean)^2 ) )/1:N
##             lines( cummean+2*esterr, col=gray(0.5) )
##             lines( cummean-2*esterr, col=gray(0.5) )
##           }

##       }

##     return(list(values=values, est.mean=est.mean, est.sd=est.sd))
##   }

