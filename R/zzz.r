
##' @import utils
##' @import coda
##' @useDynLib BMAmevt, .registration=TRUE
.onAttach=function(libname,pkgname = "BMAmevt")
{
## data(frechetdat)
## data(Leeds)
## data(nl.Hpar)
## data(nl.MCpar)
## data(pb.Hpar)
## data(pb.MCpar)
## data(dm.expar.D3k3)
## data(winterdat)
set.seed(1)
}




.onUnload = function(libpath)
  {
    library.dynam.unload(chname="BMAmevt",libpath)
  }
