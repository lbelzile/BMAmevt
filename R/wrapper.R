
wrapper <-
function(x, y, FUN,...)    #internal
      {
       sapply(seq_along(x), FUN = function(i) FUN(x[i], y[i],...))
      }

