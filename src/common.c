#include "header.h"

int is_on_simplex(double *x, int dim)
{
  double sum, *xend;
  sum = 0;
  for (xend = x+dim ;  x != xend; x++)
	{
	  if (*x <= 0) return 0;
	  sum += *x;
	}
  if (fabs(sum - 1) > TOL_BMAMEVT) return 0;
  return 1;
}


void transf_to_rect(double *vectequi)
{
  double x= vectequi[0] ;
  double y= vectequi[1] ;
  *vectequi = 1./sqrt(3.) * (sqrt(3./2.) * x - sqrt(2.)/2. * y  ) ;
  *(vectequi+1) = sqrt(2./3.) * y ;
}
