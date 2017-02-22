#include <assert.h>

#include "header.h"


double baseV_trinestlog(double *x, double alpha, double *beta)
{
  double U1 = exp( beta[0]*log(pow(x[0], -1./(alpha*beta[0]) ) + 
				pow(x[1], -1./(alpha*beta[0]) )) )  ;

  double U2 = exp( beta[1]*log(pow(x[0], -1./(alpha*beta[1]) ) + 
				pow(x[2], -1./(alpha*beta[1]) )) ) ;

  double U3 = exp( beta[2]*log(pow(x[1], -1./(alpha*beta[2]) ) + 
				pow(x[2], -1./(alpha*beta[2]) )) ) ;
 
  return U1+U2+U3 ;
}

void expfunction_trinestlog(double *x, double *par, double *result)
{
  double unpowered = baseV_trinestlog(x, par[0], par+1);
  *result = pow(unpowered, *par)/ pow(2,*par);
}



double dx_baseV_trinestlog(double *x, double alpha, double *beta, int coord)
{
  int i = coord-1 ; 
  int j,k, ij, ik ;
  double T1, T2,total ;
  if(i==0) /* d/dw1 ; j=w2 k = w3 */
    {
      j=1;  k=2; 
      ij = 0 ;  ik = 1 ;
    }
  if(i==1) /* d/dw2 ; j=w3 k = w1 */
    {
      j=2; k=0;
      ij = 2;  ik = 0 ;
    }
  if(i==2) /* d/dw3 ; j=w1 k = w2 */
    {
      j=0; k=1;
      ij = 1 ; ik = 2 ;
    }
  T1 = pow(x[i], -1./(alpha*beta[ij])-1.) * 
    exp( (beta[ij]-1)*log (pow(x[i], -1./(alpha*beta[ij]))  + 
			  pow(x[j], -1./(alpha*beta[ij])) ) ) ;

  T2 = pow(x[i], -1./(alpha*beta[ik])-1.) * 
    exp( (beta[ik]-1.)*log (pow(x[i], -1./(alpha*beta[ik]))  + 
			  pow(x[k], -1./(alpha*beta[ik])) ) ) ;
    
  total = -1./alpha * (T1 +  T2 ) ;
  return total ;
}

double dxy_baseV_trinestlog(double *x, double alpha, double *beta, int *coord)
{
  int i = coord[0]-1 ; 
  int j = coord[1] - 1 ;
  assert(i != j);

  int  ij=0 ;
  double T1,F1,total ;
  if(i==0) /* d/dw1 */
    {
      if(j == 1) /*; d/dw2  ij = beta_1 */
	{
	  ij = 0 ; 
	}
      if(j==2) /*; d/dw3  ij = beta_2 */
	{
	  ij = 1;
	}
    }
  if(i==1) /* d/dw2 */
    {
      if(j == 0) /*; d/dw1  ij = beta_1 */
	{
	  ij = 0 ; 
	}
      if(j==2) /*; d/dw3  ij = beta_3 */
	{
	  ij = 2;
	}

    }
  if(i==2) /* d/dw3 */
    {
      if(j == 0) /*; d/dw1  ij = beta_2 */
	{
	  ij = 1 ; 
	}
      if(j==1) /*; d/dw2  ij = beta_3 */
	{
	  ij = 2;
	}
    }
  T1 =  /*pow(x[i], -1./(alpha*beta[ij]-1)) * */
    exp( (beta[ij]-2.)*log (pow(x[i], -1./(alpha*beta[ij]))  + 
			  pow(x[j], -1./(alpha*beta[ij])) ) ) ;

  F1 = pow(x[i]*x[j], -1./(alpha*beta[ij]) - 1.) ;
  total = ( (beta[ij] - 1. ) / (pow(alpha,2.)*beta[ij]) )* F1 * T1 ;
  return total ;
}




double d_trinestlog_point(double *x, double alpha, double *beta, 
			  int take_logs)
{
  double constant = pow(2., - alpha)/3. * alpha * (alpha-1.);
  double G =  baseV_trinestlog(x,  alpha, beta);
  double dxG= dx_baseV_trinestlog(x, alpha, beta, 1);
  double dyG= dx_baseV_trinestlog(x, alpha, beta, 2);
  double dzG= dx_baseV_trinestlog(x, alpha, beta, 3);
  int xy[2]; 
  xy[0]=1;
  xy[1]=2;
  int xz[2];
  xz[0]=1;
  xz[1]=3;
  int yz[2];
  yz[0]=2;
  yz[1]=3;
  double dxdyG = dxy_baseV_trinestlog(x, alpha, beta, xy);
  double dxdzG = dxy_baseV_trinestlog(x, alpha, beta, xz);
  double dydzG = dxy_baseV_trinestlog(x, alpha, beta, yz);
    
  double F1 = pow(G, alpha-2.);
  double T1 = dxG*dydzG + dyG * dxdzG + dzG * dxdyG;
  double F2 = (alpha-2.)*pow(G, alpha-3.);
  double T2 = dxG * dyG * dzG ;

  double result = - constant * (F1*T1 + F2*T2 );
  if(take_logs)
    return log(result);
  else
    return result;
}

void test1(double *x, double *alpha, double *beta, 
int *take_logs, double *result)
{
  /*  double x[3];
  x[0]=0.1;
  x[1]= 0.3;
  x[2]=0.6;
  double alpha=0.1;
  double beta[3]; 
  beta[0]=0.1;
  beta[1]= 0.3;
  beta[2]=0.6;
  int coord[2];
  coord[0]=1;
  coord[1]=2;*/
  *result  = d_trinestlog_point(x, *alpha, beta, *take_logs);
}


void d_trinestlog(double *x, int*pnx, 
		  double *alpha, double *beta, int *take_logs, 
		  int *return_vector, double*result)
{
  int n ;
  double *index_x = x ;
  double *index_res = result ;
  int isOnSimplex ;
  double point_lDensity ;
  int onePointOut = 0 ;
  *index_res = 0 ;
  for( n=0 ; n < *pnx ; n++ )
    {
      isOnSimplex = is_on_simplex(index_x, 3) ;
      if( ! isOnSimplex)
	{
	  (*index_res ) = (*take_logs) ? log(ZERO_BMAMEVT) : 0 ;
	  onePointOut = 1 ;
	}

      else
	{
	  point_lDensity =
	    d_trinestlog_point(index_x, *alpha, beta, 1);
	  
	  if(*return_vector)
	    {
	      (*index_res ) = (*take_logs) ? 
		point_lDensity : exp(point_lDensity)  ;
	      /* index_res ++ ;*/
	    }
	  else
	    (*index_res) += point_lDensity ;
	    
	}

      index_x +=  3 ;
      if(*return_vector) index_res++ ;
    }
  
  if( ( !(*take_logs) ) && (! (*return_vector)) )
    *index_res = exp( *index_res) ;

  if (onePointOut &&  (! (*return_vector) ))
    (*index_res ) = (*take_logs) ?  log(ZERO_BMAMEVT) : 0  ;
}
