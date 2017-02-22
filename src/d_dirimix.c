#include "header.h"


/* static int is_on_simplex(double *x, int dim, int n) */
/* { */
  
/*   int i; */
/*   double sum, *xend; */

/*   for (i = 0; i < n; i++) */
/*     { */
/*       sum = 0; */
      
/*       for (xend = x + dim; x != xend; x++) */
/* 	{ */
/* 	  if (*x < 0) return 0; */
/* 	  sum += *x; */
/* 	} */

/*       if (fabs(sum - 1) > TOL_DIRIREV ) return 0; */
/*     } */
  
/*   return 1; */
/* } */

double ddirimix_point(double *mu, double *nu, double *x,
			  double *w, int dim, int nmu,
			  int take_logs)
 /* Returns the density itself (not the log) */
{
  double density=0;
  int m, i, j ;
  double log_const, log_dens;
  int isOnSimplex;

  isOnSimplex = is_on_simplex(x,  dim) ;
  if( ! isOnSimplex)
    {
      density = take_logs ? ZERO_BMAMEVT : 0 ;
      return density ;
    }
      
  for (m = 0; m < nmu; m++)
    {
      log_const = lgammafn(nu[m]);

      for (i = dim * m; i < dim * (m + 1); i++)
	log_const -= lgammafn(mu[i] * nu[m]);
      /*
	for (ix = 0; ix < nx; ix++)
	{*/
      log_dens = 0;
      for (j = 0; j < dim; j++)
	{
	  log_dens += 
	    (mu[j + dim * m] * nu[m] - 1) * log(x[j]);
	}
	      
      density += w[m] * exp(log_const + log_dens);
    }
      /*}*/
    
 
  return density;
}





 double *compute_ddirimix(double *mu, double *nu, double *x,
			  double *w, int dim, int nmu, int nx,
			  int take_logs)
 /* Returns the densities themselves (not the log) as a vector*/
{
  double *density = calloc(nx+1, sizeof(double));
  if (!density) return NULL;

  int m, i, j, ix;
  double log_const, log_dens;
  int isOnSimplex;
  double  oneOut = 0;
  for (m = 0; m < nmu; m++)
    {
      log_const = lgammafn(nu[m]);

      for (i = dim * m; i < dim * (m + 1); i++)
	log_const -= lgammafn(mu[i] * nu[m]);

      for (ix = 0; ix < nx; ix++)
	{
	  isOnSimplex = is_on_simplex(x + ix * dim,  dim) ;
	  if( ! isOnSimplex)
	    {
	      density[ix]= take_logs ? ZERO_BMAMEVT : 0 ;
	      oneOut = 1;
	    }
	  else
	    {
	      log_dens = 0;
	      
	      for (j = 0; j < dim; j++)
		{
		  log_dens += 
		    (mu[j + dim * m] * nu[m] - 1) * log(x[j + ix * dim]);
		}
	      
	      density[ix] += w[m] * exp(log_const + log_dens);
	    }
	}
    }
  density[nx] = oneOut;
  return density;
}


void d_dirimix(double *x, int *pn, int *pk, int *pp,
	      double *wei, double *mu, double *nu, 
	      int* take_logs, int* return_vector, int *is_error,
	      double *result)
{
  double *density;

  /* if (!is_on_simplex(x, *pp, *pn) || */
  /*       !is_on_simplex(mu, *pp, *pk) ) */
  /*   { */
  /*     *is_error = 1; */
  /*     return; */
  /*   } */
  if(!(density = compute_ddirimix(mu, nu, x, wei, *pp, 
				  *pk, *pn, *take_logs)))
    {
      *is_error=1;
      return;
    }
  int i;

  if (*return_vector)
    {
      for (i = 0; i < *pn; i++)
	  result[i] = *take_logs ? log(density[i]) : density[i];
    }

  else
    {
      *result = 0;
      for (i = 0; i < *pn; i++) 
	*result += log(density[i]);
      
      if (!(*take_logs)) *result = exp(*result);

      if(density[*pn]>0.5)
	*result = *take_logs ? log(ZERO_BMAMEVT) : 0  ;

    }

  free(density);
}
