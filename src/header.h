#ifndef BMAMEVT_HEADER_H
#define BMAMEVT_HEADER_H

#include <stdlib.h>
#include <Rmath.h>

#define TOL_BMAMEVT  1e-8
#define ZERO_BMAMEVT 1e-200


int is_on_simplex(double *x, int dim); 
void transf_to_rect(double *vectequi);

double compute_logConst_pairbeta (double alpha, int dim) ;

double unNorm_density_point_pairbeta(double alpha, double *beta, 
			    double *x, int dim );

double *compute_dpairbeta(double alpha, double *beta, double *x,
			int dim,  int nx, int take_logs);


void expfunction_nestlog(double *x, double *par, double *result); 

double logd_nestlog_point(double *x, double alpha0, double alpha12 );

void d_nestlog( double *x, int *pnx, int *pdim, 
		double *palpha0, double *palpha12,
		int *take_logs, int *return_vector,
		double *result) ;

void d_nestlog_grid(double *Xgrid, double *Ygrid, 
		      int *npoints, double *palpha0, double *palpha12,
		    int *equi , double *result);

double baseV_trinestlog(double *x, double alpha, double *beta);

void expfunction_trinestlog(double *x, double *par, double *result);

double dx_baseV_trinestlog(double *x, double alpha, double *beta, 
			   int coord);
double dxy_baseV_trinestlog(double *x, double alpha, double *beta, 
			    int *coord);

double d_trinestlog_point(double *x, double alpha, double *beta, int log);

void d_trinestlog(double *x, int*pnx, 
		  double *alpha, double *beta, int *take_logs, 
		  int *return_vector, double*result);

void d_trinestlog_grid(double *Xgrid, double *Ygrid, 
		      int *npoints, double *palpha,
		      double *pbeta, int *equi , 
		       double *result);


double ddirimix_point(double *mu, double *nu, double *x,
		      double *w, int dim, int nmu,
		      int take_logs);


double *compute_ddirimix(double *mu, double *nu, double *x,
			 double *w, int dim, int nmu, int nx,
			 int take_logs);


void d_dirimix(double *x, int *pn, int *pk, int *pp,
	       double *wei, double *mu, double *nu, 
	       int* take_logs, int* return_vector, int *is_error,
	       double *result);

void ddirimix_grid(double *Xgrid, double *Ygrid,
		   int *npoints, double *Mu, int*nmu,  double *wei,
		   double *nu, int*equi, double*result);


void ddirimix_grid1D(double *Xgrid, 
		   int *npoints, double *Mu,  int *nmu, double *wei,
		     double *nu,  double *result);


#endif 
