#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <math.h>
#include "demo_fn.h"

// gcc -lgsl -lgslcblas test.c && ./a.out


int main (void)
{
	double ka = 150.0;
	double m = 0.184;
	double z = 1.0;
	double c1 = 2.0;
	double ksi, n;    
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fdfsolver_type *T;
	gsl_root_fdfsolver *s;
	double x0, x = 5.0, r_expected = sqrt (5.0);
	gsl_function_fdf FDF;
	struct dukhin_params params = {ka, m, z, c1};
	FDF.f = &dukhin;
	FDF.df = &dukhin_deriv;
	FDF.fdf = &dukhin_fdf;
	FDF.params = &params;
	T = gsl_root_fdfsolver_newton;
	s = gsl_root_fdfsolver_alloc (T);
	gsl_root_fdfsolver_set (s, &FDF, x);
	printf ("using %s method\n",
			gsl_root_fdfsolver_name (s));
	printf ("%-5s %10s %10s %10s\n",
			"iter", "root", "err", "err(est)");
	do
	{
		iter++;
		status = gsl_root_fdfsolver_iterate (s);
		x0 = x;
		x = gsl_root_fdfsolver_root (s);
		status = gsl_root_test_delta (x, x0, 0, 1e-3);
		if (status == GSL_SUCCESS)
			printf ("Converged:\n");
		printf ("%5d %.30lf %+10.7f %10.7f\n",
				iter, x, x - r_expected, x - x0);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	printf("y^ek = %.25lf \n", x);


        double x1 = 9.9997;
        double zx4 = z*x1/4.0;
	double zx2 = z*x1/2.0;

	double sh2 = sinh(zx2);
	double sh4 = sinh(zx4);
	double ch2 = cosh(zx2);
	double ch4 = cosh(zx4);
        printf("Em=%lf \n", (ka+8.0*(1.0+(3.0*m)/)*sh4*sh4 - (24.0*m)*log(ch4))*x1 - c1*(ka+8.0*(1+(3.0*m)/(z*z))*sh4*sh4 - (24*m/(z*z))*log(ch4)) - x1*(1+3.0*m/(z*z))*sh4 + ((2.0/z)*sh2-3.0*m*x)*log(ch4));
        
	return status;

}
