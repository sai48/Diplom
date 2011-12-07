#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <math.h>
#include "demo_fn.h"

double dukhin (double x, void *params)
{
	struct dukhin_params *p
		= (struct dukhin_params *) params;
	double ka = p->ka;
	double m = p->m;
	double z = p->z;
	double c1 = p->c1;
	double zx4 = z*x/4.0;
	double zx2 = z*x/2.0;

	double sh2 = sinh(zx2);
	double sh4 = sinh(zx4);
	double ch2 = cosh(zx2);
	double ch4 = cosh(zx4);
	return ((ka+8.0*(1.0+(3.0*m)/(z*z))*sh4*sh4 - (24.0*m/(z*z))*log(ch4))*x - c1*(ka+8.0*(1+(3.0*m)/(z*z))*sh4*sh4 - (24*m/(z*z))*log(ch4)) - x*(1+3.0*m/(z*z))*sh4 + ((2.0/z)*sh2-3.0*m*x)*log(ch4));
}

double dukhin_deriv (double x, void *params)
{
	struct dukhin_params *p
		= (struct dukhin_params *) params;
	double ka = p->ka;
	double m = p->m;
	double z = p->z;
	double c1 = p->c1;
	double zx4 = z*x/4.0;
    double zx2 = z*x/2.0;

    double sh2 = sinh(zx2);
    double sh4 = sinh(zx4);
    double ch2 = cosh(zx2);
        //double ch4 = GSL_REAL(gsl_complex_cosh(gsl_complex_rect(zx4, 0)));
    double ch4 = cosh(zx4);
	return (ka + 8.0*(1.0 + (3.0*m)/(z*z))*sh4*sh4 - (24.0*m/(z*z))*log(ch4) + x*((2.0*z+6.0*m/z)*sh2 - (6.0*m*sh4)/(z*ch4)) - c1*((2.0*z+6.0*m/z)*sh2 - (6.0*m*sh4)/(z*ch4)) - (1.0+3.0*m/(z*z))*sh4 - x*ch4*(z/4.0+(3.0*m)/(4.0*z)) + log(ch4)*(ch2 - 3.0*m) + (sh4/ch4)*(0.5*sh2-(3.0*m*x*z)/4.0));
}
void dukhin_fdf (double x, void *params,
		double *y, double *dy)
{
	struct dukhin_params *p
		= (struct dukhin_params *) params;
	double ka = p->ka;
    double m = p->m;
    double z = p->z;
    double c1 = p->c1;
	double zx4 = z*x/4.0;
    double zx2 = z*x/2.0;

	double sh2 = sinh(zx2);
	double sh4 = sinh(zx4);
	double ch2 = cosh(zx2);
	double ch4 = cosh(zx4);
	*y = (ka+8.0*(1.0+(3.0*m)/(z*z))*sh4*sh4 - (24.0*m/(z*z))*log(ch4))*x - c1*(ka+8.0*(1+(3.0*m)/(z*z))*sh4*sh4 - (24*m/(z*z))*log(ch4)) - x*(1+3.0*m/(z*z))*sh4 + ((2.0/z)*sh2-3.0*m*x)*log(ch4);
	*dy = (ka + 8.0*(1.0 + (3.0*m)/(z*z))*sh4*sh4 - (24.0*m/(z*z))*log(ch4) + x*((2.0*z+6.0*m/z)*sh2 - (6.0*m*sh4)/(z*ch4)) - c1*((2.0*z+6.0*m/z)*sh2 - (6.0*m*sh4)/(z*ch4)) - (1.0+3.0*m/(z*z))*sh4 - x*ch4*(z/4.0+(3.0*m)/(4.0*z)) + log(ch4)*(ch2 - 3.0*m) + (sh4/ch4)*(0.5*sh2-(3.0*m*x*z)/4.0));
}
