#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
 #include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

int main (void)
{
	gsl_complex qq=gsl_complex_rect(0, 0);
	printf("GSL_COSH=%f\n",gsl_complex_cosh(qq));
	gsl_complex a=gsl_complex_cosh(qq);
	double dva=2;
	double b=GSL_REAL(a)*dva;
	 printf("GSL_COSH=%f\n",b);
	double result=GSL_REAL(gsl_complex_cosh(gsl_complex_rect(3.0/2.0, 0)));
	printf("GSL_RES=%f\n",result);
	printf("delenie= %f\n", log(2.7));
	return 0;


}
