struct dukhin_params
{
    double ka, m, z, c1;
};

double dukhin (double x, void *params);
double dukhin_deriv (double x, void *params);
void dukhin_fdf (double x, void *params, double *y, double *dy);

