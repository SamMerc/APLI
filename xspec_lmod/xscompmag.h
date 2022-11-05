
typedef double Real;

#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

double **dmatrix(long nrl, long nrh, long ncl, long nch);
double *dvector(long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dvector(double *v, long nl, long nh);

double bb2d (double R, double x);
double bremss_abs(double x);

double slope(double energies[], double** spectrum, int qmin, double h, double ktbb, double kte, int NQ);
double check_conv(double previous_index,double index);
double interp_funct(double *array_x, double *array_y, int N, double x);
double gammln(double xx);
double zeta(double xx);
double h_funct(double xc);
double delta(double a, double xc);
/* ################################################################### */

void xscompmag2(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux);

__END_DECLS

