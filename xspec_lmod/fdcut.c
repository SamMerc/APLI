#ifndef g77Fortran
#define g77Fortran 1
#endif
#include <math.h>

typedef double Real;

#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif


__BEGIN_DECLS

void fdcut(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux);
inline void c_fdcut(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux)
{

	return fdcut(energy, Nflux, parameter, spectrum, flux);
}

__END_DECLS

void fdcut(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux)
{
	int i;
	double Ecut = parameter[0];
	double Efold = parameter[1];
	double u=0, E=1;


	for(i=0;i<Nflux;i++)
	{
		u = ( energy[i] + energy[i+1] ) / 2.;
		E = u - Ecut;
		
		flux[i] = 1./ (1.0 + exp((E)/Efold) );
	}

}
