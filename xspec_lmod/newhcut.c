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

void newhcut(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux);
void c_newhcut(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux)
{
	newhcut(energy, Nflux, parameter, spectrum, flux);
	return;
}

__END_DECLS

void newhcut(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux)
{
	int i;
	double Ecut = parameter[0];
	double Efold = parameter[1];
	double W = parameter[2];
	double u, E, A, B;
	double aux1,aux2,aux3;
	//It implements eq. (2) in https://iopscience.iop.org/article/10.1086/320643/pdf

	aux1=W/2.0;
	aux2 = exp(-(aux1/Efold));
        aux3 = W/Efold;

        A = -((aux3+2.)*aux2)/(W*W*W)+(2./(W*W*W));
        B = ((3.+aux3)*aux2)/(W*W) - (3./(W*W));

	for(i=0;i<Nflux;i++)
	{
		u = ( energy[i] + energy[i+1] ) / 2.;
		E = u - Ecut + aux1;
		
		if(E <= 0)
			flux[i] = 1.0;
		else if( E >= W )
			flux[i] = exp((aux1-E)/Efold);
		else
			flux[i] = A*E*E*E + B*E*E + 1.0;
		
	}

}
