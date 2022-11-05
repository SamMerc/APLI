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

void gausabs(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux);
void agausabs(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux);


void c_gausabs(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux)
{
	gausabs(energy, Nflux, parameter, spectrum, flux);
	return;
}
void c_agausabs(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux)
{
	agausabs(energy, Nflux, parameter, spectrum, flux);
	return;
}


__END_DECLS

void gausabs(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux)
{
	int i;
	double ec = parameter[0];
	double sc = parameter[1];
	double dc = parameter[2];
	double abs_c, e, de;
	for(i=0;i<Nflux;i++)
	{
		e = ( energy[i] + energy[i+1] ) / 2.;
		de = (e -ec)/sc;
		abs_c=1;
		if(fabs(de) <= 6)
			abs_c -= dc / sc / sqrt(2.*M_PI) *exp(- 0.5*de*de);
		
		flux[i] = abs_c;
	}

}

void agausabs(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux)
{
	int i;
	double ec = parameter[0];
	double sc1 = parameter[1];
	double sc2 = parameter[2];
	double dc = parameter[3];
	double abs_c, e, de;
	for(i=0;i<Nflux;i++)
	{
		e = ( energy[i] + energy[i+1] ) / 2.;
		de = (e -ec);
		if(de<=0)
			de/=sc1;
		else
			de/=sc2;
		abs_c=1;
		if(de <= 6 && de>0)
			abs_c -= dc / sc2 / sqrt(2.*M_PI) *exp(- 0.5*de*de);
		else if (de >= -6 && de<=0)
			abs_c -= dc / sc1 / sqrt(2.*M_PI) *exp(- 0.5*de*de);
			
		
		flux[i] = abs_c;
	}

}
