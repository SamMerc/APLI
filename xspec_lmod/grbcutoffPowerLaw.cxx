#include <cfloat>
#include <xsTypes.h>

// prototype for routine which does the work

Real calcCutoffPowerLaw(const RealArray& energyArray, const Real& photIdx, const Real& cutoff, bool isRenorm, RealArray& fluxArray);


#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif


__BEGIN_DECLS
void
grbcutoffPowerLaw (const RealArray& energyArray,
                const RealArray& params,
                int spectrumNumber,
                RealArray& fluxArray,
                RealArray& fluxErrArray,
                const string& initString);
__END_DECLS

void
grbcutoffPowerLaw (const RealArray& energyArray,
                const RealArray& params,
                int spectrumNumber,
                RealArray& fluxArray,
                RealArray& fluxErrArray,
                const string& initString)
{
   // Power law with high energy exponential cutoff using E_peak.  
   // Number of model parameters: 2
   //   1       photIdx         powerlaw photon index
   //   2       e_peak          peak energy (in
   //                           energy units, e.g. keV). if <= 0
   //                           then no cutoff applied.
   // Intrinsic energy range:
   //   Emin = epsilon(>0), Emax = infinity
   //
   // algorithm:
   //   n(E)= E**(-photIdx) * exp(-E/E_peak * (2-photIdx)) dE
   //   This relies on an approximate incomplete gamma function
   //   calculation to perform the integral over n(E).  
   //   WARNING: The approximation loses accuracy in the region,
   //   10^-6 > (1-photIdx) > 0.

   const Real& photIdx = params[0];
   Real& w_cutoff = (Real &)params[1];
   if (photIdx != 2)
       w_cutoff /= (2.-photIdx);
   else
       w_cutoff = FLT_MAX;

   const Real& cutoff = w_cutoff;

   calcCutoffPowerLaw(energyArray, photIdx, cutoff, true, fluxArray);
   fluxErrArray.resize(0);

}

