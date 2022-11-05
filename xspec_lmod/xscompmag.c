/* ################################################################################### */
/* Numerical routine for the solution of the radiative transfer equation */
/* in the Fokker-Planck approximation for cylindrical simmetry */
/* using a finite-difference progonka method for relaxation problems */
/* Both the energy and space part of the equation are considered */
/* #################################################################################### */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>

#define NR_END 1
#define FREE_ARG char*
#define v_light 3*1e10
#define e_kin_fact 9.1094*1e-28* pow(v_light,2)*6.24*1e8
#define pi 3.14159
#define sigma 6.65*1e-25
#define M 15000
#define NSmass 1.4
#define PCradius  1
#define bulkflag 2
#define PI 3.14159
#define NTAU 40
#define NQ 150

double **dmatrix(long nrl, long nrh, long ncl, long nch);
double *dvector(long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dvector(double *v, long nl, long nh);

double bb2d (double R, double x);
static double bessi0( double x);
static double bessi1( double x);
static double bessk0( double x);
static double bessk1( double x);



void norm(double energies[], double** spectrum_up, double **spec_int, int qmin, double h, double htau, double kte, double**a, double**b, double** c, double** d, double** e, double** f, double series[]);


double interp_funct(double *array_x, double *array_y, int N, double x);
double gammln(double xx);
double zeta(double xx);
double h_funct(double xc);
double delta(double a, double xc);
double abs_val (double num);

/* ################################################################### */



double * xscompmag2(double* energy, int Nflux, double* parameter, int spectrum, double*  flux, double* fluxError, char* init) {

  int static call=0;

  double h, htau,  delta_ene, cmod, qmin, qmax, taumin, taumax, R0, r0,  norm_up, beta_max;
  double ktbb, kte, mdot, ratio_temp, theta, delta_hat, A, albedo,  deriv_vel=0, b12, dzdtau;
  double W, Z, P,Q, R_hat, G, H, csi, eta, zkm, zmax, aa, bremssabs, bremssabs_bulk, delta_val, vel_tau;
  double index=1.0, integral, dzsum, ntaumin, emin, emax;
  double p1, p2, p3, gaunt_fact, w1, w2, sigma_cyc, gauss, cyc_line, limit, ht, max_diff;
  

  double series[2];

  int bbflag, cycflag, bremssflag, geomflag, i_left_cyc=-1, i_right_cyc=-1;

/* Neutron star radius in Schwarzschild units*/
  double  z0=2.42;  

 /* Constant for the z(tau) formula for velocity profile v(z)=A*z^(-eta) */
  double k_zt=0.00219; 

  int i, j, ie, ii, m, iter_val=0, betaflag;
  

  FILE *dat;
 
/* ##########################################################*/
/* Model parameters */
/* ##########################################################*/

  
  ktbb=parameter[0];
  bbflag=parameter[1];
  kte=parameter[2];
  mdot=parameter[3];
  eta=parameter[4];
  beta_max=parameter[5];
  R0=parameter[6];
  zkm=parameter[7];
  b12=parameter[8];
  sigma_cyc=parameter[9];
  cycflag=parameter[10];
  albedo=parameter[11];
  betaflag=parameter[12];
  bremssflag=parameter[13];
  geomflag=parameter[14];
 

  /* printf("kTe=%5.3f mdot=%5.3f r0=%5.3f zmax=%5.3f sigma=%5.3f\n", kte, mdot, R0, zkm, sigma_cyc); */

/* #############################################################*/

  limit=0.001;

/* #############################################################*/
/* Convert R0 and zmax to adimensional units */

r0=R0/(2.95*NSmass);
  
zmax=zkm/(2.95*NSmass)+2.42;
 

/* #############################################################*/



H=(kte/511.)*(0.1/1e-3);
if (betaflag==1) {
taumin=0;
 } else {
  taumin=1e-4;
 }


/* zmax=2*z0; */
A=beta_max*pow(z0,eta); 

/* #############################################################*/
/*Determine mdot from the continuity equation for both*/
/*velocity profiles*/
/*#############################################################*/

 if (betaflag==1) {

/* mdot=taumax*A*pow(r0,2)*(1+eta)/(0.002174*(pow(zmax,1+eta)-pow(z0,1+eta))); */


taumax=mdot*0.001217*(pow(zmax,1+eta)-pow(z0,1+eta))/(A*pow(r0,2)*(1+eta));



 } else {


/* mdot=50*pow(r0,3./2)/(pow(z0,1./2) *pow(zmax-z0, 1./2))*taumax; */

taumax=mdot*(pow(z0,1./2) *pow(zmax-z0, 1./2))/(50*pow(r0,3./2));



}

/* ############################################################*/
/* Define logaritmic step of the adimensional energy */
/* ############################################################*/

 emin=0.1;
 emax=300; 
  /* qmin=-4; */
  /* qmax=6; */
 
 qmin=log(emin/kte);
 qmax=log(emax/kte);

  h=(qmax-qmin)/(NQ);

 /* ######################################################## */

  ratio_temp=kte/ktbb;
  theta=kte/511.;
 
/* ########################################################################## */

  G=(3/2.)*(1-albedo)/(1+albedo); 

/* #################################################################################### */
/* When using the definition J(x,tau)=R(tau) x^{-alpha} the inner boundary condition*/
/* becomes u[i][0]=L_tilde[i][0]* u[i][1] where  */
/* L_tilde[i][0]=1/(1+htau(-beta_max*(alpha+3)+G)); */
/* The necessary condition is L_tilde[i][0] > 0 */
/* #################################################################################### */
 
 /* if (-beta_max*(index+3)+G > 0) { */
/*    htau=h/sqrt(H);  */
/*    htau=0.5/(-beta_max*(index+3)+G);   */
/*  } else { */
/* htau=-0.5/(-beta_max*(index+3)+G); */
/*  } */


/* NTAU=(taumax-taumin)/htau;  */

/* #################################################################### */

ntaumin=40;  
htau=(taumax-taumin)/(ntaumin-1);


/* #################################################################### */

  double** a=dmatrix(0,NQ,0,NTAU);
  double** b=dmatrix(0,NQ,0,NTAU);
  double** c=dmatrix(0,NQ,0,NTAU);
  double** d=dmatrix(0,NQ,0,NTAU);
  double** e=dmatrix(0,NQ,0,NTAU);
  double** f=dmatrix(0,NQ,0,NTAU);
  double** s=dmatrix(0,NQ,0,NTAU);
  double** S_hat=dmatrix(0,NQ,0,NTAU);
  double** L=dmatrix(0,NQ,0,NTAU);
  double** K=dmatrix(0,NQ,0,NTAU);
  double** L_tilde=dmatrix(0,NQ,0,NTAU);
  double** K_tilde=dmatrix(0,NQ,0,NTAU);
  double** S_tilde=dmatrix(0,NQ,0,NTAU);
  double** spec_low=dmatrix(0,NQ,0,NTAU);
  double** spec_int=dmatrix(0,NQ,0,NTAU);
  double** spec_up=dmatrix(0,NQ,0,NTAU);
  double** flux_bulk=dmatrix(0,NQ,0,NTAU);

  /* double a[NQ+1][NTAU+1]; */
  /* double b[NQ+1][NTAU+1]; */
  /* double c[NQ+1][NTAU+1]; */
  /* double d[NQ+1][NTAU+1]; */
  /* double e[NQ+1][NTAU+1]; */
  /* double f[NQ+1][NTAU+1]; */
  /* double s[NQ+1][NTAU+1]; */
  /* double S_hat[NQ+1][NTAU+1]; */
  /* double L[NQ+1][NTAU+1]; */
  /* double K[NQ+1][NTAU+1]; */
  /* double L_tilde[NQ+1][NTAU+1]; */
  /* double K_tilde[NQ+1][NTAU+1]; */
  /* double S_tilde[NQ+1][NTAU+1]; */
  /* double spec_low[NQ+1][NTAU+1]; */
  /* double spec_int[NQ+1][NTAU+1]; */
  /* double spec_up[NQ+1][NTAU+1]; */
  /* double flux_bulk[NQ+1][NTAU+1]; */

  static double spec_low_before[NQ][NTAU];
  static double spec_int_before[NQ][NTAU];
  static double spec_up_before[NQ][NTAU];

 
   /* ################################################################## */

  double* flux_column=dvector(0,NQ); 
 
  for (i=0;i<=NQ;i++) {
    flux_column[i]=0;
  }



  double* flux_comp=dvector(0,Nflux); 
  double* x=dvector(0,Nflux); 
  double* flux_bb=dvector(0,Nflux);
  double* q=dvector(0,NQ);
  double* xnum=dvector(0,NQ); 
  double* flux_tot=dvector(0,NQ);
  double* vel=dvector(0,NTAU);
  double* ne=dvector(0,NTAU);
  double* tau=dvector(0,NTAU);
  double* ztau=dvector(0,NTAU);
  double* bfield=dvector(0,NTAU);
  double* xcyc=dvector(0,NTAU);
  double* phtot=dvector(0,NTAU);
  double* bbnorm=dvector(0,NTAU);
  double* bremssnorm=dvector(0,NTAU);
  double* cycnorm=dvector(0,NTAU);
  double* integ_int=dvector(0,NTAU);
  double* integ_low=dvector(0,NTAU);
  double* integ_up=dvector(0,NTAU);
  double* fact=dvector(0,NTAU);



/* ####################################################################################### */
/* NORMALIZATION OF THE MODEL */

/* The adimensional intensity for which the problem is j_ad=x^3*n */
/* To pass to physical units j_phys=2 kTe^3/(h^2 c^2) */
/* If the temperature is expressed in kev the expression becomes*/
/* j_phys=3.14*10^(31) (kTe/keV)^3 j_ad keV cm^-2 s^-1 keV^-1 ster^-1 */

/* The normalization, in similar way as the bbrad model of XSPEC is written as */
/* proportional to the emitting area divided the distance power two */
/* with the further factor PI to pass from intensity to flux for uniform brightness */
/* Here however instead of R(Km)=10^5 cm the radius is R0=2.95*10^5*m*r0 */



 if (geomflag==1) {

cmod=2.26*pow(10,-3)*pow(NSmass,2)*pow(r0,2)*pow(kte,3);

  } else if  (geomflag==2) {


cmod=4.52*pow(10,-3)*pow(NSmass,2)*r0*pow(kte,3);

  } else {

   cmod=pow(10,-3)*pow(NSmass,2)*pow(kte,3)*(4.52*r0 + 2.26*pow(r0,2)); 

  }


/* ############################################################# */
/* Define the step-size for the time variable */
/* ############################################################# */

 /* ht=0.1*3*H*htau*htau;  */
  ht=1/2. * htau*h;  
 /* ht=1e-3; */

/* printf("Value of H: %5.3e  1/H %5.3f\n", H, 1/H); */

 
csi=(15.77*r0)/mdot;
aa=0.67/z0*csi;

for (j=0; j<=NTAU; j++) {

tau[j]=taumin+j*htau;


 /* printf("j=%d tau=%5.3e  taumax=%5.3e (taumax-taumin)/htau=%5.3f \n", j, tau[j], taumax, taumax/htau);  */ 

/* ###################################################################### */
/* First velocity profile */
/* v(z)=-A*z^{-eta} */
/* ###################################################################### */

 if (betaflag==1) {

ztau[j]=pow(pow(z0,1+eta)+tau[j]*A*pow(r0,2)*(1+eta)/(mdot*k_zt),1/(1.+eta));

vel[j]=-A* pow(ztau[j],-eta);
/* deriv_vel=pow(A,2)*pow(r0,2)*eta*pow(ztau,-1+1./(1+eta))*pow(ztau,-1)/(mdot*k_zt); */

deriv_vel=pow(A,2)*pow(r0,2)*eta/(mdot*k_zt)*(pow(pow(z0,1+eta)+A*pow(r0,2)*(1+eta)*tau[j]/(mdot*k_zt),-2+1/(1+eta)));



/* ###################################################################### */
/* Define the electron density from the continuity equation */
/* ###################################################################### */

 ne[j]=(4.9*pow(10,21)*A*(1+eta)*taumax)/(pow(z0,1+eta)- pow(zmax,1+eta)*vel[j]);
 /* printf("tau=%5.3f ne=%5.3e\n", tau[j], ne[j]);     */

/* ###################################################################### */
/* Second velocity profile */
/* v(z)=-a*tau */
/* ###################################################################### */

 } else {

   ztau[j]=z0+2500*pow(r0,3)*pow(tau[j],2)/(pow(mdot,2)*z0);

  /* printf ("tau %5.3e z %5.3f mdot %5.3e\n", tau[j], ztau[j], mdot);   */

 aa=16.*sqrt(3.)/(49.*log(7./3.))*(1/z0)*csi;
 vel[j]= -aa*tau[j];
 deriv_vel=-aa;
 
 ne[j]=(1.3*pow(10,23)*pow(r0,3/2.)*pow(taumax,3))/(NSmass*pow(z0,1/2.)*pow(zmax-z0,3/2.) *tau[j]); 





 }

bfield[j]=b12*pow(z0,3)*pow(ztau[j],-3);
xcyc[j]=11.67*bfield[j]/kte;
 /* printf("j %d tau=%5.3f B=%5.3e xcyc[j]=%5.3e  z=%5.3e\n", j, tau[j], bfield[j], xcyc[j], ztau[j]);   */
 }


 

/* ########################################################################################### */
/* The Gaunt factor can be written as g(x)=Sqrt[3]/PI*Exp[x/2]*K0(x/2) */
/* where K0 is the modified Bessel function of the second kind */
/* The function g(x) can be fitted with accuracy within 5% over the range 10^(-3) < x < 50 */
/* by the function f(x)=a*Exp[-b*x**c] */
/* ########################################################################################### */


 p1=9.66588;
 p2=2.42165;
 p3=0.157505;


/* ########################################################################## */
/* Seed photon space and energy distrbution*/
/* ########################################################################### */
 
 for (j=0; j<=NTAU; j++) {
 
/* ###########################################################################*/
/* Exponential space distribution of the seed photons */
/* ###########################################################################*/

/* bbnorm[j]=exp(-tau[j])/(H); */

 bbnorm[j]=exp(-(ztau[j]-z0)/(1e-2*z0))/(H);


bremssnorm[j]=3.81*pow(10,-13)*ne[j]/(1e19)*pow(theta, -9/2.);
cycnorm[j]=1.67*pow(10,-8)*ne[j]/(1e19)*pow(bfield[j],-3/2.)*pow(theta,-4);

 /* printf("tau[j]  %5.3f bremmsnorm  %5.3e cycnorm %5.3e B12 %5.3e xcyc %5.3e\n", tau[j], bremssnorm[j], cycnorm[j], bfield[j], xcyc[j]); */    

 
for (i=0; i<=NQ; i++) {

q[i]=qmin+i*h; 

if    (i < NQ  &&  xcyc[j] < exp(qmin+(i+1)*h) && xcyc[j] > exp(q[i]))  { 

i_left_cyc=i;
i_right_cyc=i+1;

w1=1-(xcyc[j]-exp(q[i]))/(exp(qmin+(i+1)*h)- exp(q[i]));
w2=1-(exp(qmin+(i+1)*h)-xcyc[j])/(exp(qmin+(i+1)*h)- exp(q[i]));

/* printf(" j %d xcyc=%4.2e Ecyc=%4.2e   exp(q[i-1])%4.2e exp(q[i]) %4.2e B12 %4.2e w1 %4.2f w2 %4.2f\n\n", j, xcyc[j], xcyc[j]*kte,  exp(q[i]),  exp(qmin+(i+1)*h), bfield[j], w1, w2);     */
delta_val=1;

 } else {
delta_val=0;
 }


/* gaunt_fact=p1*exp(-p2*pow(exp(q[i]),p3));  */

 gaunt_fact=sqrt(3.)/PI*exp(exp(q[i])/2)*bessk0(exp(q[i])/2);

 if (i==i_left_cyc) {
   cyc_line=cycflag*w1*cycnorm[j]*h_funct(xcyc[j])*exp(-xcyc[j])*exp(q[i]);
  } else if (i==i_right_cyc) {
   cyc_line=cycflag*w2*cycnorm[j]*h_funct(xcyc[j])*exp(-xcyc[j])*exp(q[i]);
 } else {
   cyc_line=0;
 }

 
gauss=sqrt(1/(2*pi))*(1/sigma_cyc)*exp(-pow(exp(q[i])- xcyc[j],2)/(2*pow(sigma_cyc,2)));  

cyc_line=cycflag*cycnorm[j]*h_funct(xcyc[j])*exp(-xcyc[j])*exp(q[i])*gauss; 

s[i][j]=bremssflag*bremssnorm[j]*gaunt_fact*exp(-exp(q[i]))+cyc_line +bbflag*bbnorm[j]*bb2d(ratio_temp, exp(q[i])); 


 }

 }
  
/* ############################################################################## */
/* Build the initial guess of the spectrum equal to the seed spectrum */
/* ############################################################################### */

for (j=0; j<NTAU; j++) {
for (i=0; i<=NQ; i++) {

if (call < 5) {
spec_low[i][j]=s[i][j];
spec_int[i][j]=s[i][j];
spec_up[i][j]=s[i][j];
 
} else {

spec_low[i][j]=spec_low_before[i][j];
spec_int[i][j]=spec_int_before[i][j];
spec_up[i][j]=spec_up_before[i][j];

 }

 }

 }

/* ############################################################################### */
/* At the top of the column j=NTAU the spectrum is equal to zero */
/* ############################################################################### */

for (i=0; i<=NQ; i++) {

s[i][NTAU]=0;   

spec_low[i][NTAU]=0;
spec_int[i][NTAU]=0;
spec_up[i][NTAU]=0;

}


/* ############################################################################ */
/* Building the solution of the first equation of 2D progonka method */
/* For any y=j, we calculate the solution of the equation for the Lx-operator */
/* ############################################################################ */

for (j=0; j<=NTAU; j++) {


delta_hat=1./(3*H)*deriv_vel;

/* ########################################################################## */
/* Functions of the space operator */
/* ########################################################################## */

 
for (i=0; i<= NQ; i++) {


xnum[i]=exp(q[i]);


 bremssabs=2.1*1e-13* exp(-xnum[i]/2)*bessk0(xnum[i]/2)*ne[j]/(1e19)*(-1+exp(xnum[i]))/pow(theta,4.5);

 bremssabs_bulk =1/3. *ne[j]/(1e19)*vel[j]/pow(theta,4.5)* (bessk0(xnum[i]/2)* (2.10351*1e-13* xnum[i]* cosh(xnum[i]/2) - 8.41405*1e-13 *sinh(xnum[i]/2)) - 2.10351*1e-13* xnum[i]* bessk1(xnum[i]/2)* sinh(xnum[i]/2));

 /* printf(" bs1 %5.3e bs2 %5.3e\n", bremssabs,bremssabs_bulk); */

/* ########################################################################## */
/* Coefficients of the space operator */
/* ######################################################################### */


W=1./(3*H);
Z=-vel[j]/H+ bremssabs_bulk * exp(-3*q[i]);

/* ########################################################################## */
/* Coefficients of the energy operator */
/* ######################################################################### */

P= (3*kte+ (bulkflag-1)*e_kin_fact*pow(vel[j],2))/(3*kte);

Q= bremssabs_bulk*vel[j]*exp(-3*q[i]) + exp(q[i])-3+delta_hat - (bulkflag-1)*e_kin_fact *pow(vel[j],2)/kte;
 
R_hat= exp(q[i])-3*delta_hat-pow(vel[j],2)*pow(csi,2)/H - bremssabs*exp(-3*q[i])-3* vel[j]*bremssabs_bulk*exp(-3*q[i])  -1./ht;
 
/* First term of the bremsstrahlung absorption  */


 /* printf("Bremms %5.4e\n", -3* vel[j]*bremssabs_bulk*exp(-3*q[i])); */
 /* printf("CIAO P=%lf Q=%lf R=%lf W=%lf Z=%lf \n",P, Q, R_hat, W,Z); */

 

 a[i][j]=P/(h*h);
 b[i][j]=-2*P/(h*h)-Q/h+R_hat;
 c[i][j]=P/(h*h)+ Q/h;

  if (j==0) {
  
     d[i][j]=0;
     e[i][j]=-2*W/(htau*htau) - Z/htau;
     f[i][j]=W/(htau*htau)+ Z/htau;
  
   } else if (j==NTAU) {

      d[i][j]=W/(htau*htau);
      e[i][j]=-2*W/(htau*htau)- Z/htau;
      f[i][j]=0;
  
 } else {


     d[i][j]=W/(htau*htau);
     e[i][j]=-2*W/(htau*htau)- Z/htau;
     f[i][j]=W/(htau*htau)+ Z/htau;
  
  }

 }
 }


/* ########################################################################### */
/* Boundary conditions on energy (i=0,i=N) for coefficients and L,K                     */
/* ########################################################################### */
  
  for (j=0; j<= NTAU; j++) {
  
     a[0][j]=0; b[0][j]=1; c[0][j]=0; s[0][j]=0;
     a[NQ][j]=0; b[NQ][j]=1; c[NQ][j]=0; s[NQ][j]=0;
   
  L[0][j]=-c[0][j]/(b[0][j]);
  K[0][j]=0;

  }

/* ########################################################################### */
/* Starting the loop over m: at any cicle, we calculate the spectrum at three */
/* layers respect to the time (=m) */
/* ########################################################################### */

  m=0;  

 /* dat=fopen("spectau.dat","w");    */

 while (m <= M){  


/* ########################################################################## */
/* Recursively building of the coefficients L, K for any i,j                  */
/* ########################################################################## */

   for (j=0; j<= NTAU; j++) {

/* ########################################################################## */
/* Set to zero at the outer energy-boundary the right-hand side */
/* source term when solving for the energy operator  */
/* ########################################################################## */

     S_hat[0][j]=0;
     S_hat[NQ][j]=0;

/* ########################################################################## */
/* Build the source function on the right-hand side of the energy operator */
/* ########################################################################## */

   for (i=1; i< NQ; i++) {

   if (j==0) {

   S_hat[i][j]=-((e[i][j]+1./ht)*spec_low[i][j]+f[i][j]*spec_low[i][j+1]+s[i][j]);
      
  } else if (j==NTAU) {

     S_hat[i][j]=-(d[i][j]*spec_low[i][j-1]+(e[i][j]+1./ht)*spec_low[i][j]+s[i][j]);
      
  } else {
  
     S_hat[i][j]=-(d[i][j]*spec_low[i][j-1]+(e[i][j]+1./ht)*spec_low[i][j]+f[i][j]*spec_low[i][j+1]+s[i][j]);

   }

/* ########################################################################## */

   L[i][j]=-c[i][j]/(a[i][j]* L[i-1][j] + b[i][j]);
   K[i][j]=(S_hat[i][j]- a[i][j]*K[i-1][j])/(a[i][j]* L[i-1][j] + b[i][j]); 


   

 
   }
   
   }

/* ########################################################################## */
/* Find the solution of the progonka first equation and normalization      */ 
/* ########################################################################## */

/* printf("Iteration m=%d L[%d][%d]=%5.4e  K[%d][%d]=%5.4e\n", m, i,j, L[i][j], i, j, K[i][j]); */

for (j=0; j<=NTAU; j++) {

spec_int[NQ][j]=0; 
 
integ_int[j]=0;
integ_low[j]=0;
integ_up[j]=0;

/* ########################################################################## */

for (i=NQ; i>0; i--) {
spec_int[i-1][j]=L[i-1][j]*spec_int[i][j]+K[i-1][j]; 

}

}

/* ########################################################################## */
/* Boundary conditions on the optical depth for the coefficients and L,K */
/* of the second equation */
/* ########################################################################## */

for (i=0; i<= NQ; i++) {

  if (betaflag==1) {
  L_tilde[i][0]=1./(1+htau*(-beta_max*(index+3)+G));
  } else {
 L_tilde[i][0]=1./(1+htau*G);

  }

 K_tilde[i][0]=0;
 }

/* ########################################################################## */
/* Recursively building of the coefficients L_tilde, K_tilde  */
/* for the second progonka over tau */
/* ########################################################################## */
 
  for (i=1; i< NQ; i++) {
  for (j=1; j< NTAU; j++) {
 
  
 L_tilde[i][j]=-f[i][j]/(d[i][j]* L_tilde[i][j-1] + e[i][j]-1./ht);

  if (j==NTAU) {

	 S_tilde[i][j]=(d[i][j]*spec_low[i][j-1]+e[i][j]*spec_low[i][j]-(1./ht)*spec_int[i][j]);

    } else {

	 S_tilde[i][j]=(d[i][j]*spec_low[i][j-1]+e[i][j]*spec_low[i][j]+f[i][j]*spec_low[i][j+1]-(1./ht)*spec_int[i][j]);

       }

   K_tilde[i][j]=(S_tilde[i][j]- d[i][j]*K_tilde[i][j-1])/(d[i][j]* L_tilde[i][j-1] + e[i][j]-1./ht);

   }
   
}  

/* ########################################################################## */
/* Finding the solution of the progonka first equation and normalization      */ 
/* ########################################################################## */
  
for (j=0; j<=NTAU; j++) {
    spec_up[0][j]=0;  
    spec_up[NQ][j]=0; 
    integ_up[j]=0;
      }  

for (i=0; i<=NQ; i++) {
 
if (m==0)  spec_up[i][NTAU]=s[i][NTAU];
if (m > 0)  spec_up[i][NTAU]=spec_low[i][NTAU];


  }


/* ################################################################################ */

for (j=NTAU; j>0; j--) {
for (i=1; i< NQ; i++) {

 spec_up[i][j-1]=L_tilde[i][j-1]*spec_up[i][j]+K_tilde[i][j-1]; 
   if (spec_up[i][j-1] > 0 && spec_up[i][j-1] < 1e-30) spec_up[i][j-1]=0;

     /*  integ_low[j]=integ_low[j]+spec_low[i][j-1]*h; */
 /*      integ_int[j]=integ_int[j]+spec_int[i][j-1]*h; */
 /*      integ_up[j]=integ_up[j]+spec_up[i][j-1]*h; */
 }

 }
 


norm(x, spec_up, spec_low, qmin, h, htau, kte, a,b,c,d,e,f, series);



 max_diff=series[0];
 integral=series[1];

 /* double max_diff = norm(x, spec_up, spec_low, qmin, h, htau, kte, NQ, NTAU,a,b,c,d,e,f); */

 if (m > 0) {

   /* fprintf(dat, "m %d  %5.3e %5.3e\n", m, max_diff, integral);   */


  if (max_diff <  limit*ht) {
   /* printf("Iteration stopped at m=%d (ht=%7.4e)   (qmax=%5.3f) (htau=%5.4e)  (h=%7.5f)\n", m, ht, qmax, htau, h); */
      break;  

   }

 }

/* ######################################################################################### */


for (j=0; j< NTAU; j++) {

vel_tau=vel[j];

for (i=0; i<=NQ; i++) {
 
if (m < M)  spec_low[i][j]=spec_up[i][j];


/* ####################################################################################### */
/* Now write the expression for the flux at j=NTAU-1*/
/* ####################################################################################### */
  
if (j==NTAU-1)  {

if (i==NQ) {
flux_bulk[i][j]=-(spec_up[i][j]-spec_up[i][j-1])/(3*htau) + vel_tau*(spec_up[i][j]-(1./3)*(spec_up[i][j]-spec_up[i-1][j])/h);
 } else {
flux_bulk[i][j]=-(spec_up[i][j]-spec_up[i][j-1])/(3*htau) + vel_tau*(spec_up[i][j]-(1./3)*(spec_up[i+1][j]-spec_up[i][j])/h);
 }

 }

/* End of loop over i */
 }

/* End of loop over j */
 }


 
iter_val=m;
 
 m++;

 }


dzsum=0;

for(i=0;i<=NQ;i++){

 
 for (j=0; j< NTAU; j++) {


    if (betaflag==1) {

      dzdtau=A*pow(r0,2)/ (mdot*k_zt)* pow(pow(z0,1+eta)+ A*pow(r0,2)*(1+eta)*tau[j]/(mdot*k_zt),-1+1/(1+eta));
 
    } else {

      dzdtau=5000*pow(r0,3)/(pow(mdot,2)*z0)*tau[j];

    }

  flux_column[i] += spec_up[i][j]*dzdtau*htau;
  

  /* printf("column %5.3e spec %5.3e dzdtau %lf  htau %lf\n", flux_column[i], spec_up[i][j], dzdtau, htau);  */

 
}
 }


 /* fclose(dat);   */

 
/* ############################################################################################################## */


for(i=0;i<=NQ;i++){

if (geomflag==1) {
flux_tot[i]=4*PI*flux_bulk[i][NTAU-1];

} else if  (geomflag==2) { 
flux_tot[i]=flux_column[i];



  } else {
flux_tot[i]=flux_column[i]+4*PI*flux_bulk[i][NTAU-1];

  }


 /* printf("i %d flux_tot[i] %7.4e\n", i, flux_tot[i]); */

}



/* ####################################################################################### */
/* if (norm_up < 0) printf ("Norm_up < 0 %7.3e\n", norm_up); */
/* ####################################################################################### */


for (ie = 0; ie <=Nflux; ie++) {

ii=ie-1;
x[ie]=(energy[ie])/kte;


/* ####################################################################################### */

if (x[ie] < xnum[0]) {

flux_comp[ie]=flux_tot[0];

} else if (x[ie] > xnum[NQ]) {

flux_comp[ie]=flux_tot[NQ];

} else {

flux_comp[ie]=interp_funct(xnum, flux_tot, NQ, x[ie]);

}
 
/* ################################################################ */

flux_comp[ie]=cmod*flux_comp[ie]/energy[ie]; 

if (ie > 0) {

delta_ene=(energy[ie]-energy[ie-1]);

flux[ii]=0.5*(flux_comp[ie]+flux_comp[ie-1])*delta_ene ;

 }

}

/* ##################################################################### */


for (i=0;i<=NQ; i++) {
 for (j=0;j<=NTAU; j++) {

   spec_low_before[i][j]=spec_low[i][j];
   spec_int_before[i][j]=spec_int[i][j];
   spec_up_before[i][j]=spec_up[i][j];

}

 }

/* ##################################################################### */

 free_dmatrix(a,0,NQ,0,NTAU);
 free_dmatrix(b,0,NQ,0,NTAU);
 free_dmatrix(c,0,NQ,0,NTAU);
 free_dmatrix(d,0,NQ,0,NTAU);
 free_dmatrix(e,0,NQ,0,NTAU);
 free_dmatrix(f,0,NQ,0,NTAU);
 free_dmatrix(s,0,NQ,0,NTAU);
 free_dmatrix(S_hat,0,NQ,0,NTAU);
 free_dmatrix(L_tilde,0,NQ,0,NTAU);
 free_dmatrix(K_tilde,0,NQ,0,NTAU);
 free_dmatrix(S_tilde,0,NQ,0,NTAU);
 free_dmatrix(L,0,NQ,0,NTAU);
 free_dmatrix(K,0,NQ,0,NTAU);
 free_dmatrix(spec_low,0,NQ,0,NTAU);
 free_dmatrix(spec_int,0,NQ,0,NTAU);
 free_dmatrix(spec_up,0,NQ,0,NTAU);
 free_dmatrix(flux_bulk,0,NQ,0,NTAU);

 free_dvector(xnum, 0,NQ); 
 free_dvector(flux_tot,0,NQ);
 free_dvector(flux_comp,0,Nflux);  
 free_dvector(vel, 0,NTAU);
 free_dvector(tau,0,NTAU);
 free_dvector(phtot,0,NTAU);
 free_dvector(bbnorm,0,NTAU);
 free_dvector(integ_int,0,NTAU);
 free_dvector(integ_low,0,NTAU);
 free_dvector(integ_up,0,NTAU);
 free_dvector(fact,0,NTAU);
 free_dvector(flux_bb,0,Nflux); 
 free_dvector(x,0,Nflux);  
 free_dvector(q, 0,NQ); 

 call++;

  return(flux);

}


/* #################################################################################### */
/* Blackbody seed energy spectrum */
/* #################################################################################### */

double bb2d(double R, double z) {

   double value;
   value=pow(z,3)/(exp(R*z)-1);

   if (value > 1e-30) {
     return(value);
   } else{
     return(0.);
   }

}



/* ############################################################################# */
/* Allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
/* ############################################################################# */

double **dmatrix(long nrl, long nrh, long ncl, long nch)

{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	/* check_alloc(m); */
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) {
        printf("Allocation failure 2 in matrix()");
	exit(1);
	}

	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

/* ###################################################################### */
/* Allocate a double vector with subscript range v[nl..nh] */
/* ############################################################################# */

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) {
        printf("allocation failure in dvector()");
	exit(1);
	}

	return v-nl+NR_END;
}

/* ###################################################################### */

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)

{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

/* ###################################################################### */

void free_dvector(double *v, long nl, long nh)

{
	free((FREE_ARG) (v+nl-NR_END));
}


/* ############################################################################# */

double interp_funct(double *array_x, double *array_y, int N, double x) {

  int i;
  double m, q, interp_y;

  for (i=0; i<N; i++) {

    if (array_x[i] <= x && array_x[i+1] >=x) {

	m=(array_y[i+1]-array_y[i])/(array_x[i+1]-array_x[i]);	
	q=array_y[i]-m*array_x[i];
	interp_y=m*x+q;
	return(interp_y);
     }

  }

  return(0.);

}

/* ##################################################################*/
/* Logarithm of the gamma function */
/* ##################################################################*/

double gammln(double xx)

{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

/* ##################################################################*/
/* Riemann zeta function */
/* ##################################################################*/


double zeta(double xx) {

  int n;
  double value=0;

 for (n=1; n<=10; n++) {
 value=value+1/(pow(n,xx));
 }

 return(value);

}

/* ##################################################################*/
/* H-function for the cyclotron emission (eq. 115 BW07) */
/* ##################################################################*/

double h_funct(double xc) {
double f1, xt,t, b,c;

  /* if (xc > 7.5) { */
  /* return(0.41); */
  /* } else { */
  /*   return(0.15*pow(xc,0.5)); */
  /* } */

xt=7;
if (xc > xt) {

  b=0.5;
  t=(2.73/sqrt(xt)-1)*exp(xt/b);
  c=b*log(t);

  f1=0.41/(1+exp(-(xt-c)/b));
 return(f1);

 } else {
 return(0.15*pow(xc,0.5));
 }

}


  

/* ##################################################################*/
/* Delta-function for the cyclotron emission (eq. 115 BW07) */
/* ##################################################################*/

double delta (double a, double xc) {

  if (xc > 7.5) {
  return(0.41);
  } else {
    return(0.15*pow(xc,0.5));
  }

}


/* ########################################################################## */
/* Chek integration of the spectrum over energy and tau */
/* ########################################################################## */

/* double norm (double energies[], double** spectrum_up, double** spectrum_low, int qmin, double h, double htau, double kte, int NQ, int NTAU) { */



void norm(double energies[], double** spectrum_up, double **spectrum_low, int qmin, double h, double htau, double kte, double**a, double**b, double** c, double** d, double** e, double** f, double series[]) {


  double max_diff=-1e5, local_diff;
  int i_maxdiff, j_maxdiff;

  int static niter=0;
  double static first_val;

  double emin,emax, integ=0, resid=0,  diff_op=0, diff_norm=0, norm_up=0;
  int Nmin, Nmax, i,j, ntot=0;

/* ################################################################### */
/* Start from energies > 7 ktbb */

 emin=0.1;
 emax=200;

 Nmin=(log(emin/kte)-qmin)/h;
 Nmax=(log(emax/kte)-qmin)/h;

for (j=0; j<=NTAU-1; j++) {
for (i=0; i<=NQ-1; i++) {

ntot++;

integ=integ+pow(spectrum_up[i][j],2)*htau*h;

/* ############################################################################################################### */
/* Make a point-by-point check of the difference of matrix */

local_diff=abs_val(spectrum_up[i][j]-spectrum_low[i][j])/spectrum_up[i][j];

if (local_diff > max_diff) {

max_diff=local_diff;

 i_maxdiff=i;
 j_maxdiff=j;

}

/* ############################################################################################################### */
 /* Consider the convergence of the series |V_{k+1} -V_{k}|/|V_{k+1}| */

 norm_up=norm_up+pow(spectrum_up[i][j],2);
 resid=resid+pow(spectrum_up[i][j]-spectrum_low[i][j],2);



/* ############################################################################################################### */
 /* resid=resid+ht*(a[i][j]+b[i][j]+c[i][j])*(d[i][j]+e[i][j]+f[i][j])*(spectrum_up[i][j]-spectrum_low[i][j]); */

/* #################################################################################### */

diff_op=diff_op+(a[i][j]+b[i][j]+c[i][j]+d[i][j]+e[i][j]+f[i][j])*spectrum_up[i][j];


diff_norm=diff_norm+pow(spectrum_up[i][j]-spectrum_low[i][j],2);



 /* printf("spec_low[%d][%d] %5.3e\n", i,j, spectrum_low[i][j]); */
/* printf("Nmin=%d Nmax=%d Integral: %5.3e\n", Nmin, Nmax, integ); */

 }
 }

 /* printf("Residual %5.3e\n", abs_val(resid/diff_op));  */
 
 if (niter==1) first_val=sqrt(diff_norm);
 

 /* printf("Niter %d first_val %5.3e  diff_norm %5.3e ratio %5.3e\n", niter, first_val, sqrt(diff_norm), sqrt(diff_norm)/first_val);  */

/* printf("Max difference %7.5e at i=%d j=%d   |V_{k+1} -V_{k}|/|V_{k+1}|=%7.5e\n", max_diff, i_maxdiff, j_maxdiff, sqrt(resid/norm_up));   */

 series[0]=max_diff;
 series[1]=sqrt(resid/norm_up);

 niter++;
 
/* return(max_diff);  */

 /* return(sqrt(integ)); */

}



double abs_val (double num) {


if (num>=0) {
 return(num);
 } else {
  return(-num);
 }



}



/* Bessel functions */

double bessel_k0 (double x) {

  double a,b,c,f;

a=4.23279;
b=2.3363;
c=0.616657;

 f=a*exp(-b*pow(x,c));

 return(f);


}


static double bessi0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=0.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   } else {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2))))))));
   }
   return ans;
}


static double bessi1( double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
   } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (exp(ax)/sqrt(ax));
   }
   return x < 0.0 ? -ans : ans;
}



static double bessk0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=0.  */
/*------------------------------------------------------------*/
{
   double y,ans;

   if (x <= 2.0) {
      y=x*x/4.0;
      ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
         +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
         +y*(0.10750e-3+y*0.74e-5))))));
   } else {
      y=2.0/x;
      ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
         +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
         +y*(-0.251540e-2+y*0.53208e-3))))));
   }
   return ans;
}




static double bessk1( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double y,ans;

   if (x <= 2.0) {
      y=x*x/4.0;
      ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
         +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
         +y*(-0.110404e-2+y*(-0.4686e-4)))))));
   } else {
      y=2.0/x;
      ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
         +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
         +y*(0.325614e-2+y*(-0.68245e-3)))))));
   }
   return ans;
}

