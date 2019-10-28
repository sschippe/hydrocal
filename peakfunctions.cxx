/** 
 * @file peakfunctions.cxx
 * @brief Coding of peak functions 
 * 																			
 * @author Stefan Schippers
 * @date 2013-04-23
 * @version $Id$
 *										            
*/

#include <complex>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include "hydromath.h"
#include "readxsec.h"
#include "Faddeeva_w.h"
#include "peakfunctions.h"


void sum_of_peaks(void)
{
  int read_mode=0, peak_shape = PEAK_undefined, npeaks=0, number_of_levels=1, level_number=1;
  char peak_filename[200];

  printf("\n Create a spectrum as a sum of peaks.\n");
  npeaks = open_peak_file(peak_filename, read_mode, number_of_levels);
  if (npeaks<0)
    {
      printf("\n File %s not found.\n",peak_filename);
      return;
    }
  printf("\n %d peak data sets detected\n",npeaks); 

  double *energy   = new double[npeaks];
  double *strength = new double[npeaks];
  double *qFano    = new double[npeaks];
  double *wLorentz = new double[npeaks];
  double *wGauss   = new double[npeaks];
  double emin, emax;

  peak_shape = read_peak_file(peak_filename,read_mode,level_number,energy,strength,qFano,wLorentz,wGauss,npeaks,emin,emax);
  if (peak_shape==PEAK_undefined) return;

  if (number_of_levels>1)
    {
      printf("\n Data for %d different initial levels found.",number_of_levels);
      printf("\n Give the number (in the range 1 - %d) of the level to be used: ",number_of_levels);
      scanf("%d",&level_number);
      if ( (level_number<1) || (level_number>number_of_levels) )
	{
	  printf("\n\n ERROR: Level number of of range!\n\n");
	  exit(0);
	}
    }            

  double wConvolution =0.0;
  if (peak_shape==PEAK_delta) 
    {
      printf("\n Convolute delta peaks with Gaussian?");
      printf("\n  (no convolution for zero width)");
      printf("\n Give Gaussian width ..........: ");
      scanf("%lf",&wConvolution);
      if (wConvolution>0.0)  peak_shape = PEAK_Gauss;
    }
  else if ((peak_shape==PEAK_Lorentz)||(peak_shape==PEAK_Lorentz_Steih))
    {
      printf("\n Convolute Lorentzian peaks with Gaussian?");
      printf("\n  (no convolution for zero width)");
      printf("\n Give Gaussian width ..........: ");
      scanf("%lf",&wConvolution);
      if (wConvolution>0.0)  peak_shape = PEAK_Voigt;
    }
  else if (peak_shape==PEAK_Fano)
    {
      printf("\n Convolute Fano peaks with Gaussian?");
      printf("\n  (no convolution for zero width)");
      printf("\n Give Gaussian width ..........: ");
      scanf("%lf",&wConvolution);
      if (wConvolution>0.0)  peak_shape = PEAK_FanoVoigt;
    }
  else
    {
      printf("\n Convolute peaks with Gaussian?");
      printf("\n  (no convolution for zero width)");
      printf("\n Give Gaussian width ..........: ");
      scanf("%lf",&wConvolution);
    }
  if (wConvolution > 0.0)
    {
      for (int i=0; i< npeaks; i++ ) 
	{
	  wGauss[i] = sqrt(wGauss[i]*wGauss[i]+wConvolution*wConvolution);
	  //printf("%12.4f %12.4g %12.4f %12.4f %12.4f\n",energy[i],strength[i],qFano[i],wLorentz[i],wGauss[i]);
	}
    }
  double xx, yy, xdelta, xmin = 9E99, xmax = -9E99;
  int nlines=0, n, itest;
  FILE *fout;
  char filenameout[200];

  printf("\n Peak positions range from %g to %g\n",emin,emax);
  printf("\n Give emin, emax, and edelta .............................: ");
  scanf("%lf %lf %lf",&xmin,&xmax,&xdelta);
  int npts = int(fabs(xmax-xmin)/xdelta+1.1);

  printf("\n Give name of output file ................................: ");
  scanf("%s",filenameout);
  fout = fopen(filenameout,"w");  

  if (peak_shape == PEAK_Gauss)
    {
      printf("\n Calculating sum of Gauss peaks ...");
      for(int n=0; n < npts; n++) 
	{
	  xx = xmin + n*xdelta; 
	  yy = 0.0;
	  for (int i=0; i<npeaks; i++)
	    {
	      yy += gauss(xx,energy[i],strength[i],wGauss[i]); 
	    }
	  fprintf(fout,"%20g %20g\n",xx,yy);
	}
    } 
  else if ((peak_shape == PEAK_Lorentz)||(peak_shape == PEAK_Lorentz_Steih))
    {
      printf("\n Calculating sum of Lorentz peaks ...");
      for(int n=0; n < npts; n++) 
	{
	  xx = xmin + n*xdelta; 
	  yy = 0.0;
	  for (int i=0; i<npeaks; i++)
	    {
	      yy += lorentz(xx,energy[i],strength[i],wLorentz[i]); 
	    }
	  fprintf(fout,"%20g %20g\n",xx,yy);
	}
    }
  else if (peak_shape == PEAK_Voigt)
    {
      printf("\n Calculating sum of Voigt peaks ...");
      for(int n=0; n < npts; n++) 
	{
	  xx = xmin + n*xdelta; 
	  yy = 0.0;
	  for (int i=0; i<npeaks; i++)
	    {
	      yy += voigt(xx,energy[i],strength[i],wLorentz[i],wGauss[i]); 
	    }
	  fprintf(fout,"%20g %20g\n",xx,yy);
	}
    }
  else if (peak_shape == PEAK_Fano)
    {
      printf("\n Calculating sum of Fano peaks ...");
      for(int n=0; n < npts; n++) 
	{
	  xx = xmin + n*xdelta; 
	  yy = 0.0;
	  for (int i=0; i<npeaks; i++)
	    {
	      yy += fano(xx,energy[i],strength[i],qFano[i],wLorentz[i]); 
	    }
	  fprintf(fout,"%20g %20g\n",xx,yy);
	}
    }
  else if (peak_shape == PEAK_FanoVoigt)
    {
      printf("\n Calculating sum of FanoVoigt peaks ...");
      for(int n=0; n < npts; n++) 
	{
	  xx = xmin + n*xdelta; 
	  yy = 0.0;
	  for (int i=0; i<npeaks; i++)
	    {
	      yy += fanovoigt(xx,energy[i],strength[i],qFano[i],wLorentz[i],wGauss[i]); 
	    }
	  fprintf(fout,"%20g %20g\n",xx,yy);
	}
    }

  fclose(fout);

  delete[] wGauss;
  delete[] wLorentz;
  delete[] qFano;
  delete[] strength;
  delete[] energy;
}


////////////////////////////////////////////////////
/**
 * @brief Area normalized gaussian
 *
 *	@param E  energy axis
 *  @param E0 peak center energy
 *  @param wg gaussian linewidths (FWHM)
 *	@param a  peak area
 *
 *  @return profile height at the given energy
 */
double gauss(double E, double E0, double a, double wg)
{
  const double sqrtln2 = 0.832554611; 
  const double sqrtpi  = 1.772453851;
  double w = 0.5*wg/sqrtln2;
  double eps = (E-E0)/w;
  if (fabs(eps)>10.0) return 0.0;
  double eps2 = eps*eps;
  
  return a/w*exp(-eps2)/sqrtpi; 
}

////////////////////////////////////////////////////
/**
 * @brief Area normalized lorentzian
 *
 *  @param E  energy axis
 *  @param E0 peak center energy
 *  @param a  peak area
 *  @param wl lorentzian linewidths (FWHM)
 *
 *  @return profile height at the given energy
 */
double lorentz(double E, double E0, double a, double wl)
{
  const double pi = 3.141592654;
  
  double eps = 2.0*(E-E0)/wl;

  return 2.0*a/(1.0+eps*eps)/wl/pi;
}

////////////////////////////////////////////////////
/**
 * @brief Fano profile
 *
 * @param E   energy axis
 * @param E0  peak center energy
 * @param a   peak area
 * @param q   esymmetry parameter
 * @param wl  lorentzian linewidths (FWHM)
 *
 * @return height of the profile at the given energy (goes to zero for E->+-infinity)
*/
double fano(double E, double E0, double a, double q, double wl)
{
  const double pi      = 3.141592654;
  
  double eps = 2.0*(E-E0)/wl;
  return 2.0/wl/pi*((q+eps)*(q+eps)/(1.0+eps*eps)-1.0)*fabs(a/(q*q-1.0));
}


/////////////////////////////////////////////////////////////////////
/**
 * @brief Voigt profile (Convolution of a lorentzian with a gaussian)
 *
 * @param E   energy axis
 * @param E0  peak center energy
 * @param a   peak area
 * @param wl  lorentzian linewidths (FWMH)
 * @param wg  gaussian linewidths (FWHM)
 *
 * @return height of the profile at the given energy (goes to zero for E->+-infinity)
*/
double voigt(double E, double E0, double a, double wl, double wg)
{
  
  if (fabs(wg)<1e-20)
    {
      return lorentz(E,E0,a,wl);
    }
  else
    {
      const double sqrtln2 = 0.832554611; 
      const double sqrtpi  = 1.772453851;
      double x, y, re, im;
      x = 2.0*sqrtln2*(E-E0)/wg;
      y = sqrtln2*wl/wg;
      
      // complex_error_function(x,y,re,im);
      std::complex<double> z(x,y);
      std::complex<double> cerf=Faddeeva::w(z);
      return a*real(cerf)*2.0*sqrtln2/sqrtpi/wg;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Convolution of a Fano profile with a gaussian
 *
 * @param Fano profile convoluted with a gaussian
 * @param E  Energy axis
 * @param E0 Peak center energy
 * @param a  Peak area
 * @param q  Asymmetry parameter
 * @param wl Lorentzian linewidths (FWMH)
 * @param wg Gaussian linewidths (FWHM)
 *
 * @return Height of the profile at the given energy (goes to zero for E->+-infinity)
*/
double fanovoigt(double E, double E0, double a, double q, double wl, double wg)
{
  double q2 = q*q;
  if (fabs(q)<1e-10) q=1e-10;
  if (q2>=1.0e6) 
    {//symmetric peak
      return voigt(E,E0,a,wl,wg);
    }
  else
    {//asymmetric peak
      if (fabs(wg)<1e-20)
	{
	  return fano(E,E0,a,q,wl);
	}
      else
	{
	  const double sqrtln2 = 0.832554611; 
	  const double sqrtpi  = 1.772453851;
	  
	  double x, y, re, im;
	  
	  x = 2.0*sqrtln2*(E0-E)/wg; // the sign matters!
	  y = sqrtln2*wl/wg;
	  //complex_error_function(x,y,re,im);
	  std::complex<double> z(x,y);
	  std::complex<double> cerf=Faddeeva::w(z);
	  re = real(cerf);
	  im = imag(cerf);
	  
	  double qfac1 = q2-1.0;
	  double qfac2 = 2.0*q;
	  double fac   = 2.0*fabs(a/qfac1)/sqrtpi/wl;
	  
	  return fac*y*(re*qfac1-im*qfac2);
	}
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Convolution of a Fano profile with a gaussian (returning also partial derivatives with respecte to parameters)
 *
 * 2018-07-10: Attention needs to be updated because of q² -> q²-1 substitution
 *
 * @param E  energy axis
 * @param E0 peak center energy
 * @param a  peak area (in the limit q -> infinity)
 * @param q  asymmetry parameter
 * @param wl lorentzian linewidths (FWMH)
 * @param wg gaussian linewidths (FWHM)
 * @param df_de0 on exit, partial derivative of the profile with respect to E0
 * @param df_da  on exit, partial derivative of the profile with respect to a
 * @param df_dq  on exit, partial derivative of the profile with respect to q
 * @param df_dwl on exit, partial derivative of the profile with respect to wl
 * @param df_dwg on exit, partial derivative of the profile with respect to wg
 * @return profile height at the given energy (goes to zero for E->+-infinity)
*/
double fanogvoigtderiv(double E, double E0, double a, double q, double wl, double wg,
		       double &df_dE0, double &df_da, double &df_dq, double &df_dwl, double &df_dwg)
{
  const double sqrtln2 = 0.832554611; 
  const double pi      = 3.141592654;
  const double sqrtpi  = 1.772453851;
  
  if (fabs(q)<1e-10) q=1e-10;
  
  double q2 = q*q;
  if (q2>=1.0e6) 
    {//symmetric peak
      if (fabs(wg)<1e-20)
	{   // Lorentz
	  double x = 2.0*(E-E0)/wl;
	  double fac = 2.0*a/wl/pi;
	  double f = fac/(1.0+x*x);

	  double df_dx = -2.0*x*f/(1.0+x*x);
	  
	  df_dE0 = -df_dx*2.0/wl;
	  df_dq  = 0.0;
	  df_dwg = 0.0;
	  df_dwl = -(f+x*df_dx)/wl;
	  df_da = f/a;
	  
	  return f;
	}
      else
	{   // Voigt
	  double fac, x, y, re, im;
	  
	  fac = a*2.0*sqrtln2/sqrtpi/wg;
	  x = 2.0*sqrtln2*(E-E0)/wg;
	  y = sqrtln2*wl/wg;
	  
	  //	complex_error_function(x,y,re,im);
	  std::complex<double> z(x,y);
	  std::complex<double> cerf=Faddeeva::w(z);
	  re = real(cerf);
	  im = imag(cerf);
	  
	  double dre_dx = -2*(x*re-y*im);
	  double dre_dy =  2*(y*re+x*im-1/sqrtpi);
	  
	  double f = fac*re;
	  
	  double df_dx = fac*dre_dx;
	  double df_dy = fac*dre_dy;
	  
	  df_dE0 = -df_dx*2.0*sqrtln2/wg;
	  df_dq = 0.0;
	  df_dwg = -(f+x*df_dx+y*df_dy)/wg;
	  df_dwl = y*df_dy/wl;
	  df_da = f/a;
	  
	  return f;
	}
    }
  else
    {//asymmetric peak
      if (fabs(wg)<1e-20)
	{ // Fano
	  double x = 2.0*(E-E0)/wl;
	  double xq = 1.0+x/q;
	  double x2 = x*x+1.0;
	  double fac = 2.0*a/wl/pi;
	  double f = fac*(xq*xq/x2-1);
	  
	  double df_dx = 2*xq/(x2*x2)*(1.0/q-x);
	  
	  df_dE0 = -2.0*df_dx/wl;
	  df_dq  = -2.0*fac*xq/(q2*x2);
	  df_dwg = 0.0;
	  df_dwl = -(f+x*df_dx)/wl;
	  df_da = f/a;
	  
	  return f;
	}
      else
	{ // FanoVoigt
	  double fac, qfac1, qfac2, x, y, re, im;
	  
	  fac   = a*2.0*sqrtln2/sqrtpi/wg;
	  qfac1 = 1.0-1.0/q2;
	  qfac2 = 2.0/q;
	  x = 2.0*sqrtln2*(E0-E)/wg;
	  y = sqrtln2*wl/wg;
	  
	  // complex_error_function(x,y,re,im);
	  std::complex<double> z(x,y);
	  std::complex<double> cerf=Faddeeva::w(z);
	  re = real(cerf);
	  im = imag(cerf);
	  
	  double dre_dx = -2*(x*re-y*im);
	  double dre_dy =  2*(y*re+x*im-1/sqrtpi);
	  double dim_dx = -dre_dy;
	  double dim_dy =  dre_dx;
          
	  double f = fac*(re*qfac1-im*qfac2);

	  double df_dre = fac*qfac1;
	  double df_dim =-fac*qfac2;
	  double df_dx  = df_dre*dre_dx+df_dim*dim_dx;
	  double df_dy  = df_dre*dre_dy+df_dim*dim_dy;
	  
	  df_dE0 = df_dx*2.0*sqrtln2/wg;
	  df_dq  = 2.0*fac*(re/q-im)/q2;
	  df_dwg = -(f+x*df_dx+y*df_dy)/wg;
	  df_dwl = y*df_dy/wl;
	  df_da  = f/a;
          
	  return f; 
	}
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Convolution of a lorentzian with a normalized trapezoidal 
 *
 * @param E  energy axis
 * @param E0 peak center energy
 * @param a peak area
 * @param wl lorentzian width (FWHM)
 * @param wb trapezoidal base width
 * @param wt trapezoidal top width
 *
 * @return profile height at the given energy 
 */
double trapezlorentzian(double E, double E0, double a, double wl, double wb, double wt)
{
  const double pi = 3.141592654;
	
  if (wb<wt) wb=wt;
	 
  double peak =  2.0*a/pi/(wb+wt)*(atan(2.0*(E+0.5*wt-E0)/wl)-atan(2.0*(E-0.5*wt-E0)/wl));
	 
  if (wt<wb)
    {
      double wfac = 1.0/(wb*wb-wt*wt);
      peak += -4.0*a/pi*wfac*(E-E0-0.5*wb)*(atan(2.0*(E-0.5*wt-E0)/wl)-atan(2.0*(E-0.5*wb-E0)/wl));
      peak +=  4.0*a/pi*wfac*(E-E0+0.5*wb)*(atan(2.0*(E+0.5*wb-E0)/wl)-atan(2.0*(E+0.5*wt-E0)/wl));
      peak +=  wfac*a*wl/pi*log( (pow(E-0.5*wt-E0,2)+pow(0.5*wl,2)) / (pow(E-0.5*wb-E0,2)+pow(0.5*wl,2)) );
      peak += -wfac*a*wl/pi*log( (pow(E+0.5*wb-E0,2)+pow(0.5*wl,2)) / (pow(E+0.5*wt-E0,2)+pow(0.5*wl,2)) );
      
    }
  return peak;
}
////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Rydberg series of FanoGauss peaks
 *
 * Quantum defect formula for photoionization cross section, 
 * [Tully et al., Astron. Astrophys. 211 (1989) 485].
 * The experimental photon energy spread is accounted for by convolution with a gaussian.
 *
 * @param x independent variable 
 * @param dx high-n cut-off parameter (should be of the order of <i>wg</i>)
 * @param wG gaussian width
 * @param zeff effective nuclear charge
 * @param slim energy of series limit
 * @param m0,m1   quantum defect parameters of series
 * @param q0,q1   Fano asymmetry parameters of series
 * @param w0,w1   natural line width parameters of series
 *
 * @return resonance cross section at energy <i>x</i>
*/
double RydbergSeries(double x, double dx, double wG, double zeff, double slim, 
		     double m0, double m1, double q0, double q1, double w0, double w1)
{
  const double Ryd = 13.6057;
  const double pi  =  3.1415926535;
  const double sqrtln2 = 0.832554611; 
  const double sqrtpi  = 1.772453851;
  
  if (x<=0.0) return 0.0;
  
  double z2Ryd = zeff*zeff*Ryd;
  double eps = (x-slim)/z2Ryd;
  
  // if the Rydberg resonances are too closely spaced the numerical
  // evaluation becomes error prone. Therefore high Rydberg resonances are
  // cut off if two neighboring resonance positions are smaller than dx
  double epscut = dx > 0 ? -exp(2.0*log(dx/(2*z2Ryd))/3.0) : 0;
  if (eps>epscut) return 1.0;
  
  double m = m0+eps*m1;  // quantum defect
  double q = q0+eps*q1;  // asymmetry parameter
  double w = w0+eps*w1;  // width parameter
  if (w<1e-20)
    {
      w = 1e-20;
    }
  double t = tanh(pi*w);
  double xx = tan(pi*(1.0/sqrt(-eps)+m));
  
  double fano;
  
  if (wG==0.0)
    {
      fano = (xx < 1e10*t) ? pow(xx+q*t,2)/(t*t+xx*xx) : pow(1.0+q*t/xx,2);
    }
  else
    { // convolution of Fano profiles with a Gaussian
      double Rez, Imz, Rew, Imw;
      
      double swG = sqrtln2/wG;
      Rez = -xx*swG;
      Imz =  t*swG;
      // complex_error_function(Rez,Imz,Rew,Imw);
      std::complex<double> z(Rez,Imz);
      std::complex<double> cerf=Faddeeva::w(z);
      Rew = real(cerf);
      Imw = imag(cerf);
      
      fano = Imz*((q*q-1.0)*Rew-2.0*q*Imw)*sqrtpi+1.0;
    }
  return fano;
}
