/**
 * @file fele.cxx
 *
 * @brief energy distribution functions to be used in convolutions
 *
 * @author Stefan Schippers
 * @verbatim
   $Id $
 @endverbatim
*/

#include <math.h>
#include "hydromath.h"

using namespace std;

/////////////////////////////////////////////////////////////////////
/**
 * @brief flattened Maxwellian energy distribution (internal use only)
 *
 * @param e energy
 * @param eint integration variable
 * @param ktpar parallel temperature
 * @param ktperp transversal temperature
 */
double flatmax(double e, double eint, double ktpar, double ktperp)
   {
   double kt1 = 1-ktpar/ktperp;
   double kt2 = sqrt(ktpar/kt1);
   double x = sqrt(e)/kt1/kt2;
   double h = sqrt(eint)/kt2;
   double erfdiff = daerf(x,h);
   double f=0;
   if (erfdiff)
     {
     double exparg = (e/kt1-eint)/ktperp;
     if (fabs(exparg)<277.0)
        {
        f = exp(exparg)*erfdiff*0.5/sqrt(kt1)/ktperp;
        }
     }
   return f;
   }

/////////////////////////////////////////////////////////////////////
/**
 * @brief Gaussian energy distribution
 *
 * @param x variable
 * @param x0 centroid
 * @param fwhm full width at half-maximum
 * @param dummy not used
 */
double felegauss(double x, double x0, double fwhm, double dummy)
   {
     const double c0 = 2.77258872; // 4*ln2;
     const double c1 = 0.93943728; // 1/sqrt(Pi)*sqrt(4*ln2) 
     double arg = (x-x0)/fwhm;
     return c1*exp(-c0*arg*arg)/fwhm;
   }


/////////////////////////////////////////////////////////////////////
/**
 * @brief flattened Maxwellian energy distribution
 *
 * @param e energy
 * @param eint integration variable
 * @param ktpar parallel temperature
 * @param ktperp transversal temperature
 */
double fecool(double e, double eint, double ktpar, double ktperp)
   {
   const double c1 = 0.173286795; // ln(2)/4
   const double c2 = 3.330218445; // 4*sqrt(ln(2))
   double emean = c1*ktperp*ktperp/ktpar; 
   if (e>1000*emean)
     {
     double fwhm = c2*sqrt(e*ktpar);
     return felegauss(e,eint,fwhm,fwhm);
     }
   else
     {
     return flatmax(e,eint,ktpar,ktperp);
     } 
   }


/////////////////////////////////////////////////////////////////////
/**
 * @brief trapezoidal energy distribution
 *
 * @param e energy
 * @param eint integration variable
 * @param wb bottom width
 * @param wt top width
 */
double trapezoid(double e, double eint, double wb, double wt)
{
    double x = e-eint;
    double y = 0;
    if (wt<wb)
    {
	double ww = (wt+wb)*(wb-wt);
	if ( (x>-0.5*wb) && (x<-0.5*wt) )
	{
	    y = 4.0/ww*(x+0.5*wb);
	}
	else if ((x>=-0.5*wt) && (x<=0.5*wt))
	{
	    y = 2.0/(wt+wb);
	}
	else if ( (x>0.5*wt) && (x<0.5*wb) )
	{
	    y = -4.0/ww*(x-0.5*wb);
	}
	else
	{
	y = 0;
	}
    }
    else // if wt>=wb ignore wb
    {
	if ((x>=-0.5*wt) && (x<=0.5*wt))
	{
	    y = 1.0/wt;
	}
	else
	{
	    y=0;
	}
    }
    return y;
}


/////////////////////////////////////////////////////////////////////
/**
 * @brief Maxwellian energy distribution
 *
 * @param kT temperature
 * @param E integration variable
 * @param dummy1 not used
 * @param dummy2 not used
 */
double fmaxwell(double kT, double E, double dummy1, double dummy2)
{
  const double c0 = 1.128379167; // 2/sqrt(pi)
  double x = E/kT;

  // compute sqrt(x)*exp(-x)
  double exparg = 0.5*log(x)-x;
  double expfac = exp(exparg);

  return c0*expfac/kT;
}



























