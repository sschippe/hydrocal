/** 
    @file osci.cxx
    @brief Calculation of hydrogenic bound-bound and bound-free oscillator strengths
 
    @par CREATION  
    @author Stefan Schippers
    @date 1997, 1999
    
    @par VERSION
    @verbatim
    $Id: osci.cxx 281 2014-07-23 13:24:45Z iamp $
    @endverbatim   


  Hydrogenic dipole oscillator strengths see
 
 - Kazem Omidvar and Patricia T. Guimares,                        
  The Astrophysical Journal Suplement Series 73 (1990) 555-602,       
 - H.A. Bethe and E.E. Salpeter, Quantum Mechanics of One and Two 
     Electron Systems n Handbuch der Physik Vol XXXV (Springer, 1957).

  Calculation of radial hydrogenic matrix elements using recursion formulae 
given by L. Infield and T.E.Hull, Rev. Mod. Phys 23 (1951) 21.
                                                                          

*/ 

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "clebsch.h"

using namespace std;

/* since the radial matrix elements are calculated iteratively             */
/* all bound-bound (bound-continuum) matrix elements for given n1 (k) and  */
/* n2 are calculated and stored in the vectors below of maximum length     */
/* nmax. The latest values of n1 (k) and n2 are stored.                    */
/* In order to avoid unnecessary calculations, in any summation of         */
/* oscillator strength sum over l before summing over n.                   */
      
const int nmax=10000;
double Rbbp[nmax], Rbbm[nmax];
double Rbcp[nmax], Rbcm[nmax];
static int n1bb_calculated = 0;
static int n2bb_calculated = 0;
static int n2bc_calculated = 0;
static double kbc_calculated = -1;



//------------------------------------------------------------------------------------------------------
/**
 * @brief recursive calculation of bound-bound radial matrix elements for all l
 *
 * The results for l->l-1 and l->l+1 are stored in the global arrays Rbbm[] 
 * and Rbbp[], respectively.
 * The calculation will be carried out only if n1 or n2 have changed.
 *
 * @param k wave number of continuum electron
 * @param n principal quantum number of bound electron
 * 
 */
void CalcRbb(int n1, int n2)
  {
   if ((n1bb_calculated==n1) && (n2bb_calculated==n2)) return;

// start value
  int n, nii;
  int nm  = n1-n2;
  int np  = n1+n2;
  int n22 = n2*2;
  double nn2 = 2*sqrt(double(n1*n2));
  double xm  = nm/nn2;
  double xp  = nn2/np;
  double fac = n1*n2*xp/double(np*np);

  if (nm==1) fac *= nn2;

  /* substituded by Geralds numerically more stable version below (11.8.2004) 
  for (n=2;    n<nm;    n++) fac *= xm/sqrt(n);
  for (n=2;    n<n22;   n++) fac *= xp;
  for (n=n22;  n<=np;   n++) fac *= xp*sqrt(n);
  */
  for (n=0; n<=nm; n++)
  {
      if ((n>=2) && (n<nm)) fac *= xm/sqrt(double(n));
      nii = n+n22;
      if ((nii>=n22) && (nii<=np)) fac *= xp*sqrt(double(nii));
  }
  for (n=2; n<n22; n++) fac *= xp;

  Rbbm[n2] = fac; 
  Rbbp[n2] = 0.0;

  double a,ap,am;
  for (int l=n2-1;l>=0;l--) 
    {
    int l1 = l+1;
    a = 2*sqrt(double((n2+l)*(n2-l)))/double(n2);
    am = (2*l+1)*sqrt(double((n1+l1)*(n1-l1)))/double(n1*l1);
    ap = sqrt(double((n2+l1)*(n2-l1)))/double(n2*l1);
    Rbbm[l] = (am*Rbbm[l1]+ap*Rbbp[l1])/a;
    a = 2*sqrt(double((n1+l)*(n1-l)))/double(n1);
    am = sqrt(double((n1+l1)*(n1-l1)))/double(n1*l1);
    ap = (2*l+1)*sqrt(double((n2+l1)*(n2-l1)))/double(n2*l1);
    Rbbp[l] = (am*Rbbm[l1]+ap*Rbbp[l1])/a;
    }
  n1bb_calculated = n1;
  n2bb_calculated = n2;
  }

//------------------------------------------------------------------------------------------------------
/**
 * @brief recursive calculation of bound-continuum radial matrix elements for all l
 *
 * The results for l->l-1 and l->l+1 are stored in the global arrays Rbcm[] 
 * and Rbcp[], respectively.
 * The caluclation will be carried out only if k or n have changed.
 *
 * @param k wave number of continuum electron
 * @param n principal quantum number of bound electron
 * 
 */
void CalcRbc(double k, int n2)
  { 

  if ((n2bc_calculated==n2) && (kbc_calculated==k)) return;

// start values
  int n;
  double ksqr = k*k;
  double n2sqr = n2*n2;
  double n2d = n2;
  double n22 = 2*n2d;
  double x = 4.0*n2d*k/(ksqr+n2sqr);
  double fac = sqrt(n2d*0.5)*0.25*x*x*ksqr; 
  for (n=1; n<=n2; n++) fac *= x*sqrt((ksqr+n*n)/double(n*(n22-n)));
  fac *= exp(-2.0*k*atan(n2d/k));
  if (k<10.0) fac /= sqrt(1.0-exp(-8*atan(1.0)*k));

  Rbcm[n2] = fac; 
  Rbcp[n2] = 0.0;
  double a,ap,am;
  for (int l=n2-1;l>=0;l--) 
    {
    int l1 = l+1;
    a = 2*sqrt(double((n2+l)*(n2-l)))/double(n2);
    am = (2*l+1)*sqrt(k*k+l1*l1)/double(k*l1);
    ap = sqrt(double((n2+l1)*(n2-l1)))/double(n2*l1);
    Rbcm[l] = (am*Rbcm[l1]+ap*Rbcp[l1])/a;
    a = 2*sqrt(k*k+l*l)/double(k);
    am = sqrt(k*k+l1*l1)/double(k*l1);
    ap = (2*l+1)*sqrt(double((n2+l1)*(n2-l1)))/double(n2*l1);
    Rbcp[l] = (am*Rbcm[l1]+ap*Rbcp[l1])/a;
    }
  kbc_calculated = k;
  n2bc_calculated = n2;
  }


//------------------------------------------------------------------------------------------------------
/**
 * @brief hydrogenic bound-bound oscillator strength
 *
 * @param n1 principal quantum number of initial level
 * @param l1 orbital angularmomentum quantum number of initial level 
 * @param n2 principal quantum number of final level
 * @param l2 orbital angularmomentum quantum number of final level 
 *
 * @return bound-bound oscillator strength
 */
double fosciBB(int n1, int l1, int n2, int l2)
  {
  if (n1>n2) {CalcRbb(n1,n2);} else {CalcRbb(n2,n1);};
  double n1sqr = n1*n1;
  double n2sqr = n2*n2;
  double ll    = 3.0*(2.0*l1+1.0);
  if (l2==(l1-1))
     {
     double r = n1>n2 ? Rbbm[l1] : Rbbp[l2+1];
     return (1.0/n1sqr-1.0/n2sqr)*l1/ll*r*r;
     }
  else if (l2==(l1+1))
     {
     double r = n1>n2 ? Rbbp[l1+1] : Rbbm[l2];
     return (1.0/n1sqr-1.0/n2sqr)*l2/ll*r*r;
     }
  else return 0;
  }

//------------------------------------------------------------------------------------------------------
/**
 * @brief hydrogenic bound-free oscillator strength
 *
 * @param e scaled photon energy
 * @param n principal quantum number of bound electron
 * @param l orbital angularmomentum quantum number of bound electron
 *
 * @return bound-free oscillator strength
 */
double fosciBC(double e, int n, int l)
  {
  double k = 1.0/sqrt(e);

  CalcRbc(k,n);

  double fac = (e+1.0/n/n)/3/(2*l+1);
  double rm  = Rbcm[l+1];
  double rp  = l>0 ? Rbcp[l] : 0;

  return fac*((l+1)*rm*rm + l*rp*rp);
  }

//------------------------------------------------------------------------------------------------------
/**
 * @brief confluent hypergeometric function
 *
 * see http://dlmf.nist.gov/13.8 equation 13.8.1 for large c
 *
 * @param a parameter < 0
 * @param c parameter |c| >> 1
 * @param x variable
 *
 * @return function value
 */
double hypconfl(double a, double c, double x)
  {
  double k, sum, fac;

  if (a>0) return 0.0;

  k   = 1;
  fac = 1;
  sum = 1;

  while (a)
    { 
    fac = fac*a/c/k*x;
    sum += fac;
    a++;
    c++;
    k++;
    }
  return sum; 
  }

//------------------------------------------------------------------------------------------------------
/**
 * @brief bound-free matrix element squared for l->l+1 and zero energy
 *
 * @param n principal quantum number of bound electron
 * @param l orbital angularmomentum quantum number of bound electron
 *
 * @return bound-free matrix element squared
 */
double rpsqr0(int n, int l)
  {
  double n4 = n*4.0;
  double prod = pow(n4,3);
  double ll = 0;
  for(int m=-l; m<=l; m++)
    {
    ll++;
    prod *= (n+m)*n4/ll/ll;
    } 
  double hyp1 = hypconfl(-n+l+1,2*l+2,n4);
  double hyp2 = hypconfl(-n+l+2,2*l+3,n4);
  return prod*pow(hyp1+(n-l-1.0)/(l+1.0)*hyp2,2);
  }

//------------------------------------------------------------------------------------------------------
/**
 * @brief bound-free matrix element squared for l->l-1 and zero energy
 *
 * @param n principal quantum number of bound electron
 * @param l orbital angularmomentum quantum number of bound electron
 *
 * @return bound-free matrix element squared
 */
double rmsqr0(int n, int l)
  {
  double n4 = n*4.0;
  double prod = (n-l)*(n+l)*n4;
  double ll = 0;
  for(int m=-l+1; m<=l-1; m++)
    {
    ll++;
    prod *= (n+m)*n4/ll/ll;
    } 
  double hyp1 = hypconfl(-n+l+1,2*l,n4);
  double hyp2 = hypconfl(-n+l-1,2*l,n4);
  return prod*pow(hyp1-hyp2,2);
  }

//------------------------------------------------------------------------------------------------------
/**
 * @brief hydrogenic bound-free oscillator strength at threshold
 *
 * @param n principal quantum number of bound electron
 * @param l orbital angularmomentum quantum number of bound electron
 *
 * @return bound-free oscillator strength at threshold
 */
double fosciBC(int n, int l)
  {
  double r = l==0 ? rpsqr0(n,l) : l*rmsqr0(n,l)+(l+1)*rpsqr0(n,l);
  return r*exp(-4.0*n)/(2*l+1)/6.0;
  }

//------------------------------------------------------------------------------------------------------
/**
 * @brief hydrogenic bound-free oscillator strength at threshold
 *
 * @param n principal quantum number of bound electron
 * @param l orbital angularmomentum quantum number of bound electron
 * @param m magentic quantum number
 *
 * @return bound-free oscillator strength at threshold
 */
double fosciBC(int n, int l, int m)
  {
  double r;
  if (l==0)
    {
    r = rpsqr0(n,0);
    }
  else
    {
    double splus = 0.0, sminus = 0.0;
    for (int dm = -1; dm<=1; dm++)
      { 
      splus  += Strength(l,l+1,m,dm);
      sminus += Strength(l,l-1,m,dm);
      }
      r = l*rmsqr0(n,l)*splus+(l+1)*rpsqr0(n,l)*sminus;
    }
  return r*exp(-4.0*n)/(2*l+1)/6.0;
  }

//------------------------------------------------------------------------------------------------------
/**
 * @brief test subroutine for bound-bound oscillator strength
 */
void testOsciBB(void)
  {
  int   n1, l1, n2, l2;

  printf(" Give n,l and n',l' for n,l -> n',l' transition : ");
  scanf("%d %d %d %d",&n1, &l1, &n2, &l2);
  printf("\n\n fBB = %12.4G \n\n",fosciBB(n1,l1,n2,l2));
  }

//------------------------------------------------------------------------------------------------------
/**
 * @brief test subroutine for bound-free oscillator strength
 */
void testOsciBC(void)
  {   
  double x,f;
  int n,l;

  printf(" Give x (n*n*Rydberg), n, l : ");
  scanf("%lf %d %d",&x, &n, &l);
  f = x == 0 ? fosciBC(n,l) : fosciBC(x/n/n,n,l);
  printf("\n\n nfBC = %g /Ryd \n\n", f/n);
  }

















