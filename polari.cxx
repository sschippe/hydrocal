/**
 * @file polari.cxx
 *
 * @brief Depolarization factors for radiative transitions
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: polari.cxx 375 2016-02-05 16:40:12Z iamp $
   @endverbatim
 *
 */
#include <stdio.h>
#include <math.h>
#include "clebsch.h"
#include "lifetime.h"
#include "polari.h"

using namespace std;

inline int iround(double r)
  {
  return (r<0.0) ? (int) ceil(r-0.1) : (int) floor(r+0.1);
  }

double hfshift(int n, int l, double f, double j, double i, double gi, 
               double q, double z)
  {
  const double a0=0.529e-8; //  bohr radius in cm
  const double qunit=1e-24; //  nuclear quadrupole moment in 10^(-24) cm^2
  const double tunit=2.4188e-11; // atomic time unit in micro seconds
  const double mass=1836;        // proton mass / electron mass
  const double alpha=137.036;    // 1 over fine structure constant

  double a, b, c, r3, zn;

  zn = z/n;

  if (l==0)
    a = gi*zn*zn*zn*8/3.0/(tunit*alpha*alpha*2*mass);
  else
    a = gi*zn*zn*zn/(l+0.5)/j/(j+1)/(tunit*alpha*alpha*2*mass);

  if (l==0) r3 = 0; else r3 = zn*zn*zn/((l+0.5)*l*(l+1));

  if (i<1) 
    b = 0;
  else 
    b = r3/tunit*3*q*qunit/a0/a0/(16*i*(2*i-1)*j*(j+1));

  c = f*(f+1)-j*(j+1)-i*(i+1);

  return 0.5*a*c +b*(c*(c+1)-4/3*j*(j+1)*i*(i+1));
  }

double ehf(int n, int l, double f, double j, double i, double gi, 
           double q, double z)
//  energy of hyperfine-level in MHz
  {
  const double ry2mhz=3.289841951e10; // 1 Rydberg in Mhz
  const double alpha=137.036;         // 1 over fine structure constant

  double enlj, zn, zn2;

   zn   = z/n;
   zn2  = zn*zn;

   // hydrogenic fine-structure levles
   enlj = -zn2 - 0.25*zn2*zn2/alpha/alpha/(4*n/(j+0.5)-3)*ry2mhz;

   return enlj+hfshift(n,l,f,j,i,gi,q,z);
  }
/////////////////////////////////////////////////////////////////////////
/**
 * @brief depolarization factor for nuclear spin I>0
 *
 * for photon emission in ion-atom collision experiments
 *
 * involves calculation of hyperfine splittings of energy levels
 *
 * @param lt target width in mm
 * @param lp profile width in mm
 * @param v ion velocity  in millimeters per microseconds
 */
double g2i(int l, double s, double i, double gi, double q, int n, 
           double z, double v, double lt, double lp)
  {
  const double Pi = 3.141592654;
  const double maxiter = 100;  // max order of series expansion
  const double eps     = 1e-4; // accuracy of series expansion

  int k,kjmin,kjmax,kj1,kj2,kf1,kf2,f1min,f2min,f1max,f2max,converged;
  double halfj,halff,sum,j1,j2,f1,f2,sj1,sj2;
  double delta,gamma,factor,timefactor;
  double sin1, sin2, cos1, cos2, exp1, exp2, expk1, expk2, snew, sold;
  double gd, gk, gdk, t1, t2, fc, fs;

  gamma = 1e-6/2.0/Pi/hydrolife(z,n,l,0); // (natural line width in MHz)
  if (i<0.0) i=0.0; // i<0 is used for testing

  //  look whether j or f will be half integer
   halfj = 0.5*(iround(2*s) % 2);
   halff = 0.5*(iround(2*(s+i)) % 2);

   sum = 0;
   kjmin = iround(l-s-halfj);
   kjmax = iround(l+s-halfj);
   for (kj1 = kjmin;kj1 <= kjmax; kj1++)
     {
     j1 = kj1+halfj;
     f1min = iround(fabs(j1-i)-halff);
     f1max = iround(j1+i-halff);
     for (kj2 = kjmin; kj2 <= kjmax; kj2++)
       {
       j2 = kj2+halfj;
       f2min = iround(fabs(j1-i)-halff);
       f2max = iround(j1+i-halff);
       for (kf1 = f1min; kf1 <= f1max; kf1++)
           {
           f1 = kf1+halff;
           for (kf2 = f2min; kf2 <= f2max; kf2++)
             {
             f2 = kf2+halff;

             sj1 = SixJ(j1,f1,i, f2,j2,2);
             sj1 = sj1*sj1;
             sj2 = SixJ(l,j1,s ,j2,l,2);
             sj2 = sj2*sj2;
             factor = (2*j1+1)*(2*j2+1)*(2*f1+1)*(2*f2+1);
             delta  = ehf(n,l,f2,j2,i,gi,q,z)-ehf(n,l,f1,j1,i,gi,q,z);
             gd   = 1.0/(gamma*gamma+delta*delta);
             if ((lp<=0) | (v<=0)) timefactor = gamma*gamma*gd;
             else if (lt<=0)
               { // target width is zero 
               t2   = lp/v;
               exp2 = exp(-gamma*t2);
               sin2 = sin(delta*t2);
               cos2 = cos(delta*t2);
               timefactor = gamma*(1-cos2*exp2)+delta*sin2*exp2;
               timefactor = timefactor*gamma*gd/(1-exp2);
               } // end else if ...
             else
               { // finite target width
               t2   = lp/v;
               exp2 = exp(-gamma*t2);
               sin2 = sin(delta*t2);
               cos2 = cos(delta*t2);
               t1   = (lp-lt)/v;
               exp1 = exp(-gamma*t1);
               sin1 = sin(delta*t1);
               cos1 = cos(delta*t1);
               timefactor =  (t2-t1)*gamma+log( (1-exp2)/(1-exp1) );

   // series expansion of Integral dt over Cos(delta*t)/(1-Exp(-gamma*t))
               snew  = 0;
               gk    = 1;
               expk1 = 1;
               expk2 = 1;
               k     = 0;
               converged = 1;
               do
                 {
                 sold  = snew;
                 gk    = gk + gamma;
                 expk1 = expk1*exp1;
                 expk2 = expk2*exp2;
                 gdk   = 1/(delta*delta+gk*gk);
                 fc    = gdk*(gk*gamma-delta*delta);
                 fs    = gdk*(gk+gamma)*delta;
                 fc    = fc*(cos2*expk2-cos1*expk1);
                 fs    = fs*(sin2*expk2-sin1*expk1);
                 snew  = sold+fc-fs;

                 if ((snew==0) & (fabs(sold) < eps))
                     converged = 0;
                 else if ((snew!=0) & (fabs(1.0-sold/snew) < eps))
                     converged = 0;
                 else
                     converged = 1; 

                 if (k>maxiter)
                   {
                   printf(" g2i not converged within %3d iterations\n",maxiter);
                   converged = 0;
                   }
                 k++;
                 }
               while (converged);

               timefactor = (timefactor + snew)*gd*gamma*v/lt;
             } // end  else 

           sum += sj1*sj2*factor*timefactor;
           } // end for  f2
         } // end for  f1
     } // end for  j2
   }   // end for  j1

   return sum/(2*s+1)/(2*i+1);
  }

/////////////////////////////////////////////////////////////////////////
/**
 * @brief depolarization factor for zero nuclear spin
 */
double g20(int l, double s)
  {
  double j, half, g2, jsixj;
  int  kj, kjmin, kjmax;

  g2 = 0.0;
  half = 0.5*(iround(2*s) % 2); // check whether spin is half integer
  kjmin = iround(l-s+half); 
  kjmax = iround(l+s+half);
  for (kj = kjmin ; kj <= kjmax; kj++)
    {
    j = kj-half;
    jsixj = (2*j+1)*SixJ(l,j,s,j,l,2.0);
    g2 += jsixj*jsixj;
    }
  g2 /= 2.0*s+1.0;
  return g2;
  }


/////////////////////////////////////////////////////////////////////////
/**
 * @brief depolarization factor for photon emission in ion-atom collision experiments
 */
double G2depol(int l, double s, double i, double gi, 
               double q, int n, double z,
               double v, double lt, double lp)
//double G2depol(int l, double s, double i = 0 , double gi = 4.43, 
//               double q = 0.14, int n = 3, double z = 1.0,
//               double v = 0.0, double lt = 0.0, double lp =0.0)
  {
  return (i==0) ?  g20(l,s) : g2i(l,s,i,gi,q,n,z,v,lt,lp);
  }
/*
int polcoeff(int l, int dl, int m, double s, double i, double gi,
              double q, int n, double z, double v, double targ, double prof,
              double &a, double &b)
  {
  double p, p_0, p_m, depol, sign, t20, t20_0, t20_m;
  double sigma, sigma_0, sigma_m;
  
  sign  =1.0-2*( (2*l+dl) % 2);
  depol = G2depol(l,s,i,gi,q,n,z,v,targ,prof);
  depol *= sign*SixJ(1.0,1.0,2.0,l,l,l+dl);

  sign = 1-2*(l % 2);
  sigma_0 = 1;
  t20_0 = sigma_0*sign*ThreeJ(l,l,2.0,0.0,0.0,0.0);
  p_0 = depol*sqrt(7.5)*t20_0/
        (sigma_0/(1.5*(2*l+1))+depol/sqrt(1.2)*t20_0);
  a = p_0;
  b = 1.0;
  if (m==0) return 1;

  sign = 1.0-2*((l-m) % 2);
  sigma_m = 2.0;
  t20_m = sigma_m*sign*ThreeJ(l,l,2.0,m,-m,0.0);
  p_m = depol*sqrt(7.5)*t20_m/
         (sigma_m/(1.5*(2*l+1))+depol/sqrt(1.2)*t20_m);
  sigma = sigma_0 + sigma_m;
  t20 = t20_0 + t20_m;
  p   = depol*sqrt(7.5)*t20/
         (sigma/(1.5*(2*l+1))+depol/sqrt(1.2)*t20);

  if (p_m==0) 
    {
    b = a / p - 1;
    a = 0;
    }
  else
  {
  return (i==0) ?  g20(l,s) : g2i(l,s,i,gi,q,n,z,v,lt,lp);
  }
*/
void polcoeff(int l, int dl, int m, double s, double i, double gi,
              double q, int n, double z, double v, double targ, double prof,
              double &a, double &b)
  {
  double p, p_0, p_m, depol, sign, t20, t20_0, t20_m;
  double sigma, sigma_0, sigma_m;
  
  sign  =1.0-2*( (2*l+dl) % 2);
  depol = G2depol(l,s,i,gi,q,n,z,v,targ,prof);
  depol *= sign*SixJ(1.0,1.0,2.0,l,l,l+dl);

  sign = 1-2*(l % 2);
  sigma_0 = 1;
  t20_0 = sigma_0*sign*ThreeJ(l,l,2.0,0.0,0.0,0.0);
  p_0 = depol*sqrt(7.5)*t20_0/
        (sigma_0/(1.5*(2*l+1))+depol/sqrt(1.2)*t20_0);
  a = p_0;
  b = 1.0;
  if (m==0) return;

  sign = 1.0-2*((l-m) % 2);
  sigma_m = 2.0;
  t20_m = sigma_m*sign*ThreeJ(l,l,2.0,m,-m,0.0);
  p_m = depol*sqrt(7.5)*t20_m/
         (sigma_m/(1.5*(2*l+1))+depol/sqrt(1.2)*t20_m);
  sigma = sigma_0 + sigma_m;
  t20 = t20_0 + t20_m;
  p   = depol*sqrt(7.5)*t20/
         (sigma/(1.5*(2*l+1))+depol/sqrt(1.2)*t20);

  if (p_m==0) 
    {
    b = a / p - 1;
    a = 0;
    }
  else
    {
    a = (a-p) / (p/p_m-1);
    b = a /  p_m;
    }
  }

void testHF(void)
  {
  const double s = 0.5; // electron spin (one-electron system)

  int n, l, kj, kjmin, kf, kfmin, kfmax;
  float  j, i, gi, q, z; 
  double f, half, delta;

  printf(" Give z, n, i,   gi [mu_n],    q [1E-24 cm-2] : ");
  scanf("%f %d %f %f %f",&z,&n,&i,&gi,&q);
  for (l=0;l<n;l++)
    {
    if (l==0) kjmin = 1; else kjmin = 0;
    for (kj = kjmin; kj <= 1; kj++)
      {
      j = l+kj-s;
      half = 0.5*(iround(2*(j+i)) % 2);
      kfmin = 1+iround(fabs(j-i)-half);
      kfmax = iround(j+i-half);
      for (kf = kfmin; kf <= kfmax; kf++) 
        {
        f = kf+half;
        delta = ehf(n,l,f,j,i,gi,q,z)-ehf(n,l,f-1,j,i,gi,q,z);
        printf(" l=%2d j=%4.1f f=%4.1f-> f'=%4.1f",l,j,f,f-1);
        printf(" delta E : %10.4f MHz\n",delta);
        } // end for kf
      } // end for kj
     } // end for l
  }

void testPol(void)
  {
  int n,l,dl,m;
  float s,i,gi,q,z;
  double a,b,p, pmin, pmax, asum, bsum;

  printf(" Give electron spin s .....................: ");
  scanf("%f",&s);
  printf(" Give nuclear spin i ......................: ");
  scanf("%f",&i);
  if (i!=0)
    {
    printf(" Give nuclear charge z ....................: ");
    scanf("%f",&z);
    printf(" Give principal quantum number n ..........: ");
    scanf("%d",&n);
    printf(" Give nuclear magnetic moment [mu_n] ......: ");
    scanf("%f",&gi);
    if (i>0.5)
      {
      printf(" Give nuclear quadrupole moment [1E-24cm2] : ");
      scanf("%f",&q);
      }
    }
  pmin =  1;
  pmax = -1;
  asum = 0;
  bsum = 0;
  printf("\n Give l -> l' : ");
  scanf("%d %d",&l,&dl);
  dl = dl-l;
  printf("\n           enumerator  denominator\n");
  for (m=0;m<=l;m++)
    {
    polcoeff(l, dl, m, s, i, gi, q, n, z, 0.0, 0.0, 0.0, a, b);
    asum += a;
    bsum += b;
    p = a/b;
    if (p>pmax) pmax = p;
    if (p<pmin) pmin = p;
    printf("    m=%2d%13.6f%13.6f\n",m,a,b);
    }
   printf("  sum=%15.6f%13.6f\n\n",asum,bsum);
   printf(" pmin = %10.4f\n",pmin);
   printf(" pmax = %10.4f\n\n",pmax);
  }

void testDepol(void)
  {
  float s,i,z,gi,q,v,t,p,depol;
  int n,l;
  char choice='0', answer='n';

  printf(" Give electron spin s .....................: ");
  scanf("%f",&s);
  printf(" Give nuclear spin i ......................: ");
  scanf("%f",&i);
  if (i!=0)
    {
    printf(" Give nuclear charge z ....................: ");
    scanf("%f",&z);
    printf(" Give principal quantum number n ..........: ");
    scanf("%d",&n);
    printf(" Give nuclear magnetic moment [mu_n] ......: ");
    scanf("%f",&gi);
    if (i>0.5)
      {
      printf(" Give nuclear quadrupole moment [1E-24cm2] : ");
      scanf("%f",&q);
      }
    }
  printf(" Give l : ");
  scanf("%d",&l);
  v = 0;
  t = 0;
  p = 0;
  depol = G2depol(l,s,i,gi,q,n,z,v,t,p);
  printf(" g2(%2d) = %10.6f\n\n",l,depol);
  }
