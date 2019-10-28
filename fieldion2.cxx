/**
 * @file fieldion2.cxx
 *
 * @brief More rigorous calculation of field-ionization rates (unfinished)
 *                                                                         
 * @verbatim
  $Id: fieldion2.cxx 446 2017-08-28 16:08:49Z iamp $
 @endverbatim
 *
 * Stefan Schippers@physik.uni-giessen.de
*/                                                                         

#include <stdio.h>
#include <math.h>
#include "hydromath.h"
#include "numerov.h"

using namespace std;

const int terms=200; 
const double accuracy = 1e-20;

double Mexpansion_lo(double x, double Z1, double E, double F, int m)
  {
  double eps = accuracy;
  double x2 = x*x;

  double sum = 1.0, c1 = sum*x2, c2 = 0.0, c3 = 0.0;
  int n;
  for(n=1; n<=terms; n++)
    {
    double c0 = -(Z1*c1+2.0*E*c2-F*c3)/(4.0*n*(n+m));
    sum += c0;
    if ((n>3) && (fabs(c0/sum)<eps)) break;
    c3 = c2*x2; c2 = c1*x2; c1 = c0*x2;
    }
  if (n == terms) 
    {
    printf("\n WARNING: desired accuracy not achieved with %6d",terms);
    printf("\n          expansion terms in function Mexpansion_lo");
    printf("\n          Increase max. no. of terms in fieldion.cxx.\n");
    }
  return sum*pow(x,m)*sqrt(x);
  }

//////////////////////////////////////////////////////////////////////////

double Mexpansion_hi(double x, double Z1, double E, double F, int m)
  {
  double eps = accuracy;

  double xinv = 1.0/x;
  double sqrtF = sqrt(F);
  double EsqrtF = E/sqrtF;
  double exparg = -(sqrtF*x*x/3.0-EsqrtF)*x;
  if (exparg<-200) return 0.0; 

  double sum = 1.0, c1 = sum*xinv, c2 = 0.0, c3 = 0.0;
  int n;
  for(n=1; n<=terms; n++)
    {
    double c0 = -( (Z1+EsqrtF*EsqrtF)*c1 - 2.0*EsqrtF*(n-1)*c2 +
                   + ((n-2)*(n-1)-m*m+0.25)*c3  )/(2.0*sqrtF*n);
    sum += c0;
    if ((n>3) && (fabs(c0)<eps)) break;
    c3 = c2*xinv; c2 = c1*xinv; c1 = c0*xinv;
    }
  if (n == terms) 
    {
    printf("\n WARNING: desired accuracy not achieved with %6d",terms);
    printf("\n          expansion terms in function Mexpansion_hi");
    printf("\n          Increase max. no. of terms in fieldion.cxx.\n");
    }
  return sum/x*exp(exparg);
  }


//////////////////////////////////////////////////////////////////////////

int dNexpansion_hi(double x, double Z2, double E, double F, int m, int* ntrm,
             double* dNda1, double* dNdb1, double* dNdxda1, double* dNdxdb1)
  {
// calculate the asymptotic expansion of the unbound wavefunction's N(x,a1,b1)
// partial derivatives. 

  double eps = accuracy;
  double x2 = x*x, xinv = 1.0/x;
  double sqrtF = sqrt(F), EsqrtF = E/sqrtF, m4 = m*m-0.25;
  double arg = (sqrtF*x2/3.0+EsqrtF)*x, darg = sqrtF*x2+EsqrtF;

  double aa1 = xinv, ba1 = 0.0, ab1 = 0.0, bb1 = xinv; 
  
  double aasum = aa1, aansum = aa1, aa2 = 0.0, aa3 = 0.0;
  double basum = 0.0, bansum = 0.0, ba2 = 0.0, ba3 = 0.0;   
  double absum = 0.0, abnsum = 0.0, ab2 = 0.0, ab3 = 0.0;
  double bbsum = bb1, bbnsum = bb1, bb2 = 0.0, bb3 = 0.0;   

  aa1 *= xinv; bb1 *= xinv; 

  int n;

  for(n=2; n<=terms; n++)
    {
    double n1 = (n-1.0), n2 = (n-2.0), n3 = (n-3.0)*n2;
    double c3 = n3-m4, c2 = 2.0*EsqrtF*n2, c1 = EsqrtF*EsqrtF-Z2;
    double c0 = 2.0*sqrtF*n1;

    double aa0 = ( c3*ba3 - c2*aa2 - c1*ba1)/c0;
    double ba0 = (-c3*aa3 - c2*ba2 + c1*aa1)/c0;
    double ab0 = ( c3*bb3 - c2*ab2 - c1*bb1)/c0;
    double bb0 = (-c3*ab3 - c2*bb2 + c1*ab1)/c0;
               
    aasum += aa0; aansum += n*aa0;
    basum += ba0; bansum += n*ba0;
    absum += ab0; abnsum += n*ab0;
    bbsum += bb0; bbnsum += n*bb0;

    if ((n>3) && (fabs(aa0)<eps) && (fabs(ba0)<eps)
              && (fabs(ab0)<eps) && (fabs(bb0)<eps)) break;

    aa3 = aa2*xinv; aa2 = aa1*xinv; aa1 = aa0*xinv;
    ba3 = ba2*xinv; ba2 = ba1*xinv; ba1 = ba0*xinv;
    ab3 = ab2*xinv; ab2 = ab1*xinv; ab1 = ab0*xinv;
    bb3 = bb2*xinv; bb2 = bb1*xinv; bb1 = bb0*xinv;

    }
  (*ntrm) = n;
  if (n == terms) 
    {
    printf("\n WARNING: desired accuracy not achieved with %6d",terms);
    printf("\n          expansion terms in function Nexpansion_hi");
    printf("\n          Increase max. no. of terms in fieldion.cxx.\n");
    return 1;
    }

  double s=sin(arg), c = cos(arg);

  (*dNda1) = aasum*s + basum*c;
  (*dNdb1) = absum*s + bbsum*c;
  (*dNdxda1) = (aasum*darg-bansum*xinv)*c - (basum*darg+aansum*xinv)*s;
  (*dNdxdb1) = (absum*darg-bbnsum*xinv)*c - (bbsum*darg+abnsum*xinv)*s;
  return 0;
  }


//////////////////////////////////////////////////////////////////////////

double Nexpansion_hi(double x, double Z2, double E, double F, int m,
                     double a1, double b1)
  {
  double eps = accuracy;
  double x2 = x*x, xinv = 1.0/x;
  double sqrtF = sqrt(F), EsqrtF = E/sqrtF, m4 = m*m-0.25;
  double arg = (sqrtF*x2/3.0+EsqrtF)*x, darg = sqrtF*x2+EsqrtF;

  a1 *= xinv; b1 *= xinv; 

  double asum = a1, a2 = 0.0, a3 = 0.0;
  double bsum = b1, b2 = 0.0, b3 = 0.0;   

  a1 *= xinv; b1 *= xinv; 

  int n;

  for(n=2; n<=terms; n++)
    {
    double n1 = (n-1.0), n2 = (n-2.0), n3 = (n-3.0)*n2;
    double c3 = n3-m4, c2 = 2.0*EsqrtF*n2, c1 = EsqrtF*EsqrtF-Z2;
    double c0 = 2.0*sqrtF*n1;

    double a0 = ( c3*b3 - c2*a2 - c1*b1)/c0;
    double b0 = (-c3*a3 - c2*b2 + c1*a1)/c0;
               
    asum += a0;
    bsum += b0;

    if ((n>3)&&(fabs(a0)<eps)&&(fabs(b0)<eps)) break;

    a3 = a2*xinv; a2 = a1*xinv; a1 = a0*xinv;
    b3 = b2*xinv; b2 = b1*xinv; b1 = b0*xinv;

    }
  if (n == terms) 
    {
    printf("\n WARNING: desired accuracy not achieved with %6d",terms);
    printf("\n          expansion terms in function Nexpansion_hi");
    printf("\n          Increase max. no. of terms in fieldion.cxx.\n");
    }
  return  asum*sin(arg) + bsum*cos(arg);
  }

//////////////////////////////////////////////////////////////////////////
  
int roots_xi(double E, double Z1, double F, int m, 
             double *xi1, double *xi2, double *xi3)
    {
     return cuberoot(F,-2.0*E,-4.0*Z1,m*m-1.0,xi1,xi2,xi3);
    }
 
//////////////////////////////////////////////////////////////////////////
  
int roots_eta(double E, double Z2, double F, int m,
             double *eta1, double *eta2, double *eta3)
    {
    return cuberoot(-F,-2.0*E,-4.0*Z2,m*m-1.0,eta1,eta2,eta3);
    }

//////////////////////////////////////////////////////////////////////////
  

int turningpoint(double E, double Z2, double F, int m, double *etatp)
   {
   const int maxiter = 1000;
   const double eps = 1e-7;
   double xnew = *etatp;
   double xold = xnew*10.0*eps;
   double mm1 = m*m+1.0;
   int iter = 0;
   for (;;)
     {
     if (fabs(1.0-xnew/xold)<eps) break; 
     iter++;
     if (iter>maxiter) return 0;
     xold = xnew;
     double sx = sqrt(xold);
     double numerator = -2.0*mm1+4.0*Z2*xold;
     double denominator = ((-2.0*E*xold-4.0*Z2)*xold+mm1)*xold;
     xnew = xold - numerator/denominator;
     printf("%12.6g %12.6g\n",xold,xnew);
     }
   *etatp = xnew;
   return 1;
   } 

//////////////////////////////////////////////////////////////////////////
void  Potential(double xlo, double xhi, double dx, int npts, 
                double E, double F, int m, double *V)
{
  for (int n=0; n<npts; n++)
    {
    double x = xlo+n*dx;
    if (x==0.0) x = dx*dx;
    double x2 = x*x;
    V[n] = (m*m-0.25)/x2 - 2.0*E*x2 + F*x2*x2;
    }
}

//////////////////////////////////////////////////////////////////////////
double BE, BF; int Bm; // global variables for passing E, F, and m

void Boundaries(double xlo, double xhi, double dx, int npts,
                  int *nforward, int *nbackward, double *Psi, double EV)
  {
  int n;

  Psi[0] = 0.0; 
  (*nforward)=1;
  double x = 0.0;
  for(;x<=0.02;)
    {
    x+=dx;
    Psi[(*nforward)] = Mexpansion_lo(x,EV,BE,BF,Bm);
    (*nforward)++;
    }

  printf("\n boundary conditions for outward integration:\n");
  for (n=0; n< (*nforward) ; n++) printf("%5d %10.5g\n",n,Psi[n]);

  if (BF<0) return; // in case of the unbound wavfunction
  
  Psi[npts-1] = Mexpansion_hi(xhi,EV,BE,BF,Bm);
  x = xhi-dx;
  double Psihi = Mexpansion_hi(x,EV,BE,BF,Bm);  
  Psi[npts-2] = Psihi;
  (*nbackward) = 2;
  for (;fabs(Psihi)<1e-80;)
    {
    x -= dx;
    double Psihi = Mexpansion_hi(xhi-dx,EV,BE,BF,Bm);  
    (*nbackward)++;
    Psi[npts-(*nbackward)] = Psihi;
    }    

  printf("\n boundary conditions for inward integration:\n");
  for (n=npts-1; n >= npts-(*nbackward); n--) printf("%5d %10.5g\n",n,Psi[n]);

}

//////////////////////////////////////////////////////////////////////////

double calcZ1(double xlo, double xhi, double dx, int nodes,
              double E, double F, int m, int writewave)
  {
  int npts = 1+int(0.1+fabs(xhi-xlo)/dx);

  double *V = new double[npts];
  double *Psi = new double[npts];

  Potential(xlo, xhi, dx, npts, E, F, m, V);

// passing parameters to subroutine boundaries
  BE = E; BF = F; Bm = m;

  double Z1 = Eigenvalue(xlo, xhi, dx, npts, nodes, V, Psi, Boundaries);

  if (writewave)
     {
     WritePsi(xlo, xhi, dx, npts, V, Psi, Z1, Boundaries, "Mbound.dat", 0.01);
     }

  delete[] V; delete[] Psi;
  return 0.25*Z1;
  }

//////////////////////////////////////////////////////////////////////////
inline double sign(double r) { return (r<0) ? -1.0 : 1.0; }

//////////////////////////////////////////////////////////////////////////
double match(double xlo, double xhi, double dx, double *xmatch, double Z2,
             double E, double F, int m, int writewave)
  {
// Matching the asymptotic expansion of the unbound wavefunction to
// the solution of the waveequation at a = xhi-dx, 
// returns the phaseshift of the scattering wavefunction and 
// on request writes out the unbound wavefunction
  double eps = sqrt(accuracy);

  int npts = 1+int(0.1+fabs(xhi-xlo)/dx);
  int nmatch = int(0.1+fabs(*xmatch-xlo)/dx);

  double *V = new double[npts];
  double *Psi = new double[npts];
  double *DPsi = new double[npts];

  Potential(xlo, xhi, dx, npts, E, -F, m, V);

// passing parameters to subroutine boundaries
  BE = E; BF = -F; Bm = m;

  double Nmax = CalcPsi(xlo, xhi, dx, npts, V, Psi, DPsi, Z2, Boundaries);

  double xm = xlo+nmatch*dx;
  double N = Psi[nmatch], dNdx = DPsi[nmatch];
  for (int n=nmatch+1; (n<npts)&&(fabs(N/Nmax)>0.1)&&
                       (fabs(N)<eps)&&(fabs(dNdx)<eps); n++)
    {
    N = Psi[n]; dNdx = DPsi[n];
    xm += dx;
    }
  *xmatch = xm;
  if (writewave)
    {
    printf(" The matching point is at x =    %15.10g\n",xm);
    printf(" N and dN/dx at matching point : %15.10g %15.10g\n",N,dNdx);
    }

  int err, ntrm;
  double dNda1, dNdb1, dNdxda1, dNdxdb1;
  err = dNexpansion_hi(xm,Z2,E,F,m,&ntrm,&dNda1,&dNdb1,&dNdxda1,&dNdxdb1);
  if (err) return 0.0;

  double a1 =  N*dNdxdb1-dNdx*dNdb1;
  double b1 = -N*dNdxda1+dNdx*dNda1;

  double phaseshift = (a1 == 0.0) ? 2.0*atan(1.0)*sign(b1) : atan(b1/a1);

  if (writewave)
    {
    double denominator = dNda1*dNdxdb1-dNdb1*dNdxda1;  
    if (denominator != 0.0) { a1 /= denominator; b1 /= denominator; } 
    double amplitude = sqrt(a1*a1+b1*b1);
    N    = a1*dNda1   + b1*dNdb1;
    dNdx = a1*dNdxda1 + b1*dNdxdb1;
    printf(" N and dN/dx at matching point : %15.10g %15.10g\n",N,dNdx);
    printf(" %8d terms have been used in the asymptotic expansion\n",ntrm);


    printf("matching results in\n"); 
    printf("                             a1 : %15.10g\n",a1);
    printf("                             b1 : %15.10g\n",b1);
    printf("                      amplitude : %15.10g\n",amplitude);
    printf("                     phaseshift : %15.10g\n",phaseshift);
 
    FILE* fout = fopen("Nunbound.dat","w");
    double x=xlo,y,xout=xlo,dxout=0.01;
    for(int n=0; n<npts-1; n++)
      {
      if (x>=xout)
        {
        fprintf(fout,"%12.8g \t %12.8g\n",x,Psi[n]);
        xout += dxout;
        }
      x+=dx;
      }
    x = xhi;
    for (;x<=2*xhi;)
      {
      y = Nexpansion_hi(x,Z2,E,F,m,a1,b1);
      fprintf(fout,"%12.8g \t %12.8g\n",x,y);
      x += dxout;
      }
    fclose(fout);
    } // end if(writewave)

  delete[] V; delete[] Psi; delete[] DPsi;

  return phaseshift;
  }

//////////////////////////////////////////////////////////////////////////
  
double Epert(int N1, int N2, int M, double Z, double F)
   {
   // E-field shifted level energy from 4th order perturbation theory;
   double m = double(M), m2 = m*m, m4 = m2*m2;
   double n = N1+N2+fabs(m)+1.0, n2=n*n, n4 = n2*n2;
   double nz   = n/Z, nz2  = nz*nz, nz3 = nz2*nz, nz4  = nz2*nz2; 
   double nz7  = nz3*nz4, nz10 = nz7*nz3;
   double q    = N1-N2, q2   = q*q, q4   = q2*q2;

   double a0 = -0.5/nz2;
   double a1 = 1.5*nz*q;
   double a2 = 0.0625*nz4*(-17*n2+3*q2+9*m2-19);
   double a3 = 0.09375*nz7*q*(23*n2-q2+11*m2+39);
   double a4 = nz10*(-5487*n4-147*q4+549*m4-1806*n2*q2+3402*n2*m2+1134*m2*q2
                     -35182*n2-5754*q2+8622*m2-16211)/1024.0;
   return a0 + F*(a1+F*(a2+F*(a3+F*a4)));
   }

//////////////////////////////////////////////////////////////////////////

double Z1pert(int N1, int N2, int M, double E, double F)
   {
     /// Z1 in 4th order perturbation theory
   
   double m = fabs(double(M)), m2 = m*m, m4 = m2*m2;
   double n = N1+N2+m+1, nz = sqrt(0.5/fabs(E));
   double f = F*nz*nz*nz;

   double c = 0.25*(1.0-m2);
   double a = N1+0.5*(1.0+m), a2 = a*a, a3 = a*a2, a4 = a*a3, a5 = a*a4;
   double b1 = 0.25*(6*a2+2*c);
   double b2 = 0.0625*(-68*a3+9*m2*a-19*a);
   double b3 = 0.015625*(1500*a4-258*m2*a2+918*a2+2.75*m4-35.5*m2+32.75);
   double b4 = 0.00390625*(-42756*a5+8910*m2*a3-46810*a3-227.25*m4*a+
                           3694.5*m2*a-5677.25*a);
   return (a + f*(b1+f*(b2+f*(b3+f*b4))))/nz;
   }

//////////////////////////////////////////////////////////////////////////

void FIhydrogen(void)
{
  const double Fau = 5.1422082e9; // V/cm
  const double Tau = 2.4188843e-17; // s

  double eps = 1e-10;
  int MAXTRY = 10;
  int k,n1,n2,m,n,npts,absm;
  double E,F,Z1,Z2,Z;

  printf("\n Give electric field strength in kV/cm : ");
  scanf("%lf",&F);
  F *= 1000/Fau;
  printf("\n The electric field strength is %11.5g\n",F);
  printf("\n Give the nuclear charge Z ............: ");
  scanf("%lf",&Z);
  printf("\n Give quantum numbers n1, n2 and m ....: ");
  scanf("%d %d %d",&n1, &n2, &m);

  E = Epert(n1,n2,m,Z,F);
  double neff = Z*sqrt(0.5/fabs(E));
  printf("\n The level energy from 4th order pert. theory is %11.5g a.u.",E);
  printf("\n corresponding to an effective main quantum number %8.5f\n",neff);
  //absm = m>0 ? m : -m;
  //n  = n1+n2+absm+1;
  Z1 = Z1pert(n1,n2,m,E,F);
  Z2 = Z-Z1;
  printf("\n Z1 and Z2 from 4th order perturbation theory are");
  printf(" %10.5f, %10.5f\n",Z1,Z2);

  double xi1, xi2, xi3, eta1,eta2,eta3;
  int nxi = roots_xi(E,Z1,F,m,&xi1,&xi2,&xi3);
  int neta = roots_eta(E,Z2,F,m,&eta1,&eta2,&eta3);
  printf("\n %3d roots for  Pxi",nxi);
  if (nxi>=1) printf(" : %12.6f",xi1);
  if (nxi>=2) printf(" : %12.6f",xi2);
  if (nxi>=3) printf(" : %12.6f",xi3);
  printf("\n %3d roots for Peta",neta);
  if (neta>=1) printf(" : %12.6f",eta1);
  if (neta>=2) printf(" : %12.6f",eta2);
  if (neta>=3) printf(" : %12.6f",eta3);
  eta3 =60.0;
  int tp = turningpoint(E,Z1,F,m,&eta3);
  if (tp>0) { printf(" : %12.8g",eta3); }
      else { printf(" too many iterations\n"); }
  printf("\n\n");

  /*
  double xlo=0.0, xhi, dx;
// upper integration limit
  double rbohr = 1.5*n*n/Z;
  xhi = sqrt(20.0*rbohr);
  printf(" xhi = %12.5g\n",xhi);  
  dx = 1e-3;
  double Z1old = 0.0;
  double dxhi = 0.1*(xhi-xlo);
  for(k=0;(k<MAXTRY)&&(fabs(Z1-Z1old) > eps);k++)
    {
    xhi += dxhi;
    Z1old = Z1; 
    Z1 = calcZ1(xlo,xhi,dx,n1,E,F,m,0);
    }
  Z1 = calcZ1(xlo,xhi,dx,n1,E,F,m,1);
  printf("\n Solution: Z1 = %15.10f, Z2 = %15.10f\n",0.25*Z1,0.25*Z2);
  dx = 1e-4;
  double xmatch = 0.85*xhi;
  double phase = match(xlo,xhi,dx,&xmatch,4*Z2,E,F,m,1);
  printf(" xmatch = %8.3f, Phase = %15.10g\n",xmatch,phase);
  */
}







