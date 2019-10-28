/**
 * @file numerov.cxx
 *
 * @brief Numerov's ODE solver
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: numerov.cxx 446 2017-08-28 16:08:49Z iamp $
   @endverbatim
 *
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "numerov.h"

using namespace std;

  
double DlogPsi(double dx, int npts, double *V, double *Psi, 
               int nforward, int nbackward, double EV, int *nodes)
{
  // Forward and backward integration are performed
  // up to a intermediate match point. 
  // The logarithimic derivative at the match point
  // and the number of nodes detected are returned.

  double dlog, h12 = dx*dx/12.0;
  
  int n=0, nmatch = int(0.53*npts);
  (*nodes) = 0;

  for (n=nforward;n<=nmatch+1;n++)
    {
    Psi[n] = (   Psi[n-1]*(2.0+10.0*h12*(V[n-1]-EV)) 
               - Psi[n-2]*(1.0-     h12*(V[n-2]-EV)) )
             /(1.0 - h12*(V[n]-EV));
    if ((n<=nmatch) && (Psi[n]*Psi[n-1]<0.0)) (*nodes)++;
    }
  dlog = (Psi[n]-Psi[n-2])/(2.0*dx*Psi[n-1]);

  for (n=npts-nbackward-1;n>=nmatch-1;n--)
    {
    Psi[n] = (   Psi[n+1]*(2.0+10.0*h12*(V[n+1]-EV)) 
               - Psi[n+2]*(1.0-     h12*(V[n+2]-EV)) )
             /(1.0 - h12*(V[n]-EV));
    if ((n<=nmatch) && (Psi[n]*Psi[n+1]<0.0)) (*nodes)++;
    }
  dlog += (Psi[n]-Psi[n+2])/(2.0*dx*Psi[n+1]);

  return dlog;
}

/////////////////////////////////////////////////////////////////////////////

double Eigenvalue(double xlo, double xhi, double dx, int npts, 
                  int nodes, double *V, double *Psi, BOUNDARIES b,
                  double EVlo, double EVhi)
   {
   const int MAXTRY = 100;
   const double eps = 1.e-10;

   int k, n, nfwd=2, nbwd=2, is_lower_limit, is_upper_limit;
   double EV=0.0, Dhi, Dlo, Dnew;
// First find start values using simple boundary conditions
   Psi[0] = 0.0; Psi[1] = dx;
   Psi[npts-1] = 0.0; Psi[npts-2] = dx;
   for (k=0; k<MAXTRY; k++)
     {
     EVlo*=2.0;
     Dlo = DlogPsi(dx,npts,V,Psi,nfwd,nbwd,EVlo,&n);
     is_lower_limit = ( (n<nodes) || ((n==nodes) && (Dlo >= 0.0)));
     if (is_lower_limit) break;
     }
   for (k=0; k<MAXTRY; k++)
     {
     EVhi*=2.0;
     Dhi = DlogPsi(dx,npts,V,Psi,nfwd,nbwd,EVhi,&n);
     is_upper_limit = ((n>nodes) || ((n==nodes) && (Dhi < 0.0)));
     if (is_upper_limit) break;
     }

// binary chop to get closer
   for (k=0; ((k<MAXTRY) && (EVlo != EVhi)); k++)
     {
     EV = 0.5*(EVlo+EVhi);
     // now apply more accurate boundary condtions
     (*b)(xlo,xhi,dx,npts,&nfwd,&nbwd,Psi,EV);
     Dnew = DlogPsi(dx,npts,V,Psi,nfwd,nbwd,EV,&n);
     is_upper_limit = ((n>nodes) || ((n==nodes) && (Dnew < 0.0)));
     if (is_upper_limit) { EVhi = EV; Dhi = Dnew; } 
                    else { EVlo = EV; Dlo = Dnew; }
     }
   
// interpolate linearly for better convergence
   for (k=0; ((k<MAXTRY) && (fabs(Dhi-Dlo)>eps)); k++)
     {
     EV = EVhi-(EVhi-EVlo)*Dhi/(Dhi-Dlo);
     (*b)(xlo,xhi,dx,npts,&nfwd,&nbwd,Psi,EV);
     if ( (EV>EVhi) || (EV<EVlo) ) EV = 0.5*(EVhi+EVlo);
     Dnew = DlogPsi(dx,npts,V,Psi,nfwd,nbwd,EV,&n);
     is_upper_limit = ((n>nodes) || ((n==nodes) && (Dnew < 0.0)));
     if (is_upper_limit) { EVhi = EV; Dhi = Dnew; } 
                    else { EVlo = EV; Dlo = Dnew; }
               
     }   

   return EV;
   }

/////////////////////////////////////////////////////////////////////////////

double Eigenvalue(double xlo, double xhi, double dx, int npts, 
                  int nodes, double *V, double *Psi, BOUNDARIES b)
   {
   // starting values from particle in a box
   double EVlo = -(nodes+1)*(nodes+1)*10.0/((xhi-xlo)*(xhi-xlo));
   double EVhi =  (nodes+1)*(nodes+1)*10.0/((xhi-xlo)*(xhi-xlo));
   return Eigenvalue(xlo,xhi,dx,npts,nodes,V,Psi,b,EVlo,EVhi);
   }

/////////////////////////////////////////////////////////////////////////////

double CalcPsi(double xlo, double xhi, double dx, int npts,
               double *V, double *Psi, double *DPsi, double EV, BOUNDARIES b)
{
  double h12 = dx*dx/12.0;
  int n, nfwd, nbwd; 
  
  (*b)(xlo,xhi,dx,npts,&nfwd,&nbwd,Psi,EV);

  DPsi[0] = 0.0;
  double Psimax = Psi[0];
  if (fabs(Psi[1])>Psimax) Psimax = Psi[1];
  for (n=2; n<nfwd; n++)
    {
    if (fabs(Psi[n])>Psimax) Psimax = Psi[n];
    DPsi[n-1] = 0.5*(Psi[n]-Psi[n-2])/dx;
    }
  for (n=nfwd;n<npts;n++)
    {
    Psi[n] = (   Psi[n-1]*(2.0+10.0*h12*(V[n-1]-EV)) 
               - Psi[n-2]*(1.0-     h12*(V[n-2]-EV)) )
                      /(1.0 - h12*(V[n]-EV));
    if (fabs(Psi[n])>Psimax) Psimax = Psi[n];
    DPsi[n-1] = 0.5*(Psi[n]-Psi[n-2])/dx;
    }
  DPsi[npts-1] = 0.0;
  return Psimax;  
} 

/////////////////////////////////////////////////////////////////////////////

void WritePsi(double xlo, double xhi, double dx, int npts,
         double *V, double *Psi, double EV, BOUNDARIES b,
         const char *filename, double dxout)
{
  double matchfactor, h12 = dx*dx/12.0;
  int n, nfwd, nbwd, nmatch = int(0.53*npts); 
  
  (*b)(xlo,xhi,dx,npts,&nfwd,&nbwd,Psi,EV);

  for (n=nfwd;n<=nmatch;n++)
    {
    Psi[n] = (   Psi[n-1]*(2.0+10.0*h12*(V[n-1]-EV)) 
               - Psi[n-2]*(1.0-     h12*(V[n-2]-EV)) )
             /(1.0 - h12*(V[n]-EV));
    }
  matchfactor = Psi[nmatch];

  for (n=npts-nbwd-1;n>=nmatch;n--)
    {
    Psi[n] = (   Psi[n+1]*(2.0+10.0*h12*(V[n+1]-EV)) 
               - Psi[n+2]*(1.0-     h12*(V[n+2]-EV)) )
             /(1.0 - h12*(V[n]-EV));
  }
  matchfactor /= Psi[nmatch];  
  for (n=npts-1;n>=nmatch;n--) { Psi[n] *= matchfactor; }
    
  double Psimax = 0.0;
  for (n=0;n<npts;n++) {if (Psi[n]>Psimax) Psimax=Psi[n];}  
  
  ofstream of(filename);
  double xout = xlo; 
  for (n=0;n<npts;n++)
    {
    double x = xlo+n*dx;
    Psi[n] /= Psimax;
    if (x >= xout)
      {
      of << x << '\t' << Psi[n] << '\t' << V[n] << endl;
      xout += dxout;
      }
    }
  of.close();
}

