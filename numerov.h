// $Id: numerov.h 446 2017-08-28 16:08:49Z iamp $
/****************************************************************************/
/*                                                                          */
/* Integration of the Schroedinger equation like differential equation      */
/* [d^2/dx^2 - V(x)]Psi(x) = EV*Psi(x) using the Numerov algorithm.         */
/* Boundary conditions Psi(xlo)=0.0 and Psi(xhi)=0.0 are assumed with       */
/* xlo and xhi being the end points of the integration interval             */
/*                                                                          */
/* copyright 2000 by Stefan Schippers                                       */
/****************************************************************************/


#ifndef _numerov_H
#define _numerov_H

typedef void (*BOUNDARIES)(double, double, double, int, int*, int*,
                           double*, double);
//  Procedure defining starting values for the Numerov quadrature
//  which proceeds from both ends of the integration interval up
//  to a matching point where the difference between the derivatives is
//  minimized. An example for a BOUNDARIES procedure is given below.
//  
//  void  Boundaries(double xlo, double xhi, double dx, int npts,
//                   int *nforward, int *nbackward,
//                   double *Psi, double EV)
//     {
//     starting values for the forward integration:
//     *nfwd = 2; Psi[0] = 0.0; Psi[1] = dx; 
//     starting values for the backward integration:
//     *nbwd = 2;Psi[npts-2] = dx; Psi[npts-1] = 0.0;
//     }

////////////////////////////////////////////////////////////////////////////

double Eigenvalue(double xlo, double xhi, double dx, int npts, 
                  int nodes, double *V, double *Psi,
                  BOUNDARIES b, double EVlo, double EVhi);
// returns the Eigenvalue 
// Solves the DEQ over the interval [xlo,xhi] where the potential V has to
// be defined at npts equidistant points. dx is the distance between two
// adjacent points. The wavefunction Psi with nodes nodes has to be given at
// nforwardd>=2 and nbackward>=2 points at the beginning and at the end of 
// the interval, respectively. The boundary conditions are to be calculated 
// in a subroutine of type BOUNDARIES as specified above.
// An initial serach interval [EVlo,EVhi] for the eigenvalue has to be 
// supplied.
 
////////////////////////////////////////////////////////////////////////////

double Eigenvalue(double xlo, double xhi, double dx, int npts, 
                  int nodes, double *V, double *Psi, BOUNDARIES b);
// same as above with automatic trial values EVlo and EVhi


////////////////////////////////////////////////////////////////////////////

double CalcPsi(double xlo, double xhi, double dx, int npts,
             double *V, double *Psi, double *DPSI, double EV, BOUNDARIES b);

// returns the maximum amplitude and the solution 
// Psi(x) and its derivative DPsi(x) for xlo<=x<=xhi

////////////////////////////////////////////////////////////////////////////

void WritePsi(double xlo, double xhi, double dx, int npts,
         double *V, double *Psi, double EV, BOUNDARIES b, 
         const char *filename, double dxout);
// Writes out three colums containing the coordinate x the wavefunction
// Psi(x) and the potential V(x) for all npts grid points in the interval
// [xlo,xhi]. Supply the correct eigenvalue EV, a filename and the distance
// xout between two adjacent points to be written to the file.

#endif
