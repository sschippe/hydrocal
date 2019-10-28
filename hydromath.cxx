/**
 * @file hydromath.cxx
 *
 * @brief miscellaneous mathematical functions
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: hydromath.cxx 479 2017-11-30 21:44:26Z iamp $
   @endverbatim
 *
 */
#include <stdio.h>
#include <math.h>
#include <cstring>
#include "hydromath.h"

using namespace std;

//_____________________________________________________________________________________________________________________
/**
* @brief Cubic spline interpolating function from the Numerical Recipies
*
* @param n Number of dicrete points
* @param *x Ordererd array of discrete points on the abscissa
* @param *y Corresponding function values
* @param *y2 On exit contains the second derivates of *y
* @param yp1 boundary condition, first derivative at the beginning of the range (default: natural boundary conditions)
* @param ypn boundary condition, first derivative at the end of the range (default: natural boundary conditions)
*/
void spline(int n, double *x, double*y, double *y2, double yp1, double ypn)
{
	double p,qn,sig,un;

	double *u = new double[n];
	if (yp1 > 0.99e30)
	  { // lower boundary condition is natural
	    y2[0] =0.0;
	    u[0] = 0.0;
	  }
	else 
	  { // lower boundary condition from specified derivative
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	  }

	for (int i=1; i<n-1; i++) 
	  {
      	    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
	    p=sig*y2[i-1]+2.0;
	    y2[i]=(sig-1.0)/p;
	    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
	    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	    // cout << x[i] << " " << y[i] << " " << u[i] << endl;
	  }

	if (ypn > 0.99e30)
	  { // upper boundary condition is natural
	    qn = 0.0;
	    un = 0.0;
	  }
	else 
	  { // upper boundary condition from specified first derivative
	    qn=0.5;
	    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	  }

	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);

	for (int k=n-2;k>=0;k--) 
	  { 
	    y2[k] = y2[k]*y2[k+1]+u[k]; 
	  }

	delete[] u;
}

//_____________________________________________________________________________________________________________________
/**
* @brief Cubic spline interpolation from the Numerical Recipies
*
* @param x point on the abscissa for which an interpolated function value it to be returned
* @param n Number of dicrete points of interpolating function
* @param *xa Ordered array of abscissae
* @param *ya Array of corresponding y values
* @param *y2a Array of 2nd derivatives 
*
* @return The interpolated value.
*/
double splint(double x, int n, double *xa, double *ya, double *y2a)
{
  int klo,khi,k;
  double h,b,a,y;
  
  klo=0;
  khi=n-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  
  return y;
}


double e1exp(double x)
{
// Exponential integral e^x *int_x^Infinity exp(-t)/t dt
// see Abramowitz and Stegun, Handbook of Mathematical Functions, page 231
const double a[6] = {-0.57721566, 0.99999193, -0.24991055, 0.05519968,
                     -0.00976004, 0.00107857};
const double b[5] = {0.2677737343, 8.6347608925,18.0590169730,8.5733287401,1};
const double c[5] = {3.9584969228,21.0996530827,25.6329561486,9.5733223454,1};
double pa=0,pb=0,pc=0;
int i;
if (x<1)
  {
  for (i=5;i>0;i--) { pa = (pa+a[i])*x; }
  pa += a[0]-log(x);
  return pa*exp(x);
  }
else
  {
  for (i=4;i>0;i--) { pb = (pb+b[i])*x; pc = (pc+c[i])*x; }
  pb += b[0]; pc += c[0];
  return pb/pc/x;
  }  
}

///////////////////////////////////////////////////////////////////////////

int elliptic(double m, double* K, double* E)
    {
    // Calculation of elliptic integrals by the method of the
    // arithmectic-geometric mean (see Abramowitz & Stegun,
    // "Handbook of Mathematical Functions", page 598).

    // 0 <= m < 1.0

    const double Pihalf=1.570796326794896619231322;
    const double eps = 1E-16;
    double aold=1.0, bold=sqrt(1.0-m);
    double anew, bnew, cnew=2*eps;
    double two = 1.0;
    double csum = m;

    int count = 0;
    do
      {
      anew = 0.5*(aold+bold);
      bnew = sqrt(aold*bold);
      cnew = 0.5*(aold-bold);

      count ++;      
      two *= 2.0;
      csum += two*cnew*cnew;

      aold = anew; bold = bnew; 
      }
    while (fabs(cnew)>eps);

    *K = Pihalf/anew;
    *E = *K*(1.0-0.5*csum);
    return count;
    }

///////////////////////////////////////////////////////////////////////////

int cuberoot(double a, double b, double c, double d, 
             double *x1, double *x2, double *x3)
   { // solve a*x^3 + b*x^2 + c*x + d = 0 for x
     // returns the number of real solutions
     // and the real roots found 
   *x1 = 0.0; *x2 = 0.0; *x3 = 0.0;
 
   if (a==0.0)
      { // you are left with a quadratic equation
      if (b==0.0) 
         {
         if (c==0.0) return 0;
         *x1 = -d/c;
         return 1;
         }
      double p = c/b;
      double q = d/b;
      double D = 0.25*p*p-q;
      if (D<0) return 0;
      if (D==0.0) 
         {
         *x1 = -0.5*p;
         return 1;
         }
      *x1 = -0.5*p-sqrt(D);
      *x2 = -0.5*p+sqrt(D);
      return 2;
      }  // end if (a==0.0)

   const double eps = 1e-7;
   double ya, sa, yb, sb, y1=0.0, y2=0.0, y3=0.0, swap;
   int nosol=0; // number of real solutions (return value)
   double r=b/a, s=c/a, t=d/a, r2=r*r, rthird = r/3.0;  
   double p = (3.0*s-r2)/3.0;
   double q = (2.0*r2/9.0 - s)/3.0*r + t;
   double pthird = p/3.0, q2=q*q, qhalf = 0.5*q;
   double pthird3 = pthird*pthird*pthird;
   double D = pthird3 + qhalf*qhalf;
   if (D>=0.0)
      {
      nosol = 1;
      ya = -qhalf+sqrt(D);
      sa = ya >= 0.0 ? 1.0 : -1.0; // sign of ya
      if (ya != 0.0) ya = sa*exp(log(fabs(ya))/3.0);
      //      yb = -qhalf-sqrt(D);
      //      sb = yb >= 0.0 ? 1.0 : -1.0; // sign of yb
      //      if (yb != 0.0) yb = sb*exp(log(fabs(yb))/3.0);
      yb = -p/(3.0*ya);
      y1 = ya+yb;
      if (fabs(D)<eps)
        {
        nosol = 3;
        y2 = -0.5*y1;
        y3 = -0.5*y1; 
        }
      }
   else
      {
      nosol = 3;
      const double cos2Pi3 = 0.5;
      const double sin2Pi3 = 0.8660254037844386467;
      double rho = 2.0*sqrt(-pthird);
      double cosphi = -q/rho;
      double sinphi = 1 - cosphi*cosphi;
      y1 = rho*cosphi;
      y2 = rho*( cosphi*cos2Pi3 - sin2Pi3*sinphi);
      y3 = rho*(-cosphi*cos2Pi3 + sin2Pi3*sinphi);
      }
   switch (nosol)
      {
      case 1: *x1 = y1 - rthird; 
              break;
      case 3: 
              if (y1>y3) { swap = y1; y1 = y3; y3 = swap;}
              if (y1>y2) { swap = y1; y1 = y2; y2 = swap;}
              if (y2>y3) { swap = y2; y2 = y3; y3 = swap;}
              *x1 = y1-rthird; *x2 = y2-rthird; *x3 = y3-rthird;
              break;
      } // end switch  
   return nosol;
   } 


///////////////////////////////////////////////////////////////////////////

void TestElliptic(void)
   {
   double eps=1,m,K,E;

   printf("\n\n Calculation of the ellipitc integrals K and E");
   printf("\n of the first and second kind, respectively.\n");
   
   for (;;)
     {
     printf("\n Give argument m (0 <= m <1) : ");
     scanf("%lf",&m);
     if ((m<0.0) || (m>= 1.0)) break;
     
     int count=elliptic(m,&K,&E);
     printf("\n K = %20.15f, E= %20.9f",K,E);
     printf(", %7d iterations\n\n",count);     
     }
   }

///////////////////////////////////////////////////////////////////////////

void TestDAERF()
{
  double x, h, y;

   printf("\n\n Test of differential error function DEARF(x,h)\n");
   
   for (;;)
     {
     printf("\n Give x, h : ");
     scanf("%lf %lf",&x,&h);
     y=daerf(x,h);
     printf("\n x = %10.6f, h= %10.6f, DAERF(x,y) = %20.9f ",x,h,y);
     }
}
///////////////////////////////////////////////////////////////////////////

void TestCuberoot(void)
   {
   double x1,x2,x3,a,b,c,d;
   int nosol;
   char answer[2];
   printf("\n\n Test the solver the cubic equation ");
   printf("a*x^3 + b*x^2 + c*x + d = 0");
   for(;;)
      {
      printf("\n\n Input of coeffcients or roots or stop ? (c/r/s) : ");
      scanf("%s",&answer);
      if (!strcmp(answer,"s")) break;
      if (!strcmp(answer,"c")) 
         {
         printf("\n Give coefficients a, b, c, and d : ");
         scanf("%lf %lf %lf %lf",&a,&b,&c,&d);
         }
      else
         {
         printf("\n Give roots x1, x2, x3 : ");
         scanf("%lf %lf %lf",&x1,&x2,&x3);
         a = 1.0;                   // x^3
         b = -(x1+x2+x3);           // x^2
         c = (x1*x2+x2*x3+x3*x1);   // x^1
         d = -(x1*x2*x3);           // x^0
         printf("\n equation a*x^3 + b*x^2 + c*x + d = 0 with");
         printf("\n a = %12.6g",a);
         printf("\n b = %12.6g",b);
         printf("\n c = %12.6g",c);
         printf("\n d = %12.6g\n",d);
         }  
      nosol = cuberoot(a,b,c,d,&x1,&x2,&x3);
      printf("\n %2d solutions ",nosol);
      if (nosol>=1) printf(": x1 = %10.6g",x1);
      if (nosol>=2) printf(", x2 = %10.6g",x2);
      if (nosol>=3) printf(", x3 = %10.6g",x3);
      printf("\n\n");
      }
  printf("\n\n");
  }

////////////////////////////////////////////////////////////////////////////

void TestMath(void)
   {
   int choice;
   printf("\n\n Test of mathematical subroutines");
   printf("\n  1)  Complete Elliptic Integrals of the 1st and 2nd kind");
   printf("\n  2)  Roots of a cubic equation (Cardanic formulae)");
   printf("\n  3)  Difference of error functions (DAERF)");
   printf("\n\n Make a choice : ");
   scanf("%d",&choice);
   switch (choice)
     {
     case 2  : TestCuberoot(); break;
     case 3  : TestDAERF(); break;
     default : TestElliptic(); break;
     }
   }


///////////////////////////////////////////////////////////////////////////
      int ipmpar (int i)
      {
/*
C-----------------------------------------------------------------------
C
C     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
C     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
C     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...
C
C  INTEGERS.
C
C     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM
C
C               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.
C
C     IPMPAR(1) = A, THE BASE.
C
C     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.
C
C     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C
C     IT IS ASSUMED THAT THE SINGLE AND DOUBLE FLOATING
C     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
C     NONZERO NUMBERS ARE REPRESENTED IN THE FORM
C
C               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)
C
C               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
C               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.
C
C     IPMPAR(4) = B, THE BASE.
C
C  SINGLE-
C
C     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.
C
C     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.
C
C     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-
C
C     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.
C
C     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.
C
C     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.
C
C-----------------------------------------------------------------------
C
C     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED, ACTIVATE
C     THE DATA STATMENTS FOR THE COMPUTER BY REMOVING THE C FROM
C     COLUMN 1. (ALL THE OTHER DATA STATEMENTS SHOULD HAVE C IN
C     COLUMN 1.)
C
C     IF DATA STATEMENTS ARE NOT GIVEN FOR THE COMPUTER BEING USED,
C     THEN THE FORTRAN MANUAL FOR THE COMPUTER NORMALLY GIVES THE
C     CONSTANTS IPMPAR(1), IPMPAR(2), AND IPMPAR(3) FOR THE INTEGER
C     ARITHMETIC. HOWEVER, HELP MAY BE NEEDED TO OBTAIN THE CONSTANTS
C     IPMPAR(4),...,IPMPAR(10) FOR THE SINGLE AND DOUBLE 
C     ARITHMETICS. THE SUBROUTINES MACH AND RADIX ARE PROVIDED FOR
C     THIS PURPOSE.
C
C-----------------------------------------------------------------------
C
C     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY
C     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).
C     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE
C     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.
C
C-----------------------------------------------------------------------
      INTEGER IMACH(10)
C
C     MACHINE CONSTANTS FOR THE ALLIANT FX/8.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE AMDAHL MACHINES.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /   16 /
C     DATA IMACH( 5) /    6 /
C     DATA IMACH( 6) /  -64 /
C     DATA IMACH( 7) /   63 /
C     DATA IMACH( 8) /   14 /
C     DATA IMACH( 9) /  -64 /
C     DATA IMACH(10) /   63 /
C
C     MACHINE CONSTANTS FOR THE AT&T 3B SERIES, AT&T
C     PC 7300, AND AT&T 6300.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   33 /
C     DATA IMACH( 3) / 8589934591 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   24 /
C     DATA IMACH( 6) / -256 /
C     DATA IMACH( 7) /  255 /
C     DATA IMACH( 8) /   60 /
C     DATA IMACH( 9) / -256 /
C     DATA IMACH(10) /  255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   39 /
C     DATA IMACH( 3) / 549755813887 /
C     DATA IMACH( 4) /    8 /
C     DATA IMACH( 5) /   13 /
C     DATA IMACH( 6) /  -50 /
C     DATA IMACH( 7) /   76 /
C     DATA IMACH( 8) /   26 /
C     DATA IMACH( 9) /  -50 /
C     DATA IMACH(10) /   76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C     DATA IMACH( 1) /      2 /
C     DATA IMACH( 2) /     39 /
C     DATA IMACH( 3) / 549755813887 /
C     DATA IMACH( 4) /      8 /
C     DATA IMACH( 5) /     13 /
C     DATA IMACH( 6) /    -50 /
C     DATA IMACH( 7) /     76 /
C     DATA IMACH( 8) /     26 /
C     DATA IMACH( 9) / -32754 /
C     DATA IMACH(10) /  32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C     COMPUTERS, AND THE CDC CYBER 990 AND 995 (NOS
C     OPERATING SYSTEM).
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   48 /
C     DATA IMACH( 3) / 281474976710655 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   48 /
C     DATA IMACH( 6) / -974 /
C     DATA IMACH( 7) / 1070 /
C     DATA IMACH( 8) /   95 /
C     DATA IMACH( 9) / -926 /
C     DATA IMACH(10) / 1070 /
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 990 AND 995
C     (NOS/VE OPERATING SYSTEM).
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    63 /
C     DATA IMACH( 3) / 9223372036854775807 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    48 /
C     DATA IMACH( 6) / -4096 /
C     DATA IMACH( 7) /  4095 /
C     DATA IMACH( 8) /    96 /
C     DATA IMACH( 9) / -4096 /
C     DATA IMACH(10) /  4095 /
C
C     MACHINE CONSTANTS FOR THE CONVEX COMPUTERS
C     (NATIVE MODE).
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -127 /
C     DATA IMACH( 7) /   127 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1023 /
C     DATA IMACH(10) /  1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX COMPUTERS
C     (IEEE MODE).
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE CRAY 2, X-MP, AND Y-MP
C     (CFT77 COMPILER USING THE 64 BIT INTEGER ARITHMETIC).
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    63 /
C     DATA IMACH( 3) / 9223372036854775807 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    48 /
C     DATA IMACH( 6) / -8188 /
C     DATA IMACH( 7) /  8189 /
C     DATA IMACH( 8) /    96 /
C     DATA IMACH( 9) / -8188 /
C     DATA IMACH(10) /  8189 /
C
C     MACHINE CONSTANTS FOR THE CRAY 2, X-MP, AND Y-MP
C     (CFT77 COMPILER USING THE 46 BIT INTEGER ARITHMETIC).
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    46 /
C     DATA IMACH( 3) / 70368744177663 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    48 /
C     DATA IMACH( 6) / -8188 /
C     DATA IMACH( 7) /  8189 /
C     DATA IMACH( 8) /    96 /
C     DATA IMACH( 9) / -8188 /
C     DATA IMACH(10) /  8189 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   15 /
C     DATA IMACH( 3) / 32767 /
C     DATA IMACH( 4) /   16 /
C     DATA IMACH( 5) /    6 /
C     DATA IMACH( 6) /  -64 /
C     DATA IMACH( 7) /   63 /
C     DATA IMACH( 8) /   14 /
C     DATA IMACH( 9) /  -64 /
C     DATA IMACH(10) /   63 /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   23 /
C     DATA IMACH( 3) / 8388607 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   23 /
C     DATA IMACH( 6) / -127 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   38 /
C     DATA IMACH( 9) / -127 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000
C     AND DPS 8/70 SERIES.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   35 /
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   27 /
C     DATA IMACH( 6) / -127 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   63 /
C     DATA IMACH( 9) / -127 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 9000.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -126 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE IBM 3033, THE ICL 2900, THE ITEL AS/6, THE
C     XEROX SIGMA 5/7/9, AND THE SEL SYSTEMS 85/86.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /   16 /
C     DATA IMACH( 5) /    6 /
C     DATA IMACH( 6) /  -64 /
C     DATA IMACH( 7) /   63 /
C     DATA IMACH( 8) /   14 /
C     DATA IMACH( 9) /  -64 /
C     DATA IMACH(10) /   63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT
C     MACFORTRAN II.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   24 /
C     DATA IMACH( 6) / -127 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   56 /
C     DATA IMACH( 9) / -127 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS-4D
C     SERIES (MIPS R3000 PROCESSOR).
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE SUN 3.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    35 /
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    27 /
C     DATA IMACH( 6) /  -128 /
C     DATA IMACH( 7) /   127 /
C     DATA IMACH( 8) /    60 /
C     DATA IMACH( 9) / -1024 /
C     DATA IMACH(10) /  1023 /
C
C     MACHINE CONSTANTS FOR THE VAX AND MICROVAX
C     COMPUTERS - F AND D FLOATING ARITHMETICS.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -127 /
C     DATA IMACH( 7) /   127 /
C     DATA IMACH( 8) /    56 /
C     DATA IMACH( 9) /  -127 /
C     DATA IMACH(10) /   127 /
C
C     MACHINE CONSTANTS FOR THE VAX AND MICROVAX
C     COMPUTERS - F AND G FLOATING ARITHMETICS.
C
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -127 /
C     DATA IMACH( 7) /   127 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1023 /
C     DATA IMACH(10) /  1023 /
C
*/
      const int imach[11] = {         0, // dummy
                                      2, // 1: A (base for integer)
                                     31, // 2: N (number of integer digits)
                             2147483647, // 3: A**N-1 (largest integer)
                                      2, // 4: B (base for floating point) 
                                     24, // 5: number of base B digits
                                   -125, // 6: smallest exponent single prec.
                                    127, // 7: largest exponent single prec. 
                                     53, // 6: number of base B digits 
                                  -1021, // 8: smallest exponent double prec.
                                   1023};// 9: largest exponent double prec.
      return imach[i];
      }

///////////////////////////////////////////////////////////////////////////

      double dpmpar (int i)
      {
/*
C-----------------------------------------------------------------------
C
C     DPMPAR PROVIDES THE DOUBLE PRECISION MACHINE CONSTANTS FOR
C     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
C     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
C     DOUBLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
C     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN
C
C        DPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,
C
C        DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
C
C        DPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
C
C-----------------------------------------------------------------------
C     WRITTEN BY
C        ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WARFARE CENTER
C        DAHLGREN VIRGINIA
C-----------------------------------------------------------------------
*/
      int emin, emax;
      double b, binv, bm1, one, w, z, m;

      if (i == 0)
         { 
         b = ipmpar(4);
         m = ipmpar(8);
         return pow(b,1 - m);
         }
     else if (i == 1)
         {
         b = ipmpar(4);
         emin = ipmpar(9);
         binv = 1.0/b;
         w = pow(b,emin + 2);
         return ((w * binv) * binv) * binv;
         }
      else
         {
         int ibeta = ipmpar(4);
         m = ipmpar(8);
         emax = ipmpar(10);

         b = ibeta;
         bm1 = ibeta - 1;
         z = pow(b,m - 1);
         w = ((z - 1.0)*b + bm1)/(b*z);

         z = pow(b,emax - 2);
         return ((w * z) * b) * b;
         }
      }

///////////////////////////////////////////////////////////////////////////

      double depsln (int l)
      {
/*
C--------------------------------------------------------------------
C     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
C     SUCH THAT 1.D0 + EPS .GT. 1.D0 .  L IS A DUMMY ARGUMENT.
C--------------------------------------------------------------------
*/
      int b;
      double db, lnb, m;

      b = ipmpar(4);
      if (b == 2) {lnb = .693147180559945309417232121458;}
      else if (b == 8) {lnb = 2.07944154167983592825169636437;}
      else if (b == 16) {lnb = 2.77258872223978123766892848583;}
      else {db = b; lnb = log(db);}

      m = 1.0 - ipmpar(8);
      return m * lnb;
      }

      double dxparg (int l)
      {
/*
C--------------------------------------------------------------------
C     IF L = 0 THEN  DXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
C     DEXP(W) CAN BE COMPUTED.
C
C     IF L IS NONZERO THEN  DXPARG(L) = THE LARGEST NEGATIVE W FOR
C     WHICH THE COMPUTED VALUE OF DEXP(W) IS NONZERO.
C
C     NOTE... ONLY AN APPROXIMATE VALUE FOR DXPARG(L) IS NEEDED.
C--------------------------------------------------------------------
*/
      int b;
      double db, lnb, m ;

      b = ipmpar(4);
      if (b == 2) {lnb = .693147180559945309417232121458;}
      else if (b == 8) {lnb = 2.07944154167983592825169636437;}
      else if (b == 16) {lnb = 2.77258872223978123766892848583;}
      else {db = b; lnb = log(db);}

      m = l==0  ? ipmpar(10) : ipmpar(9) - 1;
      return  0.999999999999 * m * lnb;
      }

///////////////////////////////////////////////////////////////////////////

      double derfc0 (double x)
      {
/*
C-----------------------------------------------------------------------
C
C           EVALUATION OF EXP(X**2)*ERFC(X) FOR X .GE. 1
C
C-----------------------------------------------------------------------
C     WRITTEN BY ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WARFARE CENTER
C        DAHLGREN, VIRGINIA
C        APRIL 1992
C-------------------------------
*/
      const double rpinv = 0.56418958354775628694807945156077259;  // = 1/sqrt(Pi)

      const double p0 = .16506148041280876191828601e-03;
      const double p1 = .15471455377139313353998665e-03;
      const double p2 = .44852548090298868465196794e-04;
      const double p3 =-.49177280017226285450486205e-05;
      const double p4 =-.69353602078656412367801676e-05;
      const double p5 =-.20508667787746282746857743e-05;
      const double p6 =-.28982842617824971177267380e-06;
      const double p7 =-.17272433544836633301127174e-07;

      const double q1 = .16272656776533322859856317e+01;
      const double q2 = .12040996037066026106794322e+01;
      const double q3 = .52400246352158386907601472e+00;
      const double q4 = .14497345252798672362384241e+00;
      const double q5 = .25592517111042546492590736e-01;
      const double q6 = .26869088293991371028123158e-02;
      const double q7 = .13133767840925681614496481e-03;

      const double r0 = .145589721275038539045668824025e+00;
      const double r1 =-.273421931495426482902320421863e+00;
      const double r2 = .226008066916621506788789064272e+00;
      const double r3 =-.163571895523923805648814425592e+00;
      const double r4 = .102604312032193978662297299832e+00;
      const double r5 =-.548023266949835519254211506880e-01;
      const double r6 = .241432239725390106956523668160e-01;
      const double r7 =-.822062115403915116036874169600e-02;
      const double r8 = .180296241564687154310619200000e-02;

      const double a0 =-.45894433406309678202825375e-03;
      const double a1 =-.12281298722544724287816236e-01;
      const double a2 =-.91144359512342900801764781e-01;
      const double a3 =-.28412489223839285652511367e-01;
      const double a4 = .14083827189977123530129812e+01;
      const double a5 = .11532175281537044570477189e+01;
      const double a6 =-.72170903389442152112483632e+01;
      const double a7 =-.19685597805218214001309225e+01;
      const double a8 = .93846891504541841150916038e+01;

      const double b1 = .25136329960926527692263725e+02;
      const double b2 = .15349442087145759184067981e+03;
      const double b3 =-.29971215958498680905476402e+03;
      const double b4 =-.33876477506888115226730368e+04;
      const double b5 = .28301829314924804988873701e+04;
      const double b6 = .22979620942196507068034887e+05;
      const double b7 =-.24280681522998071562462041e+05;
      const double b8 =-.36680620673264731899504580e+05;
      const double b9 = .42278731622295627627042436e+05;
      const double b10= .28834257644413614344549790e+03;
      const double b11= .70226293775648358646587341e+03;

      const double c0 =-.7040906288250128001000086e-04;
      const double c1 =-.3858822461760510359506941e-02;
      const double c2 =-.7708202127512212359395078e-01;
      const double c3 =-.6713655014557429480440263e+00;
      const double c4 =-.2081992124162995545731882e+01;
      const double c5 = .2898831421475282558867888e+01;
      const double c6 = .2199509380600429331650192e+02;
      const double c7 = .2907064664404115316722996e+01;
      const double c8 =-.4766208741588182425380950e+02;
 
      const double d1 = .5238852785508439144747174e+02;
      const double d2 = .9646843357714742409535148e+03;
      const double d3 = .7007152775135939601804416e+04;
      const double d4 = .8515386792259821780601162e+04;
      const double d5 =-.1002360095177164564992134e+06;
      const double d6 =-.2065250031331232815791912e+06;
      const double d7 = .5695324805290370358175984e+06;
      const double d8 = .6589752493461331195697873e+06;
      const double d9 =-.1192930193156561957631462e+07;

      const double e0 = .540464821348814822409610122136e+00;
      const double e1 =-.261515522487415653487049835220e-01;
      const double e2 =-.288573438386338758794591212600e-02;
      const double e3 =-.529353396945788057720258856000e-03;

//    COEFFICIENTS FOR THE ASYMPTOTIC EXPANSION

      const double s1 = .75000000000000000000e+00;
      const double s2 =-.18750000000000000000e+01;
      const double s3 = .65625000000000000000e+01;
      const double s4 =-.29531250000000000000e+02;
      const double s5 = .16242187500000000000e+03;
      const double s6 =-.10557421875000000000e+04;
      const double s7 = .79180664062500000000e+04;
      const double s8 =-.67303564453125000000e+05;
      const double s9 = .63938386230468750000e+06;
      const double s10=-.67135305541992187500e+07;
      const double s11= .77205601373291015625e+08;

      double  t, u, v, z; 

      if (x>=1.0 && x<= 2.0)
         {
         u = ((((((p7*x + p6)*x + p5)*x + p4)*x + p3)*x + p2)*x + p1)*x + p0;
         v = ((((((q7*x + q6)*x + q5)*x + q4)*x + q3)*x + q2)*x + q1)*x + 1.0;
         t = (x - 3.75)/(x + 3.75);
         return (((((((((u/v)*t + r8)*t + r7)*t + r6)*t + r5)*t + r4)*t + r3)*t + r2)*t + r1)*t + r0;
         }
      else if (x>2 && x<= 4.0)
         {
         z = 1.0/(2.5 + x*x);
         u = (((((((a8*z + a7)*z + a6)*z + a5)*z + a4)*z + a3)*z + a2)*z + a1)*z + a0;
         v = ((((((((((b11*z + b10)*z + b9)*z + b8)*z + b7)*z + b6)*z + b5)*z + b4)*z + b3)*z + b2)*z + b1)*z + 1.0;
         t = 13.0*z - 1.0;
         return ((((u/v)*t + e2)*t + e1)*t + e0)/x;
         }
      else if (x>4 && x<= 50.0)
         {
         z = 1.0/(2.50 + x*x);
         u = (((((((c8*z + c7)*z + c6)*z + c5)*z + c4)*z + c3)*z + c2)*z + c1)*z + c0;
         v = ((((((((d9*z + d8)*z + d7)*z + d6)*z + d5)*z + d4)*z + d3)*z + d2)*z + d1)*z + 1.0;
         t = 13.0*z - 1.0;
         return (((((u/v)*t + e3)*t + e2)*t + e1)*t + e0)/x;
         }
      else
         {
         t = 1.0/x/x;
         z = (((((((((((s11*t + s10)*t + s9)*t + s8)*t + s7)*t + s6)*t + s5)*t + s4)*t + s3)*t + s2)*t + s1)*t -0.5)*t + 1.0;
         return rpinv*(z/x);
         }
      }
      
///////////////////////////////////////////////////////////////////////////

      double derf (double x)
      {
/*
C-----------------------------------------------------------------------
C        DOUBLE EVALUATION OF THE ERROR FUNCTION
C-----------------------------------------------------------------------
*/
      double ax, t, w;
      double a[21] = { .1283791670955125738961589031215E+00,
                      -.3761263890318375246320529677070E+00,
                       .1128379167095512573896158902931E+00,
                      -.2686617064513125175943235372542E-01,
                       .5223977625442187842111812447877E-02,
                      -.8548327023450852832540164081187E-03,
                       .1205533298178966425020717182498E-03,
                      -.1492565035840625090430728526820E-04,
                       .1646211436588924261080723578109E-05,
                      -.1636584469123468757408968429674E-06,
                       .1480719281587021715400818627811E-07,
                      -.1229055530145120140800510155331E-08,
                       .9422759058437197017313055084212E-10,
                      -.6711366740969385085896257227159E-11,
                       .4463222608295664017461758843550E-12,
                      -.2783497395542995487275065856998E-13,
                       .1634095572365337143933023780777E-14,
                      -.9052845786901123985710019387938E-16,
                       .4708274559689744439341671426731E-17,
                      -.2187159356685015949749948252160E-18,
                       .7043407712019701609635599701333E-20};

      ax = fabs(x);
      if (ax < 1.0)
         {
         t = x*x;
         w = a[20];
         for (int i=1; i<=20; i++)
            {
            w = t*w + a[20-i];
            }
         return  x*(1.0 + w);
         }
      else if (ax < 8.5)
         {
         double result = 0.5 + (0.5  - exp(-x*x)*derfc0(ax));
         return x<0 ? -result : result;
         }
      else
         {      //  limit value for large x
         return x<0 ? -1.0 : 1.0;
         }
      }

///////////////////////////////////////////////////////////////////////////

      double derfc (double x)
      {
/*
C-----------------------------------------------------------------------
C         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION
C-----------------------------------------------------------------------
*/
      double ax, t, w;
      double a[21] = {  .1283791670955125738961589031215e+00,
                       -.3761263890318375246320529677070e+00,
                        .1128379167095512573896158902931e+00,
                       -.2686617064513125175943235372542e-01,
                        .5223977625442187842111812447877e-02,
                       -.8548327023450852832540164081187e-03,
                        .1205533298178966425020717182498e-03,
                       -.1492565035840625090430728526820e-04,
                        .1646211436588924261080723578109e-05,
                       -.1636584469123468757408968429674e-06,
                        .1480719281587021715400818627811e-07,
                       -.1229055530145120140800510155331e-08,
                        .9422759058437197017313055084212e-10,
                       -.6711366740969385085896257227159e-11,
                        .4463222608295664017461758843550e-12,
                       -.2783497395542995487275065856998e-13,
                        .1634095572365337143933023780777e-14,
                       -.9052845786901123985710019387938e-16,
                        .4708274559689744439341671426731e-17,
                       -.2187159356685015949749948252160e-18,
                        .7043407712019701609635599701333e-20};

      ax = fabs(x);
      t = x*x;
      if (ax < 1.0) 
         {
         w = a[20];
         for (int i=1; i<=20; i++)
            {
            w = t*w + a[20-i];
            }
         return  0.5 + (0.5 - x*(1.0 + w));
         }
      else if (x <-8.3)      {return 2.0;}
      else if (x < 0.0)      {return 2.0 - exp(-t)*derfc0(ax);}
      else if (t<-dxparg(1)) {return exp(-t)*derfc0(x);}
      else                   {return 0.0;}
      
      }

      
///////////////////////////////////////////////////////////////////////////

      double daerf (double x, double h)
      {
/*
C-----------------------------------------------------------------------
C             COMPUTATION OF ERF(X + H) - ERF(X - H)
C-----------------------------------------------------------------------
*/
      if (h==0.0) {return 0.0;}

      double ah, ax, a2, dn2, e, eps, hf, hg, h2, h3, n, n1; 
      double s, st, t, u, v, xmh, xph, x2, z;

      const double c = 1.12837916709551257389615890312155; // 2/sqrt(Pi)
      const double p = 2.76959; // ln(9*sqrt(Pi))

/*
C     **** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
C          SMALLEST NUMBER SUCH THAT 1.D0 + EPS .GT. 1.D0 .
*/
      eps = dpmpar(1);

      ah = fabs(h);
      ax = fabs(x);
      xph = ax + ah;
      xmh = ax - ah;

      t = ah > ax ? ah : ax;
      t = t*t;

      if (1.6*t*t < eps) {return 2.0*c*h*(0.5 + (0.5 - (x*x + h*h/3.0)));}

      if (pow(ax*ah,2) < 0.5*eps) 
         {
         t = 2.0;
         x2 = x*x;
         if (x2 >= eps) t = 2.0*exp(-x2);
         if (h*h >= 3.0*eps) {return t*derf(h);}
            else {return c*h*t;}
         }
    
      if (ax <= ah) 
        {
        if (xph <  8.5) {return h>0 ? derf(xph)-derf(xmh) : derf(xmh)-derf(xph);}
        else if (xmh > -8.3) {return h>0 ? derfc(xmh) : -derfc(xmh);}
        else {return h>0 ? 2.0 : -2.0;}
        }

      if ((xmh >= 9.0) && (xmh*xmh + p > -dxparg(1))) {return 0.0;}

      if (4.0*ah*ax > -depsln(0)) {return h<0 ? -derfc(xmh) : derfc(xmh);}

      if (ax < 3.0*ah)
         {
         if (xph <  1.0) {return h<0 ? derf(xmh)-derf(xph) : derf(xph)-derf(xmh);}
         else { return h>0 ? derfc(xmh)-derfc(xph) : derfc(xph)-derfc(xmh);}
         } 

      e = eps > 1e-30 ? eps : 1e-30;
      if (ax <= 0.4)
         {
         h2 = xph*xph;
         a2 = xmh*xmh;
         x2 = ax + ax;
         st = 1.0;
         hf = xmh;
         n = 0.0;
         n1 = 1.0;
         dn2 = 1.0;
         s = 0.0;
         t = 1.0;
         while (fabs(t) > e*fabs(s))
           {
           n = n + 1.0;
           n1 = n1 + 2.0;
           dn2 = -dn2/n;
           st = h2*st + x2*hf;
           hf = a2*hf;
           t = st*dn2/n1;
           s = s + t;
           }
         s = 0.5 + (0.5 + s);
         return h>0 ? 2.0*c*ah*s : -2.0*c*ah*s;
         }
      else
         {
         n = 1.0;
         int j = 2;
         h2 = 0.0;
         z = exp(-0.5*ax*ax);
         u = 2.0*ah*c*z;
         h3 = z;
         v = 2.0*h*h;
         hf = 2.0*ax*ah;
         s = 0.0;
         while (j)
           {
           hg = 1e9;
           while (fabs(hg) > e*fabs(s))
             {
             h2 = (hf*h3 - v*h2)/n;
             n = n + 1.0;
             h3 = (hf*h2 - v*h3)/n;
             n = n + 1.0;
             hg = h3/n;
             s = s+hg;
             }
           j--;
           }
         return h>0 ? u*(s+z) : -u*(s+z);
         }
      }

///////////////////////////////////////////////////////////////////////////

const double lnfactorial[nmaxfactorial+1]={0.0,
0., 0.6931471805599453, 1.791759469228054, 3.178053830347945, 4.787491742782045, 6.579251212010101, 
8.52516136106541, 10.60460290274525, 12.80182748008147, 15.10441257307551, 17.50230784587388, 19.98721449566188, 
22.55216385312342, 25.19122118273868, 27.89927138384089, 30.67186010608067, 33.50507345013689, 36.39544520803305, 
39.33988418719949, 42.33561646075348, 45.3801388984769, 48.47118135183523, 51.60667556776437, 54.78472939811232, 
58.00360522298052, 61.26170176100201, 64.55753862700634, 67.88974313718154, 71.25703896716801, 74.65823634883017, 
78.09222355331532, 81.5579594561151, 85.0544670175815, 88.5808275421977, 92.1361756036871, 95.7196945421432, 
99.3306124547874, 102.9681986145138, 106.6317602606434, 110.3206397147574, 114.0342117814617, 117.771881399745, 
121.5330815154386, 125.3172711493568, 129.1239336391272, 132.9525750356163, 136.8027226373263, 140.6739236482342, 
144.5657439463448, 148.477766951773, 152.4095925844973, 156.3608363030787, 160.3311282166309, 164.3201122631952, 
168.3274454484276, 172.3527971391628, 176.3958484069973, 180.4562914175437, 184.5338288614495, 188.6281734236715, 
192.7390472878449, 196.8661816728899, 201.0093163992815, 205.1681994826412, 209.3425867525368, 213.5322414945632, 
217.7369341139542, 221.9564418191303, 226.1905483237276, 230.4390435657769, 234.7017234428182, 238.9783895618343, 
243.2688490029827, 247.5729140961868, 251.8904022097232, 256.2211355500095, 260.5649409718632, 264.9216497985528, 
269.2910976510198, 273.6731242856937, 278.0675734403661, 282.4742926876304, 286.893133295427, 291.3239500942703, 
295.7666013507606, 300.2209486470141, 304.6868567656686, 309.1641935801469, 313.6528299498791, 318.1526396202093, 
322.6634991267262, 327.1852877037752, 331.7178871969285, 336.2611819791985, 340.815058870799, 345.3794070622669, 
349.9541180407702, 354.5390855194408, 359.1342053695755, 363.7393755555636, 368.3544960724048, 372.9794688856891, 
377.6141978739187, 382.2585887730601, 386.9125491232177, 391.5759882173297, 396.2488170517916, 400.9309482789158, 
405.622296161145, 410.3227765269374, 415.0323067282498, 419.7508055995449, 424.4781934182571, 429.2143918666516, 
433.9593239950149, 438.7129141861213, 443.475088120919, 448.2457727453846, 453.0248962384962, 457.8123879812783, 
462.608178526875, 467.4121995716082, 472.2243839269806, 477.0446654925857, 481.872979229888, 486.7092611368395, 
491.553448223298, 496.4054784872176, 501.2652908915794, 506.132825342035, 511.0080226652361, 515.8908245878224, 
520.7811737160442, 525.6790135159951, 530.5842882944335, 535.4969431801696, 540.4169241059976, 545.3441777911549, 
550.2786517242856, 555.220294146895, 560.1690540372731, 565.1248810948743, 570.0877257251343, 575.0575390247102, 
580.0342727671308, 585.0178793888392, 590.0083119756179, 595.0055242493821, 600.0094705553275, 605.0201058494237, 
610.0373856862387, 615.061266207085, 620.0917041284775, 625.1286567308911, 630.1720818478104, 635.2219378550599, 
640.2781836604082, 645.3407786934351, 650.4096828956553, 655.4848567108892, 660.5662610758737, 665.6538574111062, 
670.747607611913, 675.847474039737, 680.9534195136376, 686.0654073019941, 691.1834011144108, 696.3073650938142, 
701.4372638087372, 706.5730622457875, 711.7147258022901, 716.8622202791036, 722.0155118736014, 727.1745671728159, 
732.3393531467393, 737.5098371417775, 742.6859868743513, 747.8677704246433, 753.0551562304841, 758.2481130813744, 
763.4466101126402, 768.650616799717, 773.8601029525584, 779.0750387101674, 784.2953945352457, 789.5211412089588, 
794.7522498258135, 799.9886917886433, 805.230438803703, 810.477462875864, 815.72973630391, 820.987231675938, 
826.249921864843, 831.517780023906, 836.79077958247, 842.0688942417, 847.352097970438, 852.640365001133, 
857.933669825857, 863.231987192405, 868.535292100465, 873.843559797866, 879.156765776908, 884.474885770752, 
889.79789574989, 895.12577191868, 900.458490711945, 905.796028791646, 911.138363043611, 916.485470574329, 
921.837328707805, 927.193914982477, 932.555207148186, 937.921183163208, 943.291821191336, 948.66709959902, 
954.04699695256, 959.43149201535, 964.820563745166, 970.214191291518, 975.612353993036, 981.015031374908, 
986.422203146369, 991.833849198224, 997.249949600428, 1002.6704845997, 1008.095434617181, 1013.524780246136, 
1018.95850224969, 1024.396581558613, 1029.838999269135, 1035.285736640801, 1040.736775094367, 1046.192096209725, 
1051.651681723869, 1057.115513528894, 1062.58357367003, 1068.055844343701, 1073.532307895633, 1079.012946818975, 
1084.497743752465, 1089.986681478622, 1095.479742921963, 1100.976911147256, 1106.478169357801, 1111.983500893733, 
1117.492889230361, 1123.006317976526, 1128.523770872991, 1134.045231790853, 1139.570684729985, 1145.100113817496, 
1150.633503306224, 1156.170837573242, 1161.712101118401, 1167.25727856288, 1172.806354647776, 1178.359314232697, 
1183.916142294397, 1189.476823925412, 1195.041344332735, 1200.609688836496, 1206.181842868674, 1211.75779197182, 
1217.337521797807, 1222.921018106588, 1228.508266764988, 1234.0992537455, 1239.693965125101, 1245.2923870841, 
1250.894505904979, 1256.500307971275, 1262.109779766461, 1267.722907872849, 1273.339678970515, 1278.960079836232, 
1284.58409734242, 1290.21171845611, 1295.842930237931, 1301.477719841101, 1307.116074510435, 1312.757981581373, 
1318.403428479016, 1324.052402717177, 1329.704891897446, 1335.360883708265, 1341.020365924025, 1346.683326404161, 
1352.349753092274, 1358.019634015254, 1363.692957282426, 1369.369711084694, 1375.049883693711, 1380.73346346105, 
1386.420438817389, 1392.110798271714, 1397.804530410516, 1403.501623897022, 1409.202067470412, 1414.905849945069, 
1420.612960209818, 1426.323387227192, 1432.037120032702, 1437.754147734108, 1443.474459510716, 1449.198044612668, 
1454.924892360255, 1460.654992143229, 1466.388333420127, 1472.124905717606, 1477.864698629785, 1483.607701817595, 
1489.353905008135, 1495.103297994043, 1500.855870632869, 1506.611612846456, 1512.370514620333, 1518.132566003113, 
1523.897757105898, 1529.666078101692, 1535.437519224822, 1541.212070770366, 1546.989723093589, 1552.770466609382, 
1558.554291791711, 1564.341189173078, 1570.131149343975, 1575.924162952359, 1581.720220703125, 1587.519313357585, 
1593.321431732962, 1599.126566701879, 1604.934709191859, 1610.745850184836, 1616.559980716661, 1622.377091876624, 
1628.197174806977, 1634.020220702459, 1639.84622080984, 1645.67516642745, 1651.507048904734, 1657.341859641797, 
1663.179590088962, 1669.020231746336, 1674.863776163367, 1680.710214938425, 1686.559539718372, 1692.411742198146, 
1698.266814120349, 1704.124747274832, 1709.985533498298, 1715.849164673896, 1721.715632730829, 1727.584929643963, 
1733.457047433438, 1739.33197816429, 1745.20971394607, 1751.090246932471, 1756.973569320959, 1762.859673352409, 
1768.748551310742, 1774.640195522568, 1780.534598356832, 1786.431752224469, 1792.331649578051, 1798.234282911453, 
1804.139644759507, 1810.047727697676, 1815.958524341717, 1821.872027347355, 1827.788229409963, 1833.707123264236, 
1839.62870168388, 1845.552957481294, 1851.479883507265, 1857.409472650655, 1863.341717838103, 1869.276612033722, 
1875.214148238804, 1881.154319491525, 1887.097118866651, 1893.042539475258, 1898.990574464439, 1904.941217017026, 
1910.894460351314, 1916.850297720779, 1922.808722413809, 1928.769727753432, 1934.733307097051, 1940.699453836174, 
1946.66816139616, 1952.63942323595, 1958.613232847819, 1964.589583757117, 1970.568469522018, 1976.549883733273, 
1982.53382001396, 1988.520272019244, 1994.509233436135, 2000.500697983242, 2006.49465941055, 2012.491111499168, 
2018.490048061115, 2024.491462939076, 2030.495350006183, 2036.501703165785, 2042.510516351227, 2048.521783525631, 
2054.535498681674, 2060.551655841373, 2066.570249055869, 2072.591272405218, 2078.614719998179, 2084.640585972004, 
2090.668864492235, 2096.699549752496, 2102.732635974295, 2108.76811740682, 2114.805988326742, 2120.84624303802, 
2126.888875871702, 2132.933881185739, 2138.981253364784, 2145.030986820017, 2151.083075988941, 2157.13751533521, 
2163.194299348439, 2169.25342254402, 2175.314879462948, 2181.378664671636, 2187.44477276174, 2193.513198349983, 
2199.583936077986, 2205.656980612086, 2211.732326643176, 2217.809968886525, 2223.88990208162, 2229.972120991997, 
2236.056620405072, 2242.143395131984, 2248.232440007431, 2254.323749889509, 2260.417319659554, 2266.513144221986, 
2272.611218504152, 2278.711537456172, 2284.814096050786, 2290.918889283201, 2297.025912170943, 2303.135159753709, 
2309.246627093211, 2315.360309273043, 2321.476201398526, 2327.594298596567, 2333.714596015519, 2339.837088825033, 
2345.961772215927, 2352.088641400041, 2358.217691610101, 2364.348918099584, 2370.482316142582, 2376.617881033663, 
2382.755608087749, 2388.895492639975, 2395.037530045563, 2401.181715679689, 2407.328044937358, 2413.476513233276, 
2419.627116001722, 2425.779848696426, 2431.934706790442, 2438.091685776028, 2444.25078116452, 2450.411988486215, 
2456.57530329025, 2462.740721144481, 2468.908237635369, 2475.077848367861, 2481.249548965272, 2487.423335069174, 
2493.59920233928, 2499.77714645333, 2505.957163106983, 2512.1392480137, 2518.323396904637, 2524.509605528538, 
2530.69786965162, 2536.888185057474, 2543.080547546948, 2549.274952938053, 2555.471397065848, 2561.669875782341, 
2567.870384956383, 2574.072920473571, 2580.27747823614, 2586.484054162864, 2592.692644188961, 2598.903244265985, 
2605.115850361737, 2611.33045846016, 2617.547064561244, 2623.765664680936, 2629.986254851035, 2636.208831119107, 
2642.433389548382, 2648.65992621767, 2654.888437221261, 2661.118918668839, 2667.35136668539, 2673.585777411109, 
2679.822147001313, 2686.060471626352, 2692.300747471523, 2698.542970736978, 2704.787137637642, 2711.033244403123, 
2717.281287277632, 2723.531262519892, 2729.783166403058, 2736.036995214633, 2742.292745256387, 2748.550412844269, 
2754.809994308335, 2761.071485992655, 2767.334884255247, 2773.600185467985, 2779.867386016526, 2786.136482300232, 
2792.40747073209, 2798.680347738637, 2804.955109759878, 2811.23175324922, 2817.510274673385, 2823.790670512345, 
2830.072937259241, 2836.357071420312, 2842.643069514821, 2848.930928074982, 2855.220643645892, 2861.512212785449, 
2867.805632064296, 2874.100898065736, 2880.39800738567, 2886.696956632526, 2892.997742427189, 2899.300361402934, 
2905.604810205356, 2911.911085492304, 2918.219183933814, 2924.52910221204, 2930.840837021193, 2937.15438506747, 
2943.469743068992, 2949.78690775574, 2956.105875869486, 2962.426644163737, 2968.749209403664, 2975.073568366045, 
2981.3997178392, 2987.72765462293, 2994.057375528452, 3000.388877378346, 3006.722157006486, 3013.057211257984, 
3019.39403698913, 3025.732631067334, 3032.072990371061, 3038.415111789782, 3044.758992223909, 3051.104628584737, 
3057.452017794393, 3063.801156785772, 3070.152042502487, 3076.504671898807, 3082.859041939604, 3089.2151496003, 
3095.572991866808, 3101.93256573548, 3108.293868213054, 3114.656896316593, 3121.021647073446, 3127.388117521177, 
3133.756304707528, 3140.126205690356, 3146.497817537588, 3152.871137327165, 3159.246162146993, 3165.622889094891, 
3172.001315278543, 3178.381437815443, 3184.763253832849, 3191.146760467733, 3197.53195486673, 3203.918834186093, 
3210.307395591639, 3216.697636258705, 3223.089553372097, 3229.483144126048, 3235.878405724163, 3242.27533537938, 
3248.673930313915, 3255.074187759224, 3261.476104955951, 3267.879679153885, 3274.284907611916, 3280.691787597985, 
3287.100316389044, 3293.510491271011, 3299.92230953872, 3306.335768495888, 3312.750865455059, 3319.167597737572, 
3325.585962673508, 3332.005957601655, 3338.427579869461, 3344.850826832995, 3351.2756958569, 3357.702184314358, 
3364.130289587043, 3370.560009065082, 3376.991340147016, 3383.424280239754, 3389.858826758542, 3396.294977126912, 
3402.732728776648, 3409.172079147748, 3415.613025688382, 3422.055565854849, 3428.49969711155, 3434.945416930935, 
3441.392722793476, 3447.841612187623, 3454.292082609767, 3460.744131564204, 3467.197756563097, 3473.652955126438, 
3480.10972478201, 3486.568063065354, 3493.027967519732, 3499.489435696086, 3505.952465153006, 3512.417053456696, 
3518.883198180934, 3525.350896907038, 3531.820147223835, 3538.290946727617, 3544.763293022118, 3551.23718371847, 
3557.712616435174, 3564.189588798063, 3570.668098440272, 3577.148143002199, 3583.629720131476, 3590.112827482933, 
3596.597462718568, 3603.083623507512, 3609.571307525997, 3616.060512457322, 3622.551235991825, 3629.043475826846, 
3635.537229666697, 3642.032495222635, 3648.52927021282, 3655.027552362297, 3661.527339402952, 3668.028629073493, 
3674.531419119408, 3681.035707292945, 3687.541491353073, 3694.048769065458, 3700.55753820243, 3707.067796542953, 
3713.579541872598, 3720.09277198351, 3726.607484674383, 3733.123677750425, 3739.641349023338, 3746.160496311278, 
3752.681117438837, 3759.203210237006, 3765.726772543156, 3772.251802201, 3778.77829706057, 3785.306254978193, 
3791.835673816455, 3798.366551444181, 3804.898885736403, 3811.432674574337, 3817.967915845351, 3824.504607442942, 
3831.042747266709, 3837.582333222327, 3844.123363221518, 3850.665835182025, 3857.209747027589, 3863.755096687923, 
3870.301882098684, 3876.850101201446, 3883.39975194368, 3889.950832278723, 3896.503340165758, 3903.057273569784, 
3909.612630461595, 3916.169408817753, 3922.727606620566, 3929.287221858059, 3935.848252523955, 3942.410696617649, 
3948.974552144181, 3955.539817114216, 3962.106489544019, 3968.674567455431, 3975.244048875846, 3981.814931838185, 
3988.38721438088, 3994.96089454784, 4001.53597038844, 4008.112439957487, 4014.690301315209, 4021.269552527219, 
4027.850191664503, 4034.432216803397, 4041.015626025556, 4047.600417417941, 4054.186589072796, 4060.774139087621, 
4067.363065565155, 4073.95336661335, 4080.54504034536, 4087.138084879502, 4093.732498339251, 4100.328278853213, 
4106.925424555099, 4113.523933583715, 4120.123804082927, 4126.725034201657, 4133.327622093846, 4139.931565918447, 
4146.536863839394, 4153.143514025592, 4159.751514650888, 4166.360863894055, 4172.971559938773, 4179.583600973606, 
4186.196985191986, 4192.81171079219, 4199.427775977324, 4206.045178955298, 4212.663917938814, 4219.283991145344, 
4225.905396797108, 4232.528133121057, 4239.152198348858, 4245.777590716866, 4252.404308466115, 4259.032349842295, 
4265.661713095731, 4272.292396481374, 4278.92439825877, 4285.557716692051, 4292.192350049911, 4298.828296605599, 
4305.465554636883, 4312.104122426051, 4318.743998259877, 4325.385180429617, 4332.027667230984, 4338.671456964132, 
4345.316547933638, 4351.962938448485, 4358.610626822049, 4365.259611372073, 4371.909890420661, 4378.561462294251, 
4385.214325323604, 4391.868477843788, 4398.523918194155, 4405.180644718333, 4411.838655764204, 4418.497949683887, 
4425.158524833728, 4431.820379574272, 4438.483512270263, 4445.147921290613, 4451.813605008396, 4458.480561800826, 
4465.148790049243, 4471.8182881391, 4478.489054459947, 4485.161087405408, 4491.834385373175, 4498.50894676499, 
4505.184769986625, 4511.861853447872, 4518.540195562527, 4525.219794748371, 4531.900649427162, 4538.582758024611, 
4545.266118970379, 4551.950730698046, 4558.636591645115, 4565.323700252981, 4572.012054966927, 4578.701654236107, 
4585.392496513525, 4592.084580256032, 4598.777903924302, 4605.472465982823, 4612.168264899882, 4618.865299147548, 
4625.563567201663, 4632.263067541824, 4638.963798651372, 4645.665759017375, 4652.368947130615, 4659.07336148558, 
4665.77900058044, 4672.485862917043, 4679.193947000896, 4685.903251341154, 4692.613774450608, 4699.325514845663, 
4706.03847104634, 4712.752641576249, 4719.468024962584, 4726.184619736105, 4732.902424431129, 4739.621437585514, 
4746.341657740649, 4753.06308344144, 4759.785713236295, 4766.509545677116, 4773.234579319283, 4779.960812721641, 
4786.688244446493, 4793.416873059578, 4800.146697130068, 4806.87771523055, 4813.609925937017, 4820.343327828854, 
4827.077919488827, 4833.81369950307, 4840.550666461073, 4847.288818955668, 4854.028155583026, 4860.768674942631, 
4867.510375637283, 4874.253256273076, 4880.997315459387, 4887.742551808871, 4894.488963937443, 4901.236550464274, 
4907.985310011765, 4914.735241205553, 4921.48634267449, 4928.238613050631, 4934.992050969229, 4941.746655068717, 
4948.502423990701, 4955.259356379949, 4962.017450884376, 4968.77670615504, 4975.537120846123, 4982.298693614928, 
4989.061423121859, 4995.825308030422, 5002.590347007202, 5009.356538721863, 5016.123881847128, 5022.892375058777, 
5029.662017035629, 5036.432806459539, 5043.204742015378, 5049.977822391034, 5056.752046277391, 5063.527412368328, 
5070.303919360701, 5077.081565954335, 5083.86035085202, 5090.640272759493, 5097.421330385428, 5104.203522441435, 
5110.98684764204, 5117.771304704676, 5124.556892349685, 5131.34360930029, 5138.131454282599, 5144.920426025591, 
5151.710523261106, 5158.501744723831, 5165.294089151302, 5172.087555283882, 5178.882141864757, 5185.677847639932, 
5192.474671358207, 5199.272611771182, 5206.07166763324, 5212.871837701543, 5219.673120736014, 5226.475515499339, 
5233.279020756947, 5240.083635277009, 5246.889357830426, 5253.696187190818, 5260.504122134518, 5267.313161440561, 
5274.123303890677, 5280.934548269279, 5287.746893363456, 5294.560337962967, 5301.374880860228, 5308.190520850301, 
5315.007256730897, 5321.825087302352, 5328.644011367627, 5335.464027732301, 5342.285135204558, 5349.107332595179, 
5355.930618717534, 5362.754992387577, 5369.580452423832, 5376.406997647388, 5383.234626881892, 5390.063338953533, 
5396.893132691046, 5403.724006925692, 5410.555960491258, 5417.388992224043, 5424.223100962857, 5431.058285549004, 
5437.894544826282, 5444.731877640967, 5451.570282841815, 5458.409759280044, 5465.250305809333, 5472.091921285811, 
5478.934604568049, 5485.778354517056, 5492.623169996265, 5499.469049871529, 5506.315993011114, 5513.16399828569, 
5520.013064568323, 5526.863190734469, 5533.714375661963, 5540.566618231015, 5547.419917324201, 5554.274271826455, 
5561.129680625066, 5567.986142609661, 5574.843656672206, 5581.702221706998, 5588.561836610652, 5595.422500282099, 
5602.28421162258, 5609.146969535633, 5616.010772927086, 5622.875620705055, 5629.74151177994, 5636.608445064402, 
5643.476419473372, 5650.345433924037, 5657.215487335836, 5664.086578630447, 5670.958706731785, 5677.831870565997, 
5684.706069061451, 5691.581301148727, 5698.457565760618, 5705.334861832115, 5712.213188300408, 5719.092544104868, 
5725.972928187054, 5732.854339490696, 5739.736776961693, 5746.620239548106, 5753.50472620015, 5760.390235870183, 
5767.276767512714, 5774.16432008438, 5781.052892543946, 5787.9424838523, 5794.833092972446, 5801.724718869499, 
5808.617360510671, 5815.511016865274, 5822.405686904708, 5829.301369602456, 5836.198063934078, 5843.095768877208, 
5849.994483411537, 5856.894206518822, 5863.794937182867, 5870.696674389523, 5877.599417126682, 5884.503164384265, 
5891.407915154228, 5898.313668430539, 5905.220423209188, 5912.128178488171}; 

///////////////////////////////////////////////////////////////////////////

double lnf(double n)
  {
  int i = n<0 ? 0 : int(n+0.1);
  return lnfactorial[i];
  }

