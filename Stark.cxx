/**
 * @file Stark.cxx
 *
 * @brief Stark shifts of hydrogenic energy levels (due to core polarization)
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: Stark.cxx 375 2016-02-05 16:40:12Z iamp $
   @endverbatim
 *
 */
#include <stdio.h>
#include <math.h>
#include "clebsch.h"

using namespace std;

/**
* @brief average 1st order Stark shift of a hydrogenic n,l level in Rydbergs 
*
* to be multiplied by the electric field strength in atomic units 
*/
double nlStark(double z, int n, int l) 
   {
   double Stark = 0.0;
   for (int m = 0; m<=n; m++)
     {
     if (m>l) break;
     double mult = (m==0) ? 1 : 2;
     int n1max = n-1-m;
     for (int n1=0; n1<=n1max; n1++)
       {
       int n2 = n1max-n1;
       double q = fabs(double(n1-n2));
       double cg = ThreeJ(0.5*(n-1),0.5*(n-1),l,0.5*(m+q),0.5*(m-q),-m);
       Stark += mult*cg*cg*3.0*n*q/z;
       }
     }
   return Stark;
   }
/////////////////////////////////////////////////////////////////////////////
/**
 * @brief  Binding energy of a Rydberg electron due to core polarization in Rydbergs
 */
double nlPol(double z, int n, int l)
   {
   double A, term1, term2;
   A = 1160.0/exp(3.52*log(2.21+z)); // polarizability for Li-like ions
   term1 = pow(z,4)*(3*n*n-l*(l+1));
   term2 = 2.0*pow(double(n),5)*(l+1.5)*(l+1.0)*(l+0.5)*l*(l-0.5);
   return A*term1/term2;
   }

void CalcminField(void)
   {
   double z;
   int n;

   char filename[200];
   FILE* fout;
   printf("\n Give ion charge .........: ");
   scanf("%lf",&z);
   printf("\n Give main quantum number : ");
   scanf("%d",&n);
   printf("\n Give filename for output : ");
   scanf("%s",&filename);
   fout = fopen(filename,"w");

   for (int l=1;l<n;l++)
     {
     double nls = fabs(nlStark(z,n,l)*13.6056);
     double nlp = nlPol(z,n,l)*13.6056;
     printf("%5d %12.6g %12.6g %12.6g\n",l,nls,nlp,nlp/nls*5.142e9);
     fprintf(fout,"%5d %12.6g %12.6g %12.6g\n",l,nls,nlp,nlp/nls*5.142e9);
     }
   fclose(fout);
   }

