/**
 * @file rydserie.cxx
 *
 * @brief Hydrogenic calculation of enegies of Rydberg levels
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: rydserie.cxx 651 2019-10-15 17:45:59Z iamp $
   @endverbatim
 *
 */
#include "hydroconst.h"
#include <cstdio>
#include <cstring>
#include <cmath>

using namespace std;

void rydbergseries(void)
  {
    const double alpha=_hydroconst_alpha;
    const double mec2=_hydroconst_mec2_eV;
    const double Ryd = _hydroconst_Ryd_eV;
    
    double limit, Zeff, energy;
    int nmax;
    char answer[2];
    
    printf("\n Give effective nuclear charge .......: ");
    scanf("%lf",&Zeff);
    printf("\n Give series limit (eV) ..............: ");
    scanf("%lf",&limit);
    printf("\n Give maximum n ......................: ");
    scanf("%d",&nmax);
    printf("\n Bohr or Dirac energies? (B/D) .......: ");
    scanf("%s",&answer); 
    bool dirac_flag = ( (!strcmp(answer,"d")) || (!strcmp(answer,"D")) );

    if (dirac_flag)
      {
	printf("\n\n    n      j    E(eV)\n\n");
      }
    else
      {
	printf("\n\n    n      E(eV)\n\n");
      }

    for (int n=1; n<=nmax; n++)
    {
      double dn = (double)n;
      if (dirac_flag)
	{ // Dirac
	  for (int k=1; k<=n; k++)
	    {
	      double dk = (double)k;
	      double za = Zeff*alpha;
	      double term = za/(dn-dk+sqrt(dk*dk-za*za));
	      energy = limit + mec2*(1.0/sqrt(1.0+term*term)-1.0);
	      if (energy>0) printf(" %4d %4d/2 %12.8g\n",n,2*k-1,energy);
	    } // end for k
	} // end Dirac
      else
	{ // Bohr 
	  double zn = Zeff/dn;
	  energy = limit-Ryd*zn*zn;
	  if (energy>0) printf(" %4d   %12.8g\n",n,energy);
	} // end Bohr
    } // end for n
  }
