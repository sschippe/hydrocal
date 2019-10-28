/** 
    @file lifetime.cxx
    @brief Calculation of hydrogenic radiative dipole transition rates and related quantities 
 
    @par CREATION  
    @author Stefan Schippers
    @date 1997, 1999
    
    @par VERSION
    @verbatim
    $Id: lifetime.cxx 446 2017-08-28 16:08:49Z iamp $
    @endverbatim   

 */ 
#include <stdio.h>
#include <math.h>
#include <cstring>
#include "osci.h"

using namespace std;

double hydrotrans(double z, int n1, int l1, int n2, int l2)
   {
   const double tpunit = 8.032e9;
   double etrans = z*z*(1.0/(n2*n2)-1.0/(n1*n1));
   return etrans*etrans*tpunit*(2.0*l2+1.0)/(2.0*l1+1.0)*fosciBB(n2,l2,n1,l1);
   }

double hydrolife(double z, int n1, int l1, int n2min=0)
   {
   double sum = 0.0;

   for (int n2=n1-1;n2>=l1;n2--)
     {
     if (n2<n2min) break;
     if (l1+1<n2) sum += hydrotrans(z,n1,l1,n2,l1+1);
     if ((l1-1<n2)&&(l1>0)) sum += hydrotrans(z,n1,l1,n2,l1-1);
     }
   return sum > 0.0 ? 1.0/sum : 1;
   }

double hydrobranch(double z, int n1, int l1, int n2, int l2)
   {
   return hydrolife(z,n1,l1)*hydrotrans(z,n1,l1,n2,l2);
   }


void fixed_n1_n2(void)
  {
   float z, tr, br, lt; 
   int   n1, l1, n2, l2, i, l1max;
   printf("\n Give Z, n and n' for n -> n transition : ");
   scanf("%f %d %d",&z,&n1,&n2);
   printf("\n                   ");
   printf("       life-");
   printf("   transition");
   printf("    branching\n");
   printf(" n, l -> n',l'  :  ");
   printf("    time [s]");
   printf("   rate [1/s]");
   printf("        ratio\n");
   l1max = (n1<=n2) ? n1 : n2+1;
   for (l1=0;l1<l1max;l1++)
     {
     lt = hydrolife(z,n1,l1);
     for (i=0;i<2;i++)
       {
       l2=l1-(1-2*i);
       if ((l2>=0) & (l2<n2))
         {
         tr = hydrotrans(z,n1,l1,n2,l2);
         br = lt*tr; 
         printf("%2d,%2d -> %2d,%2d  :  %12.4g %12.4g %12.4g\n",
                n1,l1,n2,l2,lt,tr,br);
         }
       }
     }
  printf("\n");
  }

void fixed_n2_l2_l1(void)
  {
    double z, tr; 
    int   n1, l1, l2, n2, nmax;
    printf("\n Give Z, n', l', nmax, l: ");
    do scanf("%lf %d %d %d %d",&z,&n2,&l2,&nmax,&l1); 
    while(((l1-l2)!=1)&&((l2-l1)!=1));
    printf("\n  n,  l ->  n', l' :   rate [1/s], rate*n^3 [1/s]\n");
    for (n1=n2+1; n1<=nmax; n1++)
      {
	tr = hydrotrans(z,n1,l1,n2,l2);
	printf("%3d,%3d -> %3d,%3d : %12.4g, %12.4g\n",n1,l1,n2,l2,tr,tr*n1*n1*n1);
      }
  }

void all_lifetimes(void)
  {
  double z;
  int nmax, choice;
  char filename[200];

  printf("\n Give nuclear charge Z and maximum n : ");
  scanf("%lf %d",&z,&nmax);
  printf("\n Which kind of output?");
  printf("\n 1: lifetimes in s");
  printf("\n 2: lifetimes relative to the maximum value per n");
  printf("\n 3: transition rates in 1/s");
  printf("\n Make a choice ......................: ");
  scanf("%d",&choice);
  printf("\n Give filename for output (*.tau) ...: ");
  scanf("%s",&filename);

  FILE *fout;
  strcat(filename,".tau");
  fout = fopen(filename,"w");

  double *tau = new double[nmax];
  
  for(int n=1; n<=nmax; n++)
    {
    int l;
    double taumax = 0.0;
    for(l=0; l<n; l++)
      {
      tau[l] = hydrolife(z,n,l);
      if (tau[l]>taumax) taumax=tau[l];
      }
    switch (choice)
      {
        case 2: for(l=0; l<n; l++) {fprintf(fout," %12.6g",tau[l]/taumax);}
                break;
        case 3: for(l=0; l<n; l++) 
                  {
                  double rate = tau[l] > 0 ? 1.0/tau[l] : 999;
                  fprintf(fout," %12.6g",rate);
                  }
                break;
       default: for(l=0; l<n; l++) {fprintf(fout," %12.6g",tau[l]);}
                break;
       }
    for(int k=n;k<nmax;k++) fprintf(fout,"%12.6g",0.0);
    fprintf(fout,"\n");
    }
  delete[] tau;
  fclose(fout);
  }

void lifetime_nl(void)
{
    double z;
    int n,l,nmin=0;
    printf("\n Give nuclear charge Z and n, l : ");
    scanf("%lf %d %d",&z,&n,&l); 
    printf("\n Give min n of final state      : ");
    scanf("%d",&nmin);
    double  tau = hydrolife(z,n,l,nmin);
    printf("\n        Lifetime: %14.8g s",tau);
    printf("\n Transition rate: %14.8g /s",1.0/tau);
}

void testLifetime(void)
  {
   int choice;
   printf("\n 1: Fixed n', l' and l");
   printf("\n 2: all n -> n' transitions");
   printf("\n 3: all hydrogenic lifetimes up to nmax");
   printf("\n 4: lifetime of a specific n,l level");
   printf("\n\n Make a choice ...........: ");
   scanf("%d",&choice);
   
   switch (choice) 
     {
     case 2 : fixed_n1_n2();break;
     case 3 : all_lifetimes(); break;
     case 4 : lifetime_nl(); break;
     default: fixed_n2_l2_l1(); break;
     }
  }



