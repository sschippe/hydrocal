/**
 * @file RRphotons.cxx
 *
 * @brief Hydrogenic calculation of RR photon spectra including cascades (unfinished)
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: RRphotons.cxx 375 2016-02-05 16:40:12Z iamp $
   @endverbatim
 *
 */
#include <stdio.h>
#include <math.h>
#include <cstring>
#include <iostream>
#include "radrate.h"
#include "sigmarr.h"
#include "hydromath.h"
#include "RRphotons.h"

using namespace std;

//#define printcasc 1

//////////////////////////////////////////////////////////////////////////

inline int index2(int n, int l) { return n*(n-1)/2+l;}
// index of two dimensional n,l array

inline int index4(int n1,int l1, int n2, int l2, int nmax) { return nmax*(nmax+1)/2+n1*(n1-1)/2+l1+n2*(n2-1)/2+l2; } 
// returns index of four dimensional n1,l1,n2,l2  array

//////////////////////////////////////////////////////////////////////////

double emission_profile(double t, double tL, double tau)
{
    double term1 =  t < tL ? 1 : exp(-(t-tL)/tau);
    double term2 = exp(-t/tau);
    return (term1-term2)/tL;
}
//////////////////////////////////////////////////////////////////////////

double emission_cascade(int counter, int &ncascstep, int n1, int l1, int &nmax, double &t, 
                        double &Tcool, int *n_list, int *l_list, double *pr_list, 
                        double *lt_list, const RADRATE& hydro, double *sigma, double *profile)
{ // this function calles itself recursively,
  // the recusion is terminated when either counter>=ncascstep or when n2>nmax,
  // whenever possible variables are passed by reference in order to save memory
  counter++;

  n_list[counter] = n1;
  l_list[counter] = l1;
  lt_list[counter] = hydro.life(n1,l1);
  pr_list[counter] = profile[index2(n1,l1)];

  static int i,k;
  static double prod;
  double emisXsec = 0.0;
  for (i=0; i<=counter; i++)
  {
    prod = 1.0;
    for (k=0; k<=counter; k++) 
    {
	if (k!=i) { prod *= 1.0-lt_list[k]/lt_list[i]; } 
    } 
    emisXsec += pr_list[i]/prod;
  }
  emisXsec *= sigma[index2(n1,l1)];

  if (counter<ncascstep) 
  {
      for (int n2=n1+1;n2<=nmax;n2++)
      {   
	  if (l1>=1)
	  {
         
	  emisXsec += hydro.branch(n2,l1-1,n1,l1)*
	      emission_cascade(counter, ncascstep, n2, l1-1, nmax, t, Tcool,
                               n_list, l_list, pr_list, lt_list, hydro, sigma, profile);
	  }  
	  emisXsec += hydro.branch(n2,l1+1,n1,l1)*
	      emission_cascade(counter, ncascstep, n2, l1+1, nmax, t, Tcool,
                               n_list, l_list, pr_list, lt_list, hydro, sigma, profile);  
      } 
 } 

#ifdef printcasc
  // diagnostic output
  for(i=counter;i>=0;i--)
    {
    printf("%3d %3d -> ",n_list[i],l_list[i]);
    }
  printf(" %12.4g\n",emisXsec);
#endif

  return emisXsec;
  }


//////////////////////////////////////////////////////////////////////////

void calc_emission()
  {
  const double clight=2.99792458e10;
  const double melectron=510999.6;
  const double pi = 3.1415926536;
  char answer[2], filenameroot[200], filename[200], pfn[200];
  FILE *fout;

  double z;
  printf("\n Give effective nuclear charge ...................: ");
  scanf("%lf",&z);

  int nmax;
  printf("\n Give maximum principal quantum number ...........: ");
  scanf("%d",&nmax);

  if (nmax>nmaxfactorial/2.0) // nmaxfactorial is defined in factorial.h 
       {
       nmax = int(nmaxfactorial/2.0);
       printf("\n nmax exceeds nmaxfactoiral/2");
       printf(" which is defined in factorial.h");
       printf("\n The calculation of all necessary Clebsch-Gordan");
       printf("\n coefficients is not prossible.");
       printf("\n nmax set to nmaxfactorial/2 = %4d.\n",nmax);
       }

  printf("\n Now calculating hydrogenic transition rates and");
  printf(" decay probabilities  \n");

  RADRATE hydro(nmax,z,1.0,2.0,2.0); // last three arguments are dummy values

  int nmaxdim = nmax*(nmax+1)/2;

  int ni,li,nf,lf;
  printf("\n\n Give principal quantum numbers n -> n' ..........: ");
  scanf("%d %d",&ni,&nf);
  
  int ncascstep;
  printf("\n Give number of cascade steps (<0 =nmax) .........: ");
  scanf("%d",&ncascstep);
  if (ncascstep<0) ncascstep=nmax;

  // calculation of n,l-specific RR cross sections
  double electron_energy;
  printf("\n Give electron-ion collision energy (eV) .........: ");
  scanf("%lf",&electron_energy);

  printf("\n\n Now calculating RR cross sections at E = %12.4g eV\n",electron_energy);
  double *sigmaRR = new double[nmaxdim];
  int n,l;
  for (n=1; n<=nmax; n++)
  {
      if ((n % 10)==0) { printf("|"); } else { printf("."); }
     for (l=0; l<n; l++) sigmaRR[index2(n,l)] = sigmarrqm(electron_energy,z,n,l);
  }


  double mass,ecool,Lcool;
  printf("\n\n Give atomic mass (amu) ..........................: ");
  scanf("%lf",&mass);
  printf("\n Give cooling energy (eV) ........................: ");
  scanf("%lf",&ecool);
  printf("\n Give cooler length (cm) .........................: ");
  scanf("%lf",&Lcool);

  double gamma = 1.0+ecool/melectron;
  double vel = clight*sqrt(1.0-1.0/(gamma*gamma));
  double Tcool = Lcool/vel; 
  printf("\n The ion velocity is .............................: %12.4g cm/s",vel);
  printf("\n The flight time through the cooler is ...........: %12.4g s\n",Tcool);

  double dmin,dmax,dstep;
  printf("\n Give start of observed lab. flight distance (cm).: ");
  scanf("%lf",&dmin);
  printf("\n Give end of observed lab. flight distance (cm) ..: ");
  scanf("%lf",&dmax);
  double tmin = dmin/vel/gamma;  
  double tmax = dmax/vel/gamma;  
  printf("\n Give step width (cm) ............................: ");
  scanf("%lf",&dstep);
  double dt = dstep/vel/gamma;

  printf("\n Give filename for output (*.dat) ................: ");
  scanf("%s",&filenameroot);
  strcpy(filename,filenameroot);
  strcat(filename,".dat");
  fout = fopen(filename,"w");
  fprintf(fout,"#      %3d -> %3d emission for Zeff =  %8.4f\n",ni,nf,z);
  fprintf(fout,"#  maximum principal quantum number : %3d\n",nmax);  
  fprintf(fout,"#           number of cascade steps : %3d\n",ncascstep);
  fprintf(fout,"#     electron-ion collision energy : %12.4g eV\n",electron_energy);
  fprintf(fout,"#                          ion mass : %12.4g u\n",mass);
  fprintf(fout,"#                      ion velocity : %12.4g cm/s\n",vel);
  fprintf(fout,"#                    cooling energy : %12.4g eV\n",ecool);
  fprintf(fout,"#                     cooler length : %12.4g cm\n",Lcool);
  fprintf(fout,"#------------------------------------------------------\n");
  fprintf(fout,"proper time(s) lab.dist.(cm)  em.Xsec(cm2)\n");   

  int *n_list = new int[nmax+1];
  int *l_list = new int[nmax+1];
  double *pr_list = new double[nmax+1];
  double *lt_list = new double[nmax+1];
  double *profile = new double[nmaxdim];
  if (tmin==0.0)
    {
      fprintf(fout," %13.4g %13.4g %13.4g\n",0.0,0.0,0.0);
      tmin = dt;
    }
  tmax += 0.5*dt;   
  for (double t=dt; t<=tmax; t+=dt)
  {
      printf("\n Emission cross section at z = %12.4g cm:",t*vel*gamma);
      for(n=1;n<=nmax;n++)
	{
	  for (l=0;l<n;l++)
	    {
	      profile[index2(n,l)] = emission_profile(t,Tcool,hydro.life(n,l));
    	    }
	}
      double emisXsec = 0.0;
      for (li = 0; li <= nf; li++)
      { 
	  if ((nf==1)&&(li==0)) continue;
	  double Xsec = sigmaRR[index2(ni,li)]*profile[index2(ni,li)]; 
          if (ncascstep>0)  
	  {
	      for (int n2=ni+1;n2<=nmax;n2++)
	      {
		  int counter = 0;
		  n_list[0]  = ni; 
		  l_list[0]  = li;
                  pr_list[0] = profile[index2(ni,li)]; 
		  lt_list[0] = hydro.life(ni,li);;
		  if (li>=1)
		  {
		      Xsec += hydro.branch(n2,li-1,ni,li)*
			  emission_cascade(counter, ncascstep, n2, li-1, nmax, t, Tcool,
					   n_list, l_list, pr_list, lt_list, hydro, sigmaRR, profile);
		  }  
		  Xsec += hydro.branch(n2,li+1,ni,li)*
		      emission_cascade(counter, ncascstep, n2, li+1, nmax, t, Tcool,
				       n_list, l_list, pr_list, lt_list, hydro, sigmaRR, profile);  
	      }
	  }
	  double br = (li>0) ?  hydro.branch(ni,li,nf,li-1) : 0.0;
	  if (li<nf-1) br += hydro.branch(ni,li,nf,li+1);
          emisXsec += Xsec*br;
      } 
      printf(" %12.4g cm2",emisXsec);
      fprintf(fout," %13.4g %13.4g %13.4g\n",t,t*vel*gamma,emisXsec);
  }
  fclose(fout);  
 
  delete[] profile;
  delete[] lt_list;
  delete[] pr_list;
  delete[] l_list;
  delete[] n_list;
  delete[] sigmaRR;;
  
}


/* Achtung Baustelle!
for (int m=0; m<=n; m++) // m is the number of cascade steps from state n,l down to state ni,li
  {
    int imax = over(n,m);    // number of different paths n -> n1 -> n2 -> n3-> ... -> ni
    int jmax = power(2,m);   // number of different paths l -> l1 -> l2 -> l3-> ... -> li
    for (int i=0; i<=imax; i++) 
      {
	int n2 = n, l2 = l, n1, l1;
	double tau2=hydro.life(n2,l2), tau1, sum = 0.0, prod = 1.0;
	for (int j=0; j<=jmax; j++) 
	  {
	    int dl = 2*(j & power(2,i))-1;
	    n1 = n2; l1 = l2;
	    l2 = l1+dl;
	    if (l2<0) break;
	    n2 = f(n,ni,m,i); // the function f still needs to be discovered
	    if (l2>=n2) break;
	    tau1 = tau2;
	    tau2 = hydro.life(n2,l2);
            sum += profile[index(n2,l2)]/(1.0-tau1/tau2);
            prod *= hydro.branch(n1,l1,n2,l2);
	  } // end for j
      } // end for i
  } // end for m
*/

//////////////////////////////////////////////////////////////////////////

// calculation of cascade matrix
/*
  double *branch = new double[2*nmaxdim];
  double *taufac = new double[2*nmaxdim];
  double *cascade = new double[2*nmaxdim];  

  for (n1=1; n1<=nmax; n1++)
  {
      for (l1=0;l1<n1;l1++)
      {
	  for(l2=0;l2<n1;l2++)
	  {
	      cascade[index4(n1,l1,n1,l2,nmax)] = (l1==l2) ? 1.0 : 0.0;
	  }
	  for (n2=1;n2<nmax;n2++)
	  {
	      if (n2<n1)
	      {
		  branch[index4(n1,l1,n2,l1+1,nmax)] = hydro.life(n1,l1)*hydro.trans(n1,l1,n2,l1+1);
		  branch[index4(n1,l1,n2,l1-1,nmax)] = hydro.life(n1,l1)*hydro.trans(n1,l1,n2,l1-1);
	      }
	      for (l2=0;l2<n2;l2++)
	      {
		  taufac[index4(n1,l1,n2,l2,nmax)] = (n1==n2) ? 0.0 : 1.0/(1.0-hydro.life(n2,l2)/hydro.life(n1,l1));
	      }
	  }
          for (n2=n1-1;n2>=1;n2--)
	  {
	      for (l2=0;l2<n2;l2++)
	      {
		  for (n=n2+1;n<=n1;n++)
		  {
		      double term1 = cascade[index4(n1,l1,n,l2+1,nmax)]*branch[index4(n,l2+1,n2,l2,nmax)];
		      double term2 = cascade[index4(n1,l1,n,l2-1,nmax)]*branch[index4(n,l2-1,n2,l2,nmax)];
		      cascade[index4(n1,l1,n2,l2,nmax)] += term1+term2;
		  }
	      }
	  }
      }
  }
*/
