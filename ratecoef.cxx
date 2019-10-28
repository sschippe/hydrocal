/**
 * @file ratecoef.cxx
 *
 * @brief convolution of theoretical recombination cross sections
 *                                                                         
 * Convolution of recombination cross sections with an electron energy     
 * distribution defined in FELE.CXX. RR cross sections are defined in      
 * SIGMARR.CXX. For compatibilty reasons DR cross sections, which are      
 * defined here, have the same list of parameters as RR cross sections.    
 * Note that all functions defining a cross section return SIGMA(E)*E.     
 * 
 * @verbatim
 $Id: ratecoef.cxx 596 2019-03-26 19:52:12Z iamp $
 @endverbatim
 *
 * Stefan Schippers@physik.uni-giessen.de
*/                                                                         

#include <cstdio>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include "convolute.h"
#include "fraction.h"
#include "sigmarr.h"
#include "gaussint.h"
#include "readxsec.h"
#include "fele.h"
#include "svnrevision.h"

using namespace std;

double fsigma(SIGMARR sigma, FELE fele, double e, double eV, 
              double ktpar, double ktperp, 
              double z, double* fraction, int use_fraction,
              int nmin, int nmax, int lmin, int *counter)
  // integrand of the convolution integral, i.e. f(v)*v*sigma(v)
   {
   (*counter)++;
   double vsigma = (*sigma)(e,z,fraction,use_fraction,nmin,nmax,lmin);
          vsigma *= Clight*sqrt(2.0/Melectron/e);

   return (*fele)(eV,e,ktpar,ktperp)*vsigma;
   }

//////////////////////////////////////////////////////////////////////////

int comp(const double *i1, const double *i2){return ((*i1)<(*i2)) ? -1 : 1;}
// comparison used in qsort

 
//////////////////////////////////////////////////////////////////////////
/**
 * @brief  integration routine performing Gauss-Legendre and Gauss-Laguerre integrations
 *
 * Gauss-Legendre and Gauss-Laguerre integrations integrations are for finite and semi-infinite subintervals, respectivly
 *
 * @param sigma cross section function of type  @ref SIGMARR
 * @param fele electron energy distribution function of type @ref FELE
 * @param theo_mode not used
 * @param eV electron energy (after convolution)
 * @param ktpar parallel temperature of electron beam
 * @param ktperp transversal temperature of electron beam
 * @param peaks[0] widths of electron energy distribution
 * @prama peaks[1] position of DR peak (if present)
 * @param peaks[2] width of DR peak
 * @param fraction array contining field-ionization survival fractions
 * @param use_fraction (0: fractions not used, 1: fraction used)
 * @param nmin minimum principal quantum number
 * @param nmax maximum principal quantum number
 * @param lmin minimum angular momentum quantum number in shell n=nim
 * @param subintervals number of integration intervals
 * @param delta integration step size
 * @param legendre_x abscissae of Legendre integration
 * @param legendre_w weights of Legendre integration
 * @param legendre_n number of points for Legendre integration
 * @param laguerre_x abscissae of Gauss-Laguerre integration
 * @param laguerre_w weights of Gauss-Laguerre integration
 * @param laguerrre_n number of points for Gauss-Laguerre integration
 */
double convolute(SIGMARR sigma, FELE fele, int theo_mode, double eV, 
       double ktpar, double ktperp, 
       double *peaks, double z, double* fraction, int use_fraction,
       int nmin, int nmax, int lmin, int subintervals, double delta,
       double* legendre_x, double* legendre_w, int legendre_n, 
       double* laguerre_x, double* laguerre_w, int laguerre_n)
   {
   int counter=0;

   int i;
   int intervals=4*(subintervals+1);
   double *e = new double[intervals];
// use integration intervals around peak of electron energy distribution
   e[0] = 1.0e-20;
   for (i=1;i<=subintervals; i++)
     { 
     e[i] = eV-delta*(subintervals-i+1)*peaks[0];
     }
   e[subintervals+1] = eV;
   for (i=subintervals+2; i<=2*subintervals+1; i++)
     {
     e[i] = eV+delta*(i-subintervals-1)*peaks[0];
     }  
   intervals = 2*subintervals+1;
// use additional integration intervals around DR peak (if present)
   if (peaks[1]>0) 
     { 
     for (i=intervals+1; i<=intervals+subintervals; i++)
       { 
       e[i] = peaks[1]-delta*(subintervals-i+1+intervals)*peaks[2];
       }
     e[intervals+subintervals+1] = peaks[1];
     for (i=intervals+subintervals+2; i<=intervals+2*subintervals+1; i++)
       {
       e[i] = peaks[1]+delta*(i-subintervals-1-intervals)*peaks[2];
       }  
     intervals += 2*subintervals+1;
     }

// sort integration intervals in ascending order
   qsort((void*)e,intervals+1,sizeof(double),
         (int(*)(const void*,const void*))(comp));
   
   double ee, sum, result = 0.0; 

   for (int k=1; k<=intervals; k++)
     {
     if (e[k-1]<=0) continue;
     sum = 0.0;
     for(i=legendre_n; i>0; i--)
        {
        ee = 0.5*((e[k]-e[k-1])*legendre_x[i]+e[k]+e[k-1]);
        double y = legendre_w[i]*
                 fsigma(sigma,fele,ee,eV,ktpar,ktperp,z,fraction,
                        use_fraction,nmin,nmax,lmin,&counter);
	if (theo_mode==10) { 
	  if (ee==0.0) continue;
	  y /= ee;
	}
	sum += y;
        }
     result += sum*0.5*(e[k]-e[k-1]);
     }

   sum = 0.0;
   for (i=1; i<=laguerre_n; i++)
     {
     ee = laguerre_x[i]+e[intervals];
     double y = laguerre_w[i]*exp(laguerre_x[i])*
              fsigma(sigma,fele,ee,eV,ktpar,ktperp,z,fraction,
                     use_fraction,nmin,nmax,lmin,&counter);
     sum += y;
     }
   result += sum;

   delete[] e;
   return result;
   }

//////////////////////////////////////////////////////////////////////////
/**
 * @brief interactive input of parameters for numerical integration
 *
 * @param ktpar parallel temperature of electron beam
 * @param ktperp transversal temperature of electron beam
 * @param nsub on exit, number of integration intervals
 * @param delta on exit, integration step size
 * @param on exit, legendre_x abscissae of Legendre integration
 * @param on exit, legendre_w weights of Legendre integration
 * @param on exit, legendre_n number of points for Legendre integration
 * @param on exit, laguerre_x abscissae of Gauss-Laguerre integration
 * @param on exit, laguerre_w weights of Gauss-Laguerre integration
 * @param on exit, laguerrre_n number of points for Gauss-Laguerre integration
 */
void setup_integration(double ktpar, double ktperp, int *nsub, double *delta,
            int *legendre_n, double *legendre_x, double *legendre_w,
            int *laguerre_n, double *laguerre_x, double *laguerre_w)
   {
   char answer[2];
   int flag = 1;

   printf("\n Use default integration parameters (y/n) ..............: ");
   scanf("%s",&answer);
   if (strcmp(answer,"n")) 
      {
      *legendre_n = 12;
      *laguerre_n = 4;
      *nsub = 10;
      *delta = 1.0;
      gauss_cof(*legendre_n,legendre_x,legendre_w,&flag);
      laguer_cof(*laguerre_n,laguerre_x,laguerre_w,&flag);
      return;
      }

   SIGMARR stest = &sigmarrtest;
   FELE    fele  = &fecool;
   double eV,z=1.0,t=1.0,alpha,norm;
   int selection, nmin=1, nmax=1, lmin=0, nele=0;


   printf("\n Test of the integration routine:  \n\n");
   printf("     The integration routine will integrate\n");
   printf("     the normalized electron energy distribution.\n"); 
   printf("     A correct integration should yield: norm = 1.\n");
   printf("     The number of integration points has to be chosen\n");
   printf("     accordingly. A reasonable guess is 12 Gauss-Legendre\n");
   printf("     and 4 Gauss-Laguerre integration points.\n\n");
   printf("     The integration range is subdivided into 'nsub'\n");
   printf("     subintervals on either side of the electron energy\n");
   printf("     distribution's maximum. The widths of the subintervals\n");
   printf("     is set 'delta'*FWHM, with the parameter 'delta' to be\n");
   printf("     supplied by the user. Reasonable parameter values are\n");
   printf("     nsub = 10 and delta=1.0.\n\n"); 
   for(;;) // test loop
     {
     flag = 1;
     while(flag)
        {
        printf(" Give number of Gauss-Legendre integration points\n");
        printf("    2..15, 16, 20, 24, 32, 40, 48, 64, 80, 96 ..........: ");
        scanf("%d",legendre_n);
        gauss_cof(*legendre_n,legendre_x,legendre_w,&flag);
        }

     flag = 1;
     while(flag)
        {
        printf(" Give number of Gauss-Laguerre integration points 2..15 : ");
        scanf("%d",laguerre_n);
        laguer_cof(*laguerre_n,laguerre_x,laguerre_w,&flag);
        }

     printf("\n Give number of integration subintervals (nsub) ........: ");
     scanf("%d",nsub);
     printf("\n Give width of integration subintervals (delta) ........: ");
     scanf("%lf",delta);

     selection = 0;
     while (selection==0)
        {
        printf("\n Give electron energy (eV) .............................: ");
        scanf("%lf",&eV);
        double normpeak[3];
        normpeak[0] = log(2.0)*ktperp+4*sqrt(log(2.0)*eV*ktpar);
        normpeak[1] = 0; normpeak[2] = 0;
        norm = convolute(stest,fele,0,eV,ktpar,ktperp,normpeak,z,&t,0,
                         nmin,nmax,lmin,
                         *nsub,*delta,legendre_x,legendre_w,*legendre_n,
                         laguerre_x,laguerre_w,*laguerre_n);
        printf("\n norm = %12.4G\n",norm);
        printf("\n What to do next?");
        printf("\n                            another energy : 0");
        printf("\n              other integration parameters : 1");
        printf("\n          calculation of rate coefficients : 2");
        printf("\n                               make a choice : ");
        scanf("%d",&selection);
        if (selection > 1) {return;}
        } // end while
     } // end for(;;)
  }

//////////////////////////////////////////////////////////////////////////
/**
 * @brief calculation of recombination rate coefficients from cross sections
 *
 * @param fselect=1 calculation yields merged-beams rate coefficient
 * @param fselect=2 calculation yields plasma rate coefficient
 *
 * The resulting rate coefficients are written into an output file.
 */

void calc_alpha(int fselect)
   {
   FELE fele;
   switch(fselect)
      {
      case 1:
	     fele = &fecool;
             break;
      case 2:
             fele = &fmaxwell;
             break;
     default: 
      	     fele = &fecool;
	     fselect=1;
      }  

   char answer[2];

   SIGMARR sigma, stest = &sigmarrtest;
   int sselect, npts = 1, theo_mode = 0;
   int number_of_levels = 1, level_number = 1;
   bool rr_flag = true;
   char theo_filename[200];  
   string marker;

   printf("\n\n Which cross section to use ?");
   printf("\n   1: RR semiclassical with Stobbe corrections (SCS)");
   printf("\n   2: SCS RR with high n integration (SCI)");
   printf("\n   3: RR quantum mechanical dipole approximation (QMD)");
   printf("\n   4: RR QMD for low n + SCS for high n (QMS)");
   printf("\n   5: binned or individual DR resonances from file");
   printf("\n   6: as function of energy (as x-, y- columns in file)");
   printf("\n                                    make a choice : ");
   scanf("%d",&sselect);

   switch(sselect)
      {
      case 1:
             sigma = &sigmarrscl;
             marker.assign("SCS");
             break;
      case 2:
             sigma = &sigmarrscl2;
             marker.assign("SCI");
             break;
      case 3:  
             sigma = &sigmarrqm; 
             marker.assign("QMD");
             break;
      case 4:  
             sigma = &sigmarraqm; 
             marker.assign("QMS");
             break;
      case 5:
             rr_flag=false;
             sigma = &sigmalorentzian;
             info_DRtheory();
             while (theo_mode==0)
               {
               printf("\n Give name of theory data file : ");
               scanf("%s",theo_filename);
               npts = open_DRtheory(theo_filename, theo_mode, number_of_levels);
               if (npts<0)
                 {
                 printf("\n File %s not found.\n",theo_filename);
                 theo_mode=0;
                 }
               }
             printf("\n %d DR cross section datasets.\n",npts); 
             if (npts==0) return;
	     if (number_of_levels>1)
	       {
		 printf("\n Data for %d different initial levels found.",number_of_levels);
		 printf("\n Give the number (in the range 1 - %d) of the level to be used: ",number_of_levels);
		 scanf("%d",&level_number);
		 if ( (level_number<1) || (level_number>number_of_levels) )
		   {
		     printf("\n\n ERROR: Level number out of range!\n\n");
		     exit(0);
		   }
	       }
             break;             
      case 6: 
             convolute_Xsec_MBrateCoeff(1);
             return;
             break;
      default:
             sigma = &sigmarrscl;
             marker.assign("SCS");
             sselect=1;
      } // end switch (sselect)

   double *energy   = new double[npts];
   double *strength = new double[npts];
   double *wLorentz = new double[npts];

   int use_fraction=0, nint=0, nmin, nmax, lmin, nele;
   double z, t=1.0;
   double emin,emax;

   if (rr_flag)
     {
     printf("\n Give z, nmin, lmin, nmax ..............................: ");
     scanf(" %lf %d %d %d",&z,&nmin,&lmin,&nmax);
     printf("\n Give number of electrons in min n,l-subshell ..........: ");
     scanf("%d",&nele);
     t = 1.0 - nele/(4.0*lmin+2.0);

     if (sselect==2) 
       {
       use_fraction = 0;
       printf("\n Give n where integration dn starts ....................: ");
       scanf("%d",&nint);
       }
     else
       {
       printf("\n Use surviving fractions from file ? (y/n) .............: ");
       scanf("%s",&answer);
       use_fraction = ( (!strcmp(answer,"y")) || (!strcmp(answer,"Y")) );
       } 
     }
   else
     {
     use_fraction = 0;
     read_DRtheory(theo_filename,theo_mode,level_number,energy,strength,wLorentz,npts,emin,emax);
     printf("\n DR resonances ranging from %8.4f eV to %8.4f eV\n",emin,emax);
     if (theo_mode==10)
       { // use interpolation for calculating cross section
	 sigma = &sigma_interpolated;
	 init_sigma_interpolation(energy,strength,npts);
	 rr_flag = true;   // interpolation is treated as rr
	 npts = 1;
	 t = 1.0;
       }
     }
 
   int fdim=4;
   if (use_fraction) 
     {
     // each n,l selective RR cross section is multiplied by fraction[nl]
     fdim = nmax*(nmax+1)/2;
     } 
   else 
     {
     // the vector 'fraction' may be used for passing up to 'fdim'
     // parameters to 'sigma'
     fdim=4;
     }
   double *fraction = new double[fdim];
   char header[200];
   if (rr_flag)
     {
     if (use_fraction) 
       {
       for (int k=0; k<=fdim; k++) {fraction[k] = 1.0;}
       for (int l=0; l<lmin; l++) {fraction[(nmin-1)*nmin/2+l] = 0.0;}
       fraction[(nmin-1)*nmin/2+lmin] = t;
       nmax = readfraction(fraction,header,nmax);
       printf("\n surviving fractions %.*s\n",sizeof(header),header);
       }
     else
       {
       fraction[0] = t;
       }
     if (theo_mode==0) {
       printf("\n z=%5.2f, nmin=%2d, lmin=%2d, nele=%2d, nmax=%3d",
	      z,nmin,lmin,nele,nmax);
       if (sselect==2) { printf(", nint=%3d\n",nint); }
       else { printf("\n");}
     }
    }
   
   double eV,edelta;
   int steps,itest;
   // in order to generate theoretical rate coefficients at exactly the
   // same energies as given in an experimental data file, the energies
   // from this file can be read in. It is assumed that energies are
   // listed in the first column of the experimental data file. 
   FILE *fenergy;
   char filename[200],line[200];
   strcpy(line,"0");
   printf("\n Read energies from file? Give filename (0=no file) ....: ");
   scanf("%s",&filename);
   int read_energy = (strcmp(filename,line));
   if (read_energy)
     {
     fenergy = fopen(filename,"r");
     if (!fenergy)
       {
       read_energy = 0;
       printf("\n file %s not found!",filename);
       }
     else
       { // count the number of energies given in the data file
       steps = 0;
       while ((itest=fgetc(fenergy))!=EOF) {fgets(line,200,fenergy); steps++;}
       steps--;
       fclose(fenergy);
      }
     }
   if (!read_energy) {
     do {
       printf("\n Give energy range (emin,emax,delta) ...................: ");
       scanf("%lf %lf %lf",&emin,&emax,&edelta);
       if (edelta>0) { // create equidistant steps
         steps = int((emax-emin)/edelta);
       } else { // create energydependant steps 
         if (emin<0) emin=0;
         double ecm=emin;
         double estep=1e-5;
         steps=0;
         while (ecm<=emax) {
           steps++;
           double estep;
           estep = -edelta * sqrt(ecm);
           if ((estep<1e-5) || (ecm<1e-5)) {estep=1e-5;}
           ecm += estep;
         }
	 cout<<endl<<" Number of steps created ........................: "<<steps<<endl;
	 if (steps>=1e6) {
	   cout<<" Number of steps to large!"<<steps<<endl;
         }
       }
     } while (steps>=1e6); 
     steps --;
   eV=emin;
   }

   double ktpar, ktperp, delta;
   int nsub, legendre_n, laguerre_n;
   double legendre_x[100], legendre_w[100], laguerre_x[100], laguerre_w[100];

   new_temperatures:

   ktpar = 0.0;
   ktperp = 0.0;
   if (fselect==1)
   {
       printf("\n Give electron beam temperatures kTpar and kTperp (meV) : ");
       scanf("%lf %lf", &ktpar,&ktperp);
       ktpar *= 0.001;   // meV -> eV
       ktperp *= 0.001;  // meV -> eV
   }

   setup_integration(ktpar, ktperp, &nsub, &delta,
                     &legendre_n, legendre_x, legendre_w,
                     &laguerre_n, laguerre_x, laguerre_w);

   FILE *fout;
   char fn[200];
   printf("\n Give filename for output ..............................: ");
   scanf("%s",&fn);
   printf("\n");

   fout = fopen(fn,"w");
   fprintf(fout,"###########################################################\n");
   if (fselect==1)
     {
       fprintf(fout,"### hydcrocal cooler recombination rate-coefficient\n");
     }
   else
     {
       fprintf(fout,"### hydcrocal plasma recombination rate-coefficient\n");
     }
   fprintf(fout,"### SVN revision       : %s\n",SVNrevision);
   fprintf(fout,"### filename           : %s\n",fn);
   if (fselect==1)
     {
       fprintf(fout,"### kTpar (meV)        : %8.3f\n",ktpar*1000.0);
       fprintf(fout,"### kTperp (meV)       : %8.3f\n",ktperp*1000.0);
     }
  if (rr_flag)
     {
       fprintf(fout,"### RR calculation type: %s\n",marker.data());
       if (use_fraction)
	 {
	   fprintf(fout,"### surviving fract.: %.*s\n",sizeof(header),header);
	 }
       fprintf(fout,"### nuclear charge z   : %5.2f\n",z);
       fprintf(fout,"### lowest occupied subshell\n");
       fprintf(fout,"###      nmin          : %d\n",nmin);
       fprintf(fout,"###      lmin          : %d\n",lmin);
       fprintf(fout,"###      nele          : %d\n",nele);
       fprintf(fout,"### highest occupied subshell\n");
       fprintf(fout,"###      nmax          : %d\n",nmax);
       if (sselect==2)
	 {
	   use_fraction = nint; // parameter nint passed instead of use_fraction
	   fprintf(fout,"###  nint            : %3\n",nint);
	 }
     }
  else     
    {
      fprintf(fout,"### DR theory file     : %s\n",theo_filename);
      fprintf(fout,"### level number       : %d\n",level_number);
    }
  fprintf(fout,"### parameters of numerical integration\n");
  fprintf(fout,"###    Legendre points : %d\n",legendre_n);
  fprintf(fout,"###    Laguerre points : %d\n",laguerre_n);
  fprintf(fout,"###########################################################\n");
  if (fselect==1) 
    {
      fprintf(fout,"###  E (eV) [1]   alpha (cm3/s) [2]   norm [3]\n");
    } 
  else 
    {
      fprintf(fout,"### kT (eV) [1]   alpha (cm3/s) [2]   norm [3]\n");
    }
  fprintf(fout,"###--------------------------------------------------------\n");
  
  if (read_energy)
    {
      fenergy = fopen(filename,"r");
     }

   for (int i=0; i<=steps; i++)
     { // loop over energies
     if (read_energy) {
       fgets(line,200,fenergy);
       sscanf(line,"%lf",&eV);
     } else {
       if (edelta>0) { // create equidistant steps
         eV = emin +i*edelta;
         if (emin<0) eV = exp(log(10.0)*eV);
       } else { // create energy dependent steps 
         double estep=1e-5;
         estep = -edelta * sqrt(eV);
         if ((estep<1e-5) || (eV<1e-5)) {estep=1e-5;}
         eV += estep;
       }
     }
     double peaks[3],normpeak[3];

     // peaks[0] stores the width of the electron energy distribution
     if (fselect==1) 
     {
	 peaks[0] = log(2.0)*ktperp+4*sqrt(log(2.0)*eV*ktpar);
     }
     else
     {
         peaks[0] = log(2.0)*eV;
     }
     peaks[1] = 0.0;
     peaks[2] = 0.0;
     for (int k=0; k<3;k++) {normpeak[k] = peaks[k];}
     double alpha = 0.0; 
     double norm = 0.0;
     for(int n=0; n<npts; n++)
       { // loop over DR resonances, in case of RR npts=1
       bool numerical_flag = true;
       if (!rr_flag) // i.e. if DR calculation 
         {
	 if (energy[n]<=0.0) continue; 
         if (fabs(wLorentz[n])<1E-9)
           { // narrow peaks are treated as delta like resonances
           numerical_flag = false;
           if (fabs(strength[n])>0) alpha += deltapeak(fele,energy[n],eV,ktpar,ktperp,fabs(strength[n]));
           norm += 1.0;
           }
         else
           { // broad peaks are convoluted numerically
	 // printf("res. no. %4d: E=%12.4g eV, S=%12.4g eVcm^2, W=%12.4g eV\n",
	 //         n,energy[n],strength[n],width[n]);
	   z = energy[n]; // pass resonance energy as parameter "z"
           // further resonance parameters are stored in the array "fraction"
           fraction[0] = strength[n];
           fraction[1] = wLorentz[n];
           // pass resonance widths and resonance strength to 
           // convolution routine via the arra "peaks"
           // used for setting up the energy mesh of the convolution
           peaks[1] = energy[n];
           peaks[2] = fabs(wLorentz[n]);
           } // end else 
         } // end if (!rr_flag)

       if (numerical_flag)
         {
	   alpha += convolute(sigma,fele,theo_mode,fabs(eV),ktpar,ktperp,peaks,z,
                            fraction,use_fraction,nmin,nmax,lmin,nsub,delta,
                            legendre_x,legendre_w,legendre_n,
                            laguerre_x,laguerre_w,laguerre_n);
	   norm += convolute(stest,fele,0,fabs(eV),ktpar,ktperp,normpeak,z,
                            fraction,use_fraction,nmin,nmax,lmin,nsub,delta,
                            legendre_x,legendre_w,legendre_n,
                            laguerre_x,laguerre_w,laguerre_n);
         }
       } // end for (n...)
     norm /= npts;
     if (alpha<1e-99) alpha=0.0;
     if (rr_flag && ((i % 100)==0)) printf("%12.6g %17.6g %12.6g\n",eV,alpha,norm);
     fprintf(fout,"%12.6g, %17.6g, %12.6g\n",eV,alpha,norm);
     } // end for (i...)
   fclose(fout);
   if (read_energy) fclose(fenergy);

   printf("\n Use other temperatures ? (y/n) : ");
   scanf("%s",&answer);
   if (!strcmp(answer,"y")) goto new_temperatures;

   delete[] fraction;
   delete[] wLorentz;
   delete[] strength;
   delete[] energy;
   if ((theo_mode==2)||(theo_mode==10)) delete_sigma_interpolation();
   }

//////////////////////////////////////////////////////////////////////////

