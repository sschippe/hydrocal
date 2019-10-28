/**
 * @file fieldion.cxx
 *
 * @brief Calculation of hydrogenic field-ionization rates and probabilities
 *                                                                         
 * @verbatim
  $Id: fieldion.cxx 446 2017-08-28 16:08:49Z iamp $
 @endverbatim
 *
 * Stefan Schippers@physik.uni-giessen.de
*/                                                                         

#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include "clebsch.h"
#include "hydromath.h"
#include "fieldion.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////
/**
 * @brief Hydrogenic rate of field ionization
 *
 * Rate of field ionization for a hydrogenic state
 * characterized by the parabolic quantum numbers N1, N2, M. 
 * On exit E0 contains the Stark shifted binding energy of the N1,N2,M-level. 
 * The formula for E0 has been derived by
 * Alliluev & Malkin, Sov.Phys.JETP 39, 627 (1974) and the formula for the field
 * inization rate by Damburg & Kolosov, JPB 12, 2637 (1979). 
 *
 * @param field electric field in V/cm
 * @param z nuclear charge
 * @param N1 parabolic quantum number n1
 * @param N2 parabolic quantum number n2
 * @param M parabolic quantum number m
 * @param E0 on exit Stark shifted binding energy
 *
 * @return field-ionization rate in 1/s
 */
double FIrate(double field, double z, int N1, int N2, int M, double *E0)
   {

   long double f    = field/Fau;

   long double m    = M<0 ? -M : M; 
   long double m2   = m*m; 
   long double m4   = m2*m2;
   long double n    = N1+N2+m+1.0; 
   long double n2   = n*n; 
   long double n3   = n2*n; 
   long double n4   = n2*n2;
   long double nz   = n/z; 
   long double nz2  = nz*nz;
   long double nz3  = nz2*nz; 
   long double nz4  = nz2*nz2; 
   long double nz7  = nz3*nz4; 
   long double nz10 = nz7*nz3;
   long double q    = N1-N2; 
   long double q2   = q*q; 
   long double q4   = q2*q2;

   long double a0 = -0.5/nz2;
   long double a1 = 1.5*nz*q;
   long double a2 = 0.0625*nz4*(-17*n2+3*q2+9*m2-19);
   long double a3 = 0.09375*nz7*q*(23*n2-q2+11*m2+39);
   long double a4 = nz10*(-5487*n4-147*q4+549*m4-1806*n2*q2+3402*n2*m2+1134*m2*q2
                     -35182*n2-5754*q2+8622*m2-16211)/1024.0;

   long double e0  = a0 + f*(a1+f*(a2+f*(a3+f*a4)));
   *E0 = double(e0);

   long double r   = -2.0*e0; r = r*sqrt(r)/f; 
   long double fac = 3*log(n)+lnf(N2)+lnf(N2+m)-(2*N2+m+1)*log(4.0*r);
   long double sum = 34.0*N2*(N2+m) + 46.0*N2 + 7.0*m2 + 23.0*m + 53.0/3.0;
   long double arg = fac + 0.25*nz3*f*sum + 2.0*r/3.0 + log(Tau);
   return double(exp(-arg));
   }

//////////////////////////////////////////////////////////////////////////

inline int AFIindex(int n, int m, int n1)
  {
  return int( (n-1)*n*(n+1)/6.0 + 0.5*(2*n-m+1)*m + n1 + 0.1 );
  }

//////////////////////////////////////////////////////////////////////////

void calcAFIarray(double field, double z, int nmax, double *AFIarray)
   {
   int n,m,n1,n2,index;
   double E0;
   for (n=1; n<=nmax; n++)
     for (m=0; m<n; m++)
       for (n1=0; n1<n-m; n1++)
          {
          index = AFIindex(n,m,n1);
          n2 = n-n1-m-1;
          AFIarray[index] = FIrate(field,z,n1,n2,m,&E0);
          }

// find out wrong solutions and replace them by a lower limit
    double AFIold, AFInew;
    for (n1=0; n1<nmax; n1++)
      for (m=0; m<nmax; m++)
        {  
        AFIold = 0.0;
        for (n=1+n1+m; n<=nmax; n++)
          {                  
          index = AFIindex(n,m,n1);
          AFInew = AFIarray[index];
          if (AFInew < AFIold) { AFIarray[index] = AFIold; }
             else { AFIold = AFInew; } 
          }
        }
   }
          
//////////////////////////////////////////////////////////////////////////

double nsurvival(double dt, double f, double z, int n, double *AFIarray)
   {
   // survival probability of a hydrogenic (nuclear charge z) n-state
   // which spends a time dt (in s) in an electric field f (in V/cm)
   double prob = 0.0;
   for (int m = 0; m<n; m++)
     {
     double mult = (m==0) ? 1 : 2;
     int n1max = n-1-m;
     for (int n1=0; n1<=n1max; n1++)
       {
       int index = AFIindex(n,m,n1);
       double rate = AFIarray[index];
       prob += mult*exp(-rate*dt);
       }
     }
   return prob/n/n;
   }

//////////////////////////////////////////////////////////////////////////

double nlsurvival(double dt, double f, double z, int n, int l, 
                  double *AFIarray)
   {
   // survival probability of a hydrogenic (nuclear charge z) n,l-state
   // which spends a time dt (in s) in an electric field f (in V/cm)
   long double l2 = (2*l+1), prob = 0.0, e0;
   for (int m = 0; m<=n; m++)
     {
     if (m>l) break;
     double mult = (m==0) ? 1 : 2;
     int n1max = n-1-m;
     for (int n1=0; n1<=n1max; n1++)
       {
       int n2 = n1max-n1;
       long double q = n1-n2;
       long double cg = ThreeJ(0.5*(n-1),0.5*(n-1),l,0.5*(m+q),0.5*(m-q),-m);
       int index = AFIindex(n,m,n1);
       long double rate = AFIarray[index];
       prob += mult*cg*cg*exp(-rate*dt);
       }
     }
   return double(prob);
   }

//////////////////////////////////////////////////////////////////////////

void testFIrate(void)
   {
   double F, z, e0, rate, *AFIarray;
   int n,n1,n2,m,nmax,AFIdim,index;
   char fn[200];
   FILE *fout;

   printf("\n as function of");
   printf("\n 1: n1");
   printf("\n 2: n2");
   printf("\n 3: m");
   printf("\n 4: field strength");
   printf("\n 5: nmax (writes FI rates up to nmax to a file)\n");
   int choice;
   printf("\n Make a choice : ");
   scanf("%d",&choice);
   
   switch (choice)
      {
      case 1:
              printf("\n Give max n1 ........................: ");
              scanf("%d",&n1);
              printf("\n Give nuclear charge z and n2 and m .: ");
              scanf("%lf %d %d",&z,&n2,&m);
              printf("\n Give electric field strength in a.u.: ");
              scanf("%lf",&F);
              for(n=0; n<=n1; n++)
                 {
                 rate = FIrate(F*Fau,z,n,n2,m,&e0)*Tau;
                 printf("n1 = %3d : level = %12.5g, rate = %12.5g\n",n,e0,rate);
                 }
              break;               
      case 2:
              printf("\n Give max n2 ........................: ");
              scanf("%d",&n2);
              printf("\n Give nuclear charge z and n1 and m .: ");
              scanf("%lf %d %d",&z,&n1,&m);
              printf("\n Give electric field strength in a.u.: ");
              scanf("%lf",&F);
              for(n=0; n<=n2; n++)
                 {
                 rate = FIrate(F*Fau,z,n1,n,m,&e0)*Tau;
                 printf("n2 = %3d : level = %12.5g, rate = %12.5g\n",n,e0,rate);
                 }
              break;               
      case 3:
              printf("\n Give max m .........................: ");
              scanf("%d",&m);
              printf("\n Give nuclear charge z and n1 and n2 : ");
              scanf("%lf %d %d",&z,&n1,&n2);
              printf("\n Give electric field strength in a.u.: ");
              scanf("%lf",&F);
              for(n=0; n<=m; n++)
                 {
                 rate = FIrate(F*Fau,z,n1,n2,n,&e0)*Tau;
                 printf("m = %3d : level = %12.5g, rate = %12.5g\n",n,e0,rate);
                 }
              break;               
      case 5:
              printf("\n Give nuclear charge and max n ......: ");
              scanf("%lf %d",&z,&nmax);
              printf("\n Give electric field strength in V/cm: ");
              scanf("%lf",&F);
              printf("\n Give filename for output ...........: ");
              scanf("%s",&fn);
              fout = fopen(fn,"w");
              fprintf(fout,"     i   n1   n2    m    n          AFI\n");
              AFIdim = int(nmax*(nmax+1)*(nmax+2)/6+0.1);
              printf("\n calculating %3d FI rates up to n=%4d ...\n\n",
                     AFIdim,nmax);
              AFIarray = new double[AFIdim];
              calcAFIarray(F, z, nmax, AFIarray);
              for (n1=0; n1<nmax; n1++)
                for (m=0; m<nmax; m++)
                   for (n=1+n1+m; n<=nmax; n++)
                    {                  
                    n2 = n-n1-m-1;
                    index = AFIindex(n,n1,m);
                    fprintf(fout,"%6d %4d %4d %4d %4d",index,n1,n2,m,n);
                    fprintf(fout," %12.5g\n",AFIarray[index]*Tau);
                }
              printf("\n FI rates writen to file %s\n\n",fn);
              fclose(fout);
              delete[] AFIarray;
              break;
      default:
              double fmin, fmax, df;
              printf("\n Give min, max and delta F...........: ");
              scanf("%lf %lf %lf",&fmin,&fmax,&df);
              printf("\n Give nuclear charge z, n1, n2 and m : ");
              scanf("%lf %d %d %d",&z,&n1,&n2,&m);
              int nmax = 1+int(0.1+(fmax-fmin)/df);
              for(n=0; n<=nmax; n++)
                 {
                 F = fmin+n*df;
                 rate = FIrate(F*Fau,z,n1,n2,m,&e0)*Tau;
                 printf("F = %12.5g : level = %12.5g, rate =  %12.5g\n",
                        F,e0,rate);
                 }
      }   

   }


//////////////////////////////////////////////////////////////////////////

void test_survival()
   {
   double z,f,dt;
   int nf;

   printf("\n Give nuclear charge Z ..................: ");
   scanf("%lf",&z);
    printf("\n Give approximate cut off quantum number");
   printf("\n   (if =0 input of E-field will follow) .: ");
   scanf("%d",&nf);
   if (nf<1)
     {
     printf("\n Give electric field in V/cm ............: ");
     scanf("%lf",&f);
     nf = int(sqrt(sqrt(pow(double(z),3)/(9*f/Fau))));
     printf("\n The approximate cut off quantum number is %3d\n",nf);
     }
   else
     {
     f = pow(double(z),3)*pow(double(nf),-4)/9*Fau;
     printf("\n The electric field is %12.5g V/cm\n",f);
     }
   printf("\n Give dwell time in s ...................: ");
   scanf("%lf",&dt);

   int lcut;   
   printf("\n summation over l up to");
   printf("\n 1 : n-1");
   printf("\n 2 : fixed l");
   printf("\n 3 : l_{1/2}"); 
   int lmode;
   printf("\n make a choice ..........................: ");
   scanf("%d",&lmode);
   if (lmode==2)
     {
     printf("\n Give cut off l-quantum number...........: ");
     scanf("%d",&lcut);
     }    

   char filenameroot[200], filename[200];
   printf("\n Give filename for output (*.spn/*.spl)..: ");
   scanf("%s",&filenameroot);
   FILE *fn, *fl;
   strcpy(filename,filenameroot);
   strcat(filename,".spn");
   fn = fopen(filename,"w");
   strcpy(filename,filenameroot);
   strcat(filename,".spl");
   fl = fopen(filename,"w");
   
   int lmax, nmax = 2*nf;
   int lmaxout = nmax-1;
   switch(lmode)
       {
        case 2: lmaxout = lcut; break;
        case 3: lmax = nmax<2 ? 0 : int(1.673*exp(0.68*log(double(nmax-1)))); break;
       default: lmaxout = nmax-1;
       }

   int AFIdim = nmax*(nmax+1)*(nmax+2)/6;
   double *AFI = new double[AFIdim];
   calcAFIarray(f,z,nmax,AFI);

   int n, l;
   double lsum_old = 2.0;
   for (n=1; n<=nmax; n++)
     {
     double lsum = 0.0;
     double wsum = 0.0;

     switch(lmode)
         {
          case 2: lmax = lcut < n ? lcut : n-1; break;
          case 3: lmax = n<2 ? 0 : int(1.673*exp(0.68*log(double(n-1)))); break;
         default: lmax = n-1;
         }

     for (l=0; l<=lmax; l++)
       {
       double pnl = nlsurvival(dt,f,z,n,l,AFI);
       fprintf(fl," %10.4g",pnl);
       lsum += (2*l+1)*pnl;
       wsum += (2*l+1);
       }
     for (l=lmax+1;l<=lmaxout;l++) fprintf(fl,"%10.4g",0.0);
     fprintf(fl,"\n");

     lsum /= wsum;
     if (lsum>lsum_old+1e-6) break;
     lsum_old=lsum;
     double pn = nsurvival(dt,f,z,n,AFI);

     printf("n = %3d: sum over l up to lmax = %10.4g",n,lsum);
     printf(" sum over all Stark states = %10.4g\n",pn); 
     fprintf(fn,"%5d %12.5g %12.5g\n",n,lsum,pn); 
     }
   delete[] AFI;
   fclose(fn);
   fclose(fl);
   printf("\n   n-selective survival probabilities written to %s.spn",
          filenameroot);
   printf("\n n,l-selective survival probabilities written in matrix");
   printf(" format to %s.spl",filenameroot);
   }

//////////////////////////////////////////////////////////////////////////

void testFI(void)
   {
   printf("\n Calculate what ?\n");
   printf("\n 1: field ionization rates");
   printf("\n 2: survival probabilities");
   
   int choice;
   printf("\n\n Make a choice ...........: ");
   scanf("%d",&choice);

   switch(choice)
      {
       case 2 : test_survival(); break;
      default : testFIrate(); break;
      }
   }


//////////////////////////////////////////////////////////////////////////

void calc_survival(double dt, double f, double z, int nmax, 
                   double *psurv, int *n_one, int *n_zero)
   {
   const double eps = 1e-4;

   cout << "\n Calculating survival probabilities \n"; cout.flush();

   int AFIdim = int(nmax*(nmax+1)*(nmax+2)/6+0.1);
   double *AFI = new double[AFIdim];
   calcAFIarray(f,z,nmax,AFI);

// The Damburg and Kolosov formula becomes inaccurate for high n and l.
// This can lead to inconsistencies in the survival probabilities. 
// For example, they can exhibit an increase for high n values which
// is unphysical. Therefore we first look for such an increase for
// each l separately and then for the l-averaged survival probabilities.
// High-n survival probabilities are set to zero, either if they exhibit an
// increase or if they are smaller than eps defined above.
 
   int n1  = nmax;
   int n0  = 0;
   double pold, pnew;
   int index, l, n, oneflag, zeroflag, correctflag;
   for (l=0;l<nmax;l++)
     {
     oneflag = 1; zeroflag = 0; correctflag = 1; pnew = 1.0;
     for (n=l+1;n<=nmax;n++)
       {
       if (correctflag)
          {
          pold = pnew;
          pnew = nlsurvival(dt,f,z,n,l,AFI);
          }
       if ( (oneflag) && (pnew<(1.0-eps)) ) 
          { 
          if (n < n1) n1 = n;
          oneflag = 0; 
          zeroflag = 1;
          } 
       if ((zeroflag) && ((pnew>pold) || (pnew<eps)) )
          { 
          if (n > n0) n0 = n;
          zeroflag = 0;
          correctflag = 0; 
          } 
       index = (n-1)*n/2 + l;
       if (correctflag) { psurv[index] = pnew; } 
                   else { psurv[index] = pold; }
       } // end for (n...)
     if ( ((l+1) % 10) == 0) { cout << "|"; } 
                        else { cout << "."; } 
     cout.flush();
     } // end for(l ...)


   double *pn = new double[nmax+1];

   for (n=1; n<=nmax; n++)
     {
     pn[n] = 0.0;
     int n12 = (n-1)*n/2;
     for (l=0; l<n; l++)
       {
       if (n < n1) psurv[n12+l] = 1.0;
       if (n > n0) psurv[n12+l] = 0.0;
       pn[n] += (2*l+1)*psurv[n12+l];
       }
     pn[n] /= (n*n);
     }

   zeroflag = 0;
   for (n=2; n<=nmax; n++)
     {
     if ( (zeroflag==0) && ((pn[n] > pn[n-1]) || (pn[n]<eps)) ) 
        { zeroflag = 1; n0 = n-1; } 
     if (zeroflag) for (l=0; l<n; l++) psurv[(n-1)*n/2+l] = 0.0; 
     } 

   printf("\n\n 1 > survival probabilities > 0 for %3d <= n <= %3d\n\n",n1,n0);
   for (n=n1; n<=n0; n++) { printf("\n n = %4d : psurv = %12.5g",n,pn[n]); }
   printf("\n");
   
   delete[] pn;
   delete[] AFI;

   *n_one = n1;
   *n_zero = n0;
   }

//////////////////////////////////////////////////////////////////////////

