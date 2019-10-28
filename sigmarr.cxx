/** 
    @file sigmarr.cxx
    @brief Collection of recombination cross sections for use in convolutions etc.
 
    @par CREATION  
    @author Stefan Schippers
    @date 1997
    
    @par VERSION
    @verbatim
    $Id: sigmarr.cxx 372 2016-02-04 18:43:22Z iamp $
    @endverbatim   
*/

#include <stdio.h>
#include <math.h>
#include <cstring>
#include <iostream>
#include "osci.h"
#include "fraction.h"
#include "sigmarr.h"
#include "hydromath.h"
#include "fele.h"

using namespace std;


// interpolation arrays
int _sigma_n = 0;
double *_sigma_x = NULL, *_sigma_y = NULL, *_sigma_dy=NULL;


//////////////////////////////////////////////////////////////////////////
/**
 * @brief initializes interpolation arrays
 *
 * @param energy array of energy values
 * @param xesc array of cross section values
 * @param npts array size
 */
void init_sigma_interpolation(double *energy, double *xsec, int npts)
{
  _sigma_n = npts;
  _sigma_x = new double[npts];
  _sigma_y = new double[npts];
  _sigma_dy = new double[npts];
  for (int i=0; i<npts; i++)
    {
      _sigma_x[i] = energy[i];
      _sigma_y[i] = xsec[i];
      //  cout << _sigma_x[i] << " " << _sigma_y[i] << endl;
    }
  spline(_sigma_n, _sigma_x, _sigma_y, _sigma_dy);
}


//////////////////////////////////////////////////////////////////////////
/**
 * @brief cleans up interpolation arrays
 */
void delete_sigma_interpolation()
{
  if (_sigma_x) delete[] _sigma_x;
  if (_sigma_y) delete[] _sigma_y;
  if (_sigma_dy) delete[] _sigma_dy;
  _sigma_n = 0;
}

//////////////////////////////////////////////////////////////////////////
/**
 * @brief returns interpolated cross section
 *
 * @param e energy 
 * @param . all other parameters are not used
 * 
 * @return interpolated cross section times energy
 * @return 0 if the energy is outside the interpolation range
 *
 */
double sigma_interpolated(double e, double z, double *fraction, 
                   int use_fraction, int nmin, int nmax, int lmin)
{
  if ((e<_sigma_x[0]) || (e > _sigma_x[_sigma_n-1])) {
    return 0;
  }
  else {
    return e*splint(e, _sigma_n, _sigma_x, _sigma_y, _sigma_dy);
  }
}

//////////////////////////////////////////////////////////////////////////

double deltapeak(FELE fele, double er, double eV, double ktpar, 
                 double ktperp, double a)
  // convolution of a delta peak like cross section
  // returns alpha_DR(E)
   {
   double vsigma = a*Clight*sqrt(2.0*er/Melectron);
   return (*fele)(eV,er,ktpar,ktperp)*vsigma;
   }

//////////////////////////////////////////////////////////////////////////

double sigmalorentzian(double e, double er, double *fraction,
                   int use_fraction, int nmin, int nmax, int lmin)
  // sigma_DR(E) ~ Lorentzian or ~ Lorentzian/E depending on the definition
  // of the parameter fraction[3]. returns E*sigma_DR(E).
  {
  double a  = fraction[0];   // peak area
  double g  = fraction[1];   // peak width
  double gg = 0.25*g*g;      // (width/2)^2
  double esigma = 0.5*fabs(a)*e/3.1415926535*fabs(g)/((e-er)*(e-er)+ gg);
  if (g<0) 
  {
      if (a<0)
      {
	  esigma *= (er+gg/er)/e;
      }
      else
      {
	  esigma *= er/e;
      }
  }
  return esigma;
  }

//////////////////////////////////////////////////////////////////////////

double sigmarrtest(double e, double z, double *fraction, 
                   int use_fraction, int nmin, int nmax, int lmin)
  // test cross section yielding alpharr = 1.0
  {   
  return 1.0/Clight/sqrt(2.0/Melectron/e);
  }


//////////////////////////////////////////////////////////////////////

double sigmarrqm(double z, int n, int l)
   { // lim e->0  e*sigma (eV cm^2)
   double nfac = z/n/n*Rydberg;
   double sigma = nfac*nfac*(2.0*l+1.0)*fosciBC(n,l);
   return 2*Pi*Pi*A0*A0*sigma/Rydberg/Alpha/Alpha/Alpha;
// Alpha is here approximately the inverse of 1/137
   }

//////////////////////////////////////////////////////////////////////

double sigmarrqm(double eV, double z, int n, int l)
   { // e*sigma (eV cm^2)
// eV= (sqrt(Melectron*Melectron+eV*eV) - Melectron);
   double e = eV/Rydberg/z/z;
   double nfac = (eV+z*z/n/n*Rydberg)/z;
   double sigma = nfac*nfac*(2.0*l+1.0)*fosciBC(e,n,l);
   return 2*Pi*Pi*A0*A0*sigma/Rydberg/Alpha/Alpha/Alpha;
// Alpha is here approximately the inverse of 1/137
   }

//////////////////////////////////////////////////////////////////////

double sigmaPIqm(double eV, double z, int nmin, int n, int l, double IP)
   { // photoionization cross section (sigma*E)
   double gf = 1;
   double gi = 4*l+2;
   double Erel = eV-13.606*z*z*(1.0/(n*n)-1.0/(nmin*nmin))+IP;
   return 1.022e6*gf/gi/eV*sigmarrqm(Erel,z,n,l);
   }

//////////////////////////////////////////////////////////////////////

double sigmaPIqm(double eV, double z, double *fraction, 
                 int use_fraction, int nmin, int nmax, int lmin, double IP)
   { // photoionization cross section (sigma*E)
   double lterm=sigmaPIqm(eV,z,nmin,nmin,lmin,IP);
   int n12=(nmin-1)*nmin/2;
   if (use_fraction) { lterm *= fraction[n12+lmin]; }
                else { lterm *= fraction[0];}
   double lsum = lterm;
   for (int l=lmin+1; l<nmin;l++)
     {
     lterm = sigmaPIqm(eV,z,nmin,nmin,l,IP);
     if (use_fraction) { lterm *= fraction[n12+l]; }
     lsum += lterm;
     }
   double nsum = lsum;
   for (int n=nmin+1; n<= nmax; n++)
      {
      lsum = 0.0;
      n12=(n-1)*n/2;
      for (int l=0; l<n; l++) 
         {
         lterm = sigmaPIqm(eV,z,nmin,n,l,IP);
         if (use_fraction) { lterm *= fraction[n12+l]; }
         lsum += lterm;
         }
      nsum += lsum;
      }
   return nsum;
   }

//////////////////////////////////////////////////////////////////////

double sigmarrqm(double eV, double z, double *fraction, 
                 int use_fraction, int nmin, int nmax, int lmin)
   { // e*sigma (eV cm^2)
// eV= (sqrt(Melectron*Melectron+eV*eV) - Melectron);
   double e = eV/Rydberg/z/z;
   double nfac = (eV+z*z/nmin/nmin*Rydberg)/z;
   double lterm = (2.0*lmin+1.0)*fosciBC(e,nmin,lmin);
   int n12 = (nmin-1)*nmin/2;
   if (use_fraction) { lterm *= fraction[n12+lmin]; }
                else { lterm *= fraction[0];}
   double lsum = lterm;
   for (int l=lmin+1; l<nmin;l++)
     {
     lterm = (2.0*l+1.0)*fosciBC(e,nmin,l);
     if (use_fraction) { lterm *= fraction[n12+l]; }
     lsum += lterm;
     }
   double nsum = nfac*nfac*lsum;
   for (int n=nmin+1; n<= nmax; n++)
      {
      nfac = (eV+z*z/n/n*Rydberg)/z;
      lsum = 0.0;
      n12=(n-1)*n/2;
      for (int l=0; l<n; l++) 
         {
         lterm = (2.0*l+1.0)*fosciBC(e,n,l);
         if (use_fraction) { lterm *= fraction[n12+l]; }
         lsum += lterm;
         }
      nsum += nfac*nfac*lsum;
      }
   return 2*Pi*Pi*A0*A0*nsum/Rydberg/Alpha/Alpha/Alpha;
   }


//////////////////////////////////////////////////////////////////////

double stobbe(int n) 
   {
   double lsum = 0.0; 
   for (int l=0; l<n; l++) { lsum += (2.0*l+1.0)*fosciBC(n,l); }
   return 3.0*sqrt(3.0)*Pi/16/n/n/n*lsum;
   }


void calc_stobbe(void)
{

    int nmax;
    FILE *fout;

    printf(" Give nmax : ");
    scanf("%d",&nmax);
    fout = fopen("Stobbecorr.dat","w");
    for (int n=0; n<= nmax; n++) fprintf(fout,"%d %g\n",n,stobbe(n));
    fclose(fout);
    printf("\n Stobbe correction factors written to Stobbecorr.dat.\n\n");
}


//////////////////////////////////////////////////////////////////////

double kstobbe(int n) 
  {
const double stobbetab[100] = 
{0.7973014601263097, 0.876185137748039, 0.907508212449926, 0.924743401639527, 
  0.93581638846388, 0.943608658751098, 0.949431439764705, 0.953971669087505, 
  0.957626026642256, 0.960640572468532, 0.963176553031456, 0.965344334771262, 
  0.967222177180456, 0.96886722316249, 0.970322239758068, 0.971619916920603, 
  0.972785701488747, 0.973839719665624, 0.974798114035793, 0.975673993973899, 
  0.976478124451407, 0.97721943395042, 0.97790539484623, 0.978542312298236, 
  0.979135546461289, 0.979689685399302, 0.980208681071948, 0.980695957327253, 
  0.981154496436549, 0.981586909013667, 0.981995490945762, 0.982382270081938, 
  0.982749044779071, 0.983097415924431, 0.983428813695211, 0.983744520043171, 
  0.984045687685223, 0.984333356221232, 0.984608465876629, 0.984871869270909, 
  0.985124341537149, 0.985366589057629, 0.985599257032808, 0.985822936062623, 
  0.986038167888209, 0.986245450417208, 0.986445242135483, 0.98663796599147, 
  0.986824012825752, 0.98700374440719, 0.987177496127627, 0.987345579399426, 
  0.987508283793644, 0.987665878951198, 0.987818616294869, 0.987966730566065, 
  0.988110441207099, 0.9882499536069, 0.988385460225757, 0.988517141612672, 
  0.988645167327187, 0.988769696776052, 0.988890879973833, 0.989008858235477, 
  0.989123764807859, 0.989235725446557, 0.989344858943334, 0.989451277609236, 
  0.989555087717599, 0.989656389910845, 0.989755279574496, 0.989851847181457, 
  0.989946178609314, 0.990038355433108, 0.990128455195772, 0.990216551658218, 
  0.990302715030843, 0.99038701218806, 0.990469506867306, 0.990550259853813, 
  0.990629329152346, 0.990706770146964, 0.990782635749768, 0.990856976539537, 
  0.990929840891035, 0.99100127509572, 0.991071323474532, 0.991140028483344, 
  0.991207430811654, 0.991273569474995, 0.991338481901563, 0.991402204013441, 
  0.991464770302852, 0.991526213903768, 0.991586566659212, 0.991645859184559, 
  0.99170412092711, 0.991761380222187, 0.991817664345996, 0.991872999565475};
  return n<= 100 ? stobbetab[n-1] : 1.0;
  }

//////////////////////////////////////////////////////////////////////

double sigmarrscl(double eV, double z, double *fraction, 
                 int use_fraction, int nmin, int nmax, int lmin)
   { // e*sigma (eV cm^2)
// eV = sqrt(Melectron*Melectron+eV*eV) - Melectron;
   double z2 = z*z*Rydberg;
   double meanfraction = 1.0;
   int n12 = (nmin-1)*nmin/2; 
   if (use_fraction)
      {
      meanfraction = 0.0;
      for(int l=0; l<nmin; l++) {meanfraction += fraction[n12+l]*(2*l+1);}
      meanfraction /= (nmin*nmin);
      }
    else 
      {
      // calculate fraction of holes of nmin-shell from fraction of holes
      // of nmin,lmin-subshell stored in fraction[0]
      meanfraction=(nmin*nmin-lmin*lmin-(1-fraction[0])*(2*lmin+1))/nmin/nmin;
      }
   double nsum = meanfraction*z2/(1.0+nmin*nmin*eV/z2)/nmin;
   if (z>0) nsum *= kstobbe(nmin);
   for (int n=nmin+1; n<=nmax; n++)
     {
     n12 = (n-1)*n/2; 
     if (use_fraction)
        {
        meanfraction = 0.0;
        for(int l=0; l<n; l++) {meanfraction += fraction[n12+l]*(2*l+1);}
        meanfraction /= (n*n);
        }
     else
        {
        meanfraction = 1.0;
        }
     double nterm = meanfraction*z2/(1.0+n*n*eV/z2)/n;
     if (z>0) nterm *= kstobbe(n);
     nsum += nterm;
     }
   return 32*Pi/Alpha/Alpha/Alpha*A0*A0/3/sqrt(3.0)*nsum;
   }

//////////////////////////////////////////////////////////////////////

double sigmarraqm(double eV, double z, double *fraction, 
                 int use_fraction, int nmin, int nmax, int lmin)
   { // e*sigma (eV cm^2)
   int nqm = 49;
   nqm = nmax > nqm ? nqm : nmax;
   double sigma = sigmarrqm(eV,z,fraction,use_fraction,nmin,nqm,lmin);
   if (nmax>nqm) 
     { 
     double ratio = Pi*3*sqrt(3.0)/16*
                    sigmarrqm(eV,z,fraction,use_fraction,nqm+1,nqm+1,lmin)/
                    stobbe(nqm+1)/
                    sigmarrscl(eV,z,fraction,use_fraction,nqm+1,nqm+1,lmin);
     sigma += ratio*sigmarrscl(eV,z,fraction,use_fraction,nqm+1,nmax,lmin);
     }
   return sigma;
   }

//////////////////////////////////////////////////////////////////////

double sigmarrscl2(double eV, double z, double *fraction, 
                   int n2, int nmin, int nmax, int lmin)
   { // e*sigma (eV cm^2)
   double sigma = sigmarrscl(eV,z,fraction,0,nmin,n2-1,lmin);
   sigma += 2.1E-22*13.6056*z*z*log((double(nmax)/double(n2)));
   return sigma;
   }

//////////////////////////////////////////////////////////////////////

double alphaRRplasmaSCL(double kt, double z, int nmin, int nmax, int nele)
{
 double z2r = Rydberg*z*z;
 double x   = z2r/kt;
 double fac = 32*Pi/Alpha/Alpha/Alpha*A0*A0/3/sqrt(3.0);
 fac *= 4*Clight*z2r*z2r/sqrt(2*Pi*kt*Melectron)/kt;
 double t   = 1-nele/(2*nmin*nmin);
 double ratecoeff = fac/pow(double(nmin),3.0)*e1exp(x/nmin/nmin)*t*kstobbe(nmin);
 for (int n=nmin+1; n<=nmax; n++)
   {
   ratecoeff += kstobbe(n)*fac/pow(double(n),3.0)*e1exp(x/n/n);
   } 
 return ratecoeff;
}  
//////////////////////////////////////////////////////////////////////

double betaRRplasmaSCL(double kt, double z, int nmin, int nmax, int nele)
{
 double z2r = Rydberg*z*z;
 double x   = z2r/kt;
 double fac = 32*Pi/Alpha/Alpha/Alpha*A0*A0/3/sqrt(3.0);
 fac *= 4*Clight*z2r*z2r/sqrt(2*Pi*kt*Melectron)/kt;
 double t   = 1-nele/(2*nmin*nmin);
 double ratecoeff = fac/pow(double(nmin),3.0)*(1-x/nmin/nmin*e1exp(x/nmin/nmin))*t*kstobbe(nmin);
 for (int n=nmin+1; n<=nmax; n++)
   {
   ratecoeff += kstobbe(n)*fac/pow(double(n),3.0)*(1-x/n/n*e1exp(x/n/n));
   } 
 return ratecoeff;
}  

//////////////////////////////////////////////////////////////////////

void calc_sigmaRR(void)
  {
  double z, emin, emax, edelta, eV, t=0, sigma, IP;
  int n, nint, l, nmin, lmin, nmax, nele, choice=0;
  int fdim=2, use_fraction=0;
  char answer[2],filename[200],fnroot[200],header[200];
  FILE *fout, *fout2;

  printf("\n Select type of RR calculation ?\n");
  printf("\n  1:  semiclassical calculation with Stobbe corrections (SCS)");
  printf("\n  2:  semicalssical calculation without Stobbe corrections (SCL)");
  printf("\n  3:  quantum mechanical calculation, dipole approximation (QMD)");
  printf("\n  4:  QMD just one n,l for a range of energies");
  printf("\n  5:  QMD, all n,l up to n=nmax at fixed energy");
  printf("\n  6:  SCS up to nint-1, from nint to nmax integral dn (SCI)");
  printf("\n  7:  photoionization from RR QMD");
  while ((choice<1) || (choice>7)) 
    {
    printf("\n\n Make a choice .................................. : ");
    scanf("%d",&choice);
    }
  printf("\n Give effective nuclear charge ...................: ");
  scanf("%lf",&z);
  z = fabs(z);
  if (choice==2) {z = -fabs(z); } // don't use Stobbe corrections 
  if (choice==4)
    {
    printf("\n Give n,l ........................................: ");
    scanf("%d %d",&nmin,&lmin);
    }
  else if (choice==5)
    {
    printf("\n Give maximum main quantum number nmax ...........: ");
    scanf("%d",&nmax);
    printf("\n Give electron energy in eV ......................: ");
    scanf("%lf",&eV);
    }
  else
    {
    printf("\n Give nmin, lmin, nmax ...........................: ");
    scanf("%d %d %d",&nmin,&lmin,&nmax);
    printf("\n Give number of electrons in minimum n,l shell ...: ");
    scanf("%d",&nele);
    t = 1.0-nele/(4.0*lmin+2.0); /* fraction of holes in min n,l subshell */

    if (choice==6) 
      {
      use_fraction = 0;
      printf("\n Give n where integration dn starts ..............: ");
      scanf("%d",&nint);
      }
    else if (choice==7) 
      {
      use_fraction = 0;
      printf("\n Give first ionization potential in eV ...........: ");
      scanf("%lf",&IP);
      }
    else
      {
      printf("\n Use surviving fractions from file ? (y/n) .......: ");
      scanf("%s",&answer);
      if ( (!strcmp(answer,"y")) || (!strcmp(answer,"Y")) ) 
        {use_fraction=1; fdim=nmax*(nmax+1)/2;}
      } 
   }

  printf("\n Give filename for output (*.sig) ................: ");
  scanf("%s",&fnroot);
  strcpy(filename,fnroot);
  strcat(filename,".sig");
  fout = fopen(filename,"w");
  
  if (choice==5)
    {
    strcpy(filename,fnroot);
    strcat(filename,".snl");
    fout2 = fopen(filename,"w");
    fprintf(fout2,"    n   lmax lmean lhalf interpolated\n");
    double *es = new double[nmax];
    double *esmean = new double[nmax];
    int l, choice2;
    printf("\n Which kind of output ?");
    printf("\n 1: sigma(n,l) in cm^2");
    printf("\n 2: sigma(n,l) times energy in eV cm^2");
    printf("\n 3: sigma(n,l)/sigma(n)");
    printf("\n 4: sigma(n,l) relative to maximum value per n");
    printf("\n Make a choice ...................................: ");
    scanf("%d",&choice2);
    for(int n=1; n<=nmax;n++)
      {
      double esmax = 0.0, estot = 0.0, esout;
      int lmax = 0;
      for(l=0; l<n; l++)
        {
        es[l] = sigmarrqm(eV,z,n,l);
        estot += es[l];
        esmean[l] = estot/(l+1.0);
        if (es[l]>esmax) { lmax = l; esmax = es[l]; }
        }
      switch (choice2)
        {
         case 2: for(l=0; l<n; l++) 
                   {
                   esout = es[l] > 1e-99 ? es[l] : 1e-99;
                   fprintf(fout,"%12.6g",esout);
                   }
                 break;
         case 3: for(l=0; l<n; l++) {fprintf(fout,"%12.6g",es[l]/estot);}
                 break;
         case 4: for(l=0; l<n; l++) {fprintf(fout,"%12.6g",es[l]/esmax);}
                 break;
        default: for(l=0; l<n; l++) 
                   {
                   esout = es[l]/eV > 1e-99 ? es[l]/eV : 1e-99;
                   fprintf(fout,"%12.6g",esout);
                   }
                 break;
        }
      for(int k=n; k<nmax; k++)
        {
        fprintf(fout,"%12.6g",1e-99);
        }
      fprintf(fout,"\n");

      int lhalf = 0, lmean = 0;
      for (int ll=n-1; ll>=0; ll--)
        {
        if ((lhalf==0) && (es[ll]>0.5*es[lmax]))   { lhalf=ll; }
        if ((lmean==0) && (es[ll]>0.5*esmean[ll])) { lmean=ll; }
        }

      double ilhalf = lhalf > 0 ? 
               lhalf + (0.5*es[lmax]-es[lhalf])/(es[lhalf+1]-es[lhalf]) : 0;
      fprintf(fout2,"%5d %5d %5d %5d %12.5g\n",n,lmax,lmean,lhalf,ilhalf);
      } // end for(n...)
    fclose(fout);
    fclose(fout2);
    printf("\n Quantum mechanical RR cross sections at E = %10.3g eV",eV);
    printf("\n written in matrix form to file %s.sig.\n",fnroot);
    delete[] es;
    delete[] esmean;
    return;
    } 

  double *fraction = new double[fdim];
  if (use_fraction) 
    {
    for (int k=0; k<= fdim; k++) {fraction[k] = 1.0;}
    for (l=0; l< lmin; l++)  {fraction[(nmin-1)*nmin/2+l] = 0.0;}
    fraction[(nmin-1)*nmin/2+lmin] = t;
    nmax=readfraction(fraction,header,nmax);
    printf("\n surviving fractions %.*s\n",sizeof(header),header);
    fprintf(fout,"surviving fractions %.*s\n",sizeof(header),header);
    }
  else
    {
    fraction[0] = t;
    }

  switch(choice)
     {
     case 1: fprintf(fout,"semiclassical RR cross section");
             fprintf(fout,"  with Stobbe correction\n");   
             fprintf(fout,"SCS: z=%5.2f, nmin=%3d, lmin=%3d",z,nmin,lmin);
             fprintf(fout," nele=%3d, nmax=%3d\n",nele,nmax);
             break;
     case 2: fprintf(fout,"semiclassical RR cross section");
             fprintf(fout,"  without Stobbe correction\n");   
             fprintf(fout,"SCL: z=%5.2f, nmin=%3d, lmin=%3d",z,nmin,lmin);
             fprintf(fout," nele=%3d, nmax=%3d\n",nele,nmax);
             break;
     case 3: fprintf(fout,"quantum mechanical RR cross section");
             fprintf(fout," (dipole approximation)\n"); 
             fprintf(fout,"QMD: z=%5.2f, nmin=%3d, lmin=%3d",z,nmin,lmin);
             fprintf(fout," nele=%3d, nmax=%3d\n",nele,nmax);
             break;
     case 4: fprintf(fout,"quantum mechanical RR cross section");
             fprintf(fout," (dipole approximation)\n"); 
             fprintf(fout,"QMD: z=%5.2f, n=%3d, l=%3d\n",z,nmin,lmin);
             break;
     case 6: fprintf(fout,"semiclassical RR cross section");
             fprintf(fout,"  with high-n integration\n"); 
             fprintf(fout,"SCI: z=%5.2f, n=%3d, l=%3d\n",z,nmin,lmin);
             fprintf(fout," nele=%3d, nint=%3d, nmax=%3d\n",nele,nmax,nint);
             break;
     case 7: fprintf(fout,"quantum mechanical PI cross section");
             fprintf(fout,"  (dipole approximation) \n"); 
             fprintf(fout,"SCI: z=%5.2f, n=%3d, l=%3d\n",z,nmin,lmin);
             fprintf(fout," nele=%3d, nmax=%3d, IP=f6.3%\n",nele,nmax,IP);
             break;
      }
  
   int itest, steps;
   // in order to generate theoretical cross sections at exactly the
   // same energies as given in an already existing file, the energies
   // from this file can be read in. It is assumed that energies are
   // listed in the first column of the file. 
   FILE *fenergy;
   char line[200]="0";
   printf("\n Read energies from file? Give filename (0=no file): ");
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
   if (!read_energy)
     {
     printf("\n Give energy range (emin,emax,delta) ......................: ");
     scanf("%lf %lf %lf",&emin,&emax,&edelta);
     steps = int((emax-emin)/edelta);
     }

  printf("\n     Ecm (eV)        sigma (cm^2)    E*sigma (eV cm^2)\n"); 
  fprintf(fout,"     Ecm             sigma           E*sigma\n"); 
  fprintf(fout,"     (eV)            (cm^2)          (eV cm^2)\n");

  if (read_energy) fenergy = fopen(filename,"r");

  for (int i=0; i<= steps; i++)
    {
     if (read_energy)
       {
       fgets(line,200,fenergy);
       sscanf(line,"%lg",&eV);
       // cout<<eV<<endl;
       }
     else
       {
       eV = emin +i*edelta;
       if (emin<0) eV = exp(log(10.0)*eV);
       }

    switch(choice)
       {
        case 1: sigma = sigmarrscl(eV,z,fraction,use_fraction,nmin,nmax,lmin); 
                break;
        case 2: sigma = sigmarrscl(eV,z,fraction,use_fraction,nmin,nmax,lmin); 
                break;
        case 3: sigma = sigmarrqm(eV,z,fraction,use_fraction,nmin,nmax,lmin); 
                break;
        case 6: sigma = sigmarrscl2(eV,z,fraction,nint,nmin,nmax,lmin); 
                break;
        case 7: sigma = sigmaPIqm(eV,z,fraction,use_fraction,nmin,nmax,lmin,IP); 
                break;
       default: sigma = sigmarrqm(eV,z,nmin,lmin); 
                break;
       }
    printf(" %12.5G     %15.5G      %15.5G\n",eV,sigma/eV,sigma);
    fprintf(fout," %12.5G     %15.5G      %15.5G\n",eV,sigma/eV,sigma);
    }
  if (read_energy) fclose(fenergy);
  fclose(fout);
  }

//////////////////////////////////////////////////////////////////////

void compare_sigmarr(void)
  {
  int timesenergy;
  int nmin, nmax, lmin, nele;
  double z, emin, emax, edelta, eV, t, sigma_qm, sigma_scl, ratio;
  char answer[2], filename[200];
  FILE *fout;
  printf("\n Comparison of QMD and SCS RR cross sections\n");
  printf("\n QMD: quantum mechanical dipole approximation");
  printf("\n SCS: semiclassical with Stobbe corrections\n");
  printf("\n Give effective nuclear charge ................: ");
  scanf("%lf",&z);
  printf("\n Give min n,l and max n (nmin,lmin,nmax) ......: ");
  scanf("%d %d %d",&nmin,&lmin,&nmax);
  printf("\n Give number of electrons in minimum n,l shell : ");
  scanf("%d",&nele);
  t = 1.0-nele/(4.0*lmin+2.0); /* fraction of holes in min n,l subshell */
  printf("\n Give filename for output .....................: ");
  scanf("%s",&filename);
  fout = fopen(filename,"w");
  fprintf(fout,"z=%4.1f, nmin=%2d, lmin=%2d, nele=%2d, nmax=%3d\n",
          z,nmin,lmin,nele,nmax);


  char header[200];
  int fdim, use_fraction;
  printf("\n Use surviving fractions from file ? (y/n) ....: ");   
  scanf("%s",&answer);
  if ( (!strcmp(answer,"y")) || (!strcmp(answer,"Y")) ) 
    {use_fraction=1; fdim=nmax*(nmax+1)/2;}
  else 
    {use_fraction=0; fdim = 2;}

  double *fraction = new double[fdim];
  if (use_fraction) 
     {
     for (int k=0; k<= fdim; k++) {fraction[k] = 1.0;}
     for (int l=0; l< lmin; l++) {fraction[(nmin-1)*nmin/2+l] = 0.0;}
     fraction[(nmin-1)*nmin/2+lmin] = t;
     nmax=readfraction(fraction,header,nmax);
     printf("\n surviving fractions %.*s\n",sizeof(header),header);
     fprintf(fout,"surviving fractions %.*s\n",sizeof(header),header);
     }
  else
     {
     fraction[0] = t;
     }

   int itest, steps;
   // in order to generate theoretical cross sections at exactly the
   // same energies as given in an already existing file, the energies
   // from this file can be read in. It is assumed that energies are
   // listed in the first column of the file. 
   FILE *fenergy;
   char line[200]="0";

   printf("\n Read energies from file? Give filename (0=no file): ");
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
   if (!read_energy)
     {
     printf("\n Give energy range (emin,emax,delta) ...................: ");
     scanf("%lf %lf %lf",&emin,&emax,&edelta);
     steps = int((emax-emin)/edelta);
     }
  printf("\n Multiply sigma times energy ? (y/n) ..........: ");
  scanf("%d",&answer);
  if ( (!strcmp(answer,"y")) || (!strcmp(answer,"Y")) ) 
    {
    timesenergy = 1;
    printf("\n\n     Ecm (eV)  sQMD (eV cm^2)  sSCS (eV cm^2)    sSCS/sQMD\n"); 
    fprintf(fout,"     Ecm (eV)  sQMD (eV cm^2)  sSCS (eV cm^2)    sSCS/sQMD\n"); 
    }
  else
    {
    timesenergy = 0;
    printf("\n\n     Ecm (eV)     sQMD (cm^2)     sSCS (cm^2)    sSCS/sQMD\n"); 
    fprintf(fout,"     Ecm (eV)     sQMD (cm^2)     sSCS (cm^2)    sSCS/sQMD\n"); 
    }

  if (read_energy) fenergy = fopen(filename,"r");

  for (int i=0; i<= steps; i++)
    {
     if (read_energy)
       {
       fgets(line,200,fenergy);
       sscanf(line,"%lg",&eV);
       }
     else
       {
       eV = emin +i*edelta;
       if (emin<0) eV = exp(log(10.0)*eV);
       }

    sigma_qm = sigmarrqm(eV,z,fraction,use_fraction,nmin,nmax,lmin);
    sigma_scl = sigmarrscl(eV,z,fraction,use_fraction,nmin,nmax,lmin);
    ratio = sigma_scl/sigma_qm;
    sigma_qm = timesenergy ? sigma_qm : sigma_qm/eV;
    sigma_scl = timesenergy ? sigma_scl : sigma_scl/eV;
    printf(" %12.5G %15.5G %15.5G %12.5G\n",eV,sigma_qm,sigma_scl,ratio);
    fprintf(fout," %12.5G %15.5G %15.5G %12.5G\n",eV,sigma_qm,sigma_scl,ratio);
    }
  if (read_energy) fclose(fenergy);
  fclose(fout);
  delete[] fraction;
  }

//////////////////////////////////////////////////////////////////////

void calcAlphaRRplasma(void)
  {
  double z;
  int nmin,nmax,nele;
  char filename[200];
  FILE *fout, *ftemp;

  printf("\n Give z, nmin, nmax ........................................: ");
  scanf("%lf %d %d",&z,&nmin,&nmax);
  printf("\n Give number of electrons already in nmin-shell ............: ");
  scanf("%d",&nele);

  int steps, itest;
  char line[200];
  printf("\n Read temperatures from file? Give filename (0=no file) ....: ");
  scanf("%s",&filename);
  strcpy(line,"0");
  int read_temp = (strcmp(filename,line));
  if (read_temp)
    {
    ftemp = fopen(filename,"r");
    if (!ftemp)
      {
      read_temp = 0;
      printf("\n file %s not found!",filename);
      }
    else
      { // count the number of energies given in the data file
      steps = 0;
      while ((itest=fgetc(ftemp))!=EOF) {fgets(line,200,ftemp); steps++;}
      steps--;
      fclose(ftemp);
      }
    }

  double kt,ktmin,ktmax,ktdelta;
  if (read_temp)
    {
    ftemp = fopen(filename,"r");
    }
  else
    {
    printf("\n Give log min, max, delta kT in eV  ........................: ");
    scanf("%lf %lf %lf",&ktmin,&ktmax,&ktdelta);  
    steps = int((ktmax-ktmin)/ktdelta);
    }


  printf("\n Give filename for output ..................................: ");
  scanf("%s",&filename);
  fout = fopen(filename,"w");
  fprintf(fout,"z=%4.1f, nmin=%2d, nele=%2d, nmax=%3d\n",
          z,nmin,nele,nmax);

  printf("\n      kT (eV)      alpha (cm^3/s)        kT^3/2*alpha         beta (cm^3/s)      kT^(5/2)*beta\n"); 
  fprintf(fout,"      kT (eV)      alpha (cm^3/s)        kT^3/2*alpha          beta (cm^3/s)     kT^5/2*beta\n"); 

  for(int i=0;i<=steps;i++)
    {
    if (read_temp)
      { 
      fgets(line,200,ftemp);
      sscanf(line,"%lg",&kt);
      }
    else
      { 
      kt = ktmin+i*ktdelta;
      kt = exp(log(10.0)*kt);
      }
    double alpha = alphaRRplasmaSCL(kt,z,nmin,nmax,nele);
    double beta = betaRRplasmaSCL(kt,z,nmin,nmax,nele);
    double akt = alpha*kt*sqrt(kt);
    double bkt = beta*kt*kt*sqrt(kt);
    printf(" %12.5G     %15.5G      %15.5G     %15.5G     %15.5G\n",kt,alpha,akt,beta,bkt);
    fprintf(fout," %12.5G     %15.5G      %15.5G     %15.5G     %15.5G\n",kt,alpha,akt,beta,bkt);
    }
  fclose(fout);
 }

//////////////////////////////////////////////////////////////////////








