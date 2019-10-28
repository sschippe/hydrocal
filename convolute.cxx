/**
 * @file convolute.cxx
 *
 * @brief Convolutions of tabulated cross sections and merged-beams rate coefficients
 *
 * @author Stefan Schippers
 * @verbatim
   $Id $
 @endverbatim
 *
 * The convolution is with
 * - either an isotropic Maxwellian yield a plasma rate coefficient as function of electron temperature [convolute_alpha()]
 * - or a flattened or gaussian energy distribution yielding a convoluted cross section and a rate  coefficient [convflat()].
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include "hydromath.h"
#include "sigmarr.h"
#include "fele.h"

using namespace std;

/**
  * @brief Convolution yielding a plasma rate coefficient
  *
  * The convolution is carried out numerically.Input data do not have to be spaced equidistantly.
  */
void convolute_PlasmaRateCoeff(void)
  {

  char line[200], filename[200];
  istringstream buf(line,istringstream::in);
  FILE *fin, *fout, *ftemp;
  int n=0, choice=0, itest;
  double elo=0.0, ehi=0.0;
  printf("\n Convolution of cross sections with a Maxwellian electron energy");
  printf("\n distribution yielding the rate coefficient in a plasma as a");
  printf("\n function of electron temperature kT (in eV).\n");

  printf("\n Cross section data are read in from files containing at least");
  printf("\n two columns. The first one is assumed to contain energies");
  printf("\n in eV and the second either cross sections in cm^2 or rate");
  printf("\n coefficients in cm^3/s. Before being convoluted the latter");
  printf("\n are converted in to cross sections by division by the electron");
  printf("\n velocity. The convolution is carried out numerically. Energies");
  printf("\n do not have to be spaced equidistantly.\n");

  printf("\n The calculation of semiclassical RR rate coefficients with");
  printf("\n Stobbe corrections from the corresponding cross section");
  printf("\n is carried out analytically without numerical convolution.\n");

  printf("\n Convolute what?");
  printf("\n   1: rate coefficients from file");
  printf("\n   2: cross sections from file (electron ion collisions)");
  printf("\n   3: semiclassical RR cross section with Stobbe corrections");
  printf("\n   4: cross sections from file (nulear reactions) \n");
  while ( (choice<1) || (choice>4) )
    {
    printf("\n Make a choice .............................................: ");
    scanf("%d",&choice);
    }

  double m1 = Melectron; // projectile mass
  double mm = 1.0;       // mass ratio (reduced mass)/m1;
  int xsec_mode = 0;
  switch (choice)
  {
      case 4:
        xsec_mode = 1;
	double m2;
	printf("\n Give target mass in atomic mass units ....: ");
	scanf("%lf",&m2);
	printf("\n Give projectile mass in atomic mass units : ");
	scanf("%lf",&m1);
	mm = m2/(m1+m2); // mass ratio (reduced mass)/m1;
	m1 *= 931.494e6; // conversion to eV
	printf("\n Give name of file containing");
	printf(" energies (in eV) and cross sections ..: ");
	scanf("%s",filename);
	break;
      case 3:  calcAlphaRRplasma();
	return;
	break;
     case 2:
       xsec_mode = 1;
       printf("\n Give name of file containing");
       printf(" energies (in eV) and cross sections ..: ");
       scanf("%s",filename);
       break;
     case 1:
       printf("\n Give name of file containing");
       printf(" energies (in eV) and rate coefficients: ");
       scanf("%s",filename);
       break;
  }

  fin = fopen(filename,"r");
  while ((itest=fgetc(fin))!=EOF) {fgets(line,200,fin); n++;}
  fclose(fin);
  int nlines = n;
  double* ecm = new double[nlines];
  double* alpha = new double[nlines];
  fin = fopen(filename,"r");
  for(n=0;n<nlines;n++)
     {
     fgets(line,200,fin);
     sscanf(line,"%lf %lf", &(ecm[n]), &(alpha[n]));
 //    buf >> ecm[n] >> alpha[n];
     printf("%10.4G %10.4G\n",ecm[n],alpha[n]);
     }
  fclose(fin);
  printf("\n %5d lines read !\n",nlines);
  // set cross sections to 0 until last negative entry
  for(n=nlines;n>=1;n--)
     {
	    if (alpha[n-1]<0) { break; }
     }
  int negative_idx = n;
  printf("First non-zero cross section energy: %10.4G\n", ecm[negative_idx]);
  for(n=0;n<negative_idx;n++)
     {
		alpha[n] = 0;
	 }
  if (ehi<=elo) {elo = ecm[0]; ehi = ecm[nlines-1];}
  double tmin,tmax,dt,t,kt;
  printf("\n Give name of output file ..................................: ");
  scanf("%s",filename);
  fout = fopen(filename,"w");
  fprintf(fout,"         T           kT         alpha alpha*kt^3/2\n\n");

  int steps;
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

  if (read_temp)
    {
    ftemp = fopen(filename,"r");
    }
  else
    {
    printf("\n Give log min, max, delta kT in eV (separated by spaces) ...: ");
    scanf("%lf %lf %lf",&tmin,&tmax,&dt);
    steps = int((tmax-tmin)/dt);
    }

  for(int i=0;i<=steps;i++)
    {
    if (read_temp)
      {
      fgets(line,200,ftemp);
      sscanf(line,"%lf",&kt);
//      line >> kt;
//      cout<<kt<<endl;
//      buf.seekg(0);
      }
    else
      {
      t = tmin+i*dt;
      kt = exp(log(10.0)*t);
      }
    double delta = 0.5*(ecm[1]-ecm[0]);
    double f = sqrt(ecm[0])*exp(-ecm[0]*mm/kt);
    if (xsec_mode) f*=sqrt(ecm[0]);
    double sum = ecm[0]<elo ? 0.0 : delta*f*alpha[0];
    for (n=1;n<nlines-1;n++)
      {
      if ((elo<ecm[n]) && (ecm[n]<ehi))
        {
        delta = 0.5*(ecm[n+1]-ecm[n-1]);
        f = sqrt(ecm[n])*exp(-ecm[n]*mm/kt);
        if (xsec_mode) f*=sqrt(ecm[n]);
        sum += delta*f*alpha[n];
        }
      }
    if (ecm[nlines-1]<=ehi)
      {
      delta = 0.5*(ecm[nlines-1]-ecm[nlines-2]);
      f = sqrt(ecm[nlines-1])*exp(-ecm[nlines-1]*mm/kt);
      if (xsec_mode) f*=sqrt(ecm[nlines-1]);
      sum += delta*f*alpha[nlines-1];
      }
    sum *= 2.0*mm*sqrt(mm/Pi);
    if (xsec_mode) sum*=sqrt(2.0/m1)*Clight;
    fprintf(fout,"%10.4G   %10.4G    %10.4G   %10.4G\n",
            kt/8.617347e-5,kt,sum/kt/sqrt(kt),sum);
    }
  fclose(fout);
  if (read_temp) fclose(ftemp);
  delete[] alpha;
  delete[] ecm;
  }

/***************************************************************************/
/**
 * @brief Convolution yielding a convoluted cross section and a merged-beams rate coefficient
 *
 * @param data_type = 1: convolute cross section data
 * @param data_type = 2: convolute rate coefficient data
 *
 * The convolution is carried out numerically. Input data do not have to be spaced equidistantly.
 */
void convolute_Xsec_MBrateCoeff(int datatype) {
  if (datatype==1) {
    printf("\n Convolution of cross sections contained as x-, y- columns");
    printf("\n separated by a <space> in a file ");
  }
  if (datatype==2) {
    printf("\n Derivation of an experimental cross section from a rate coefficient");
    printf("\n contained as <space> separated  x-, y- columns in a file and convolution");
    printf("\n of the cross sections ");
  }
  printf("with a flattened Maxwellian, with a Gaussian, or a");
  printf("\n normalized trapezoid electron energy distribution yielding");
  printf("\n a convoluted cross section and a merged beams rate coefficient.");
  printf("\n \n");

  if (datatype==1) {
    printf("\n Cross section data are read in from files containing at least");
  }
   if (datatype==2) {
    printf("\n Rate coefficient data are read in from files containing at least");
  }
  printf("\n two columns. The first one is assumed to contain energies");
  if (datatype==1) {
      printf("\n and the second cross sections.\n");
  }
   if (datatype==2) {
    printf("\n and the second rate coefficients.\n");
  }

  printf("\n The convolution is carried out numerically. Energies");
  printf("\n do not have to be spaced equidistantly.\n");

  printf("\n The output file contains four columns: ");
  printf("\n Energy (eV), convoluted cross section (cm^2),");
  printf("\n rate coeffient (cm^3/s), and convolution width (eV).\n");

  int choice;
  printf("\n Convolute with?");
  printf("\n   1: flattened Maxwellian (or isotropic with FWHM ~ E^0.5)");
  printf("\n   2: normalized Gaussian (FWHM energy independent)");
  printf("\n   3: normalized trapezoid");
  printf("\n      Make a choice ......................................: ");
  scanf("%d",&choice);

  FILE *fin, *fout, *fenergy;
  char filename[200], filenameout[200];
  char line[200];
  strcpy(line,"0");
  double e_unit, cs_unit, tolerance;
  int fin_exists = 0;

  while(!fin_exists)
    {
    printf("\n Unconvoluted cross sections are now read from a file.");
    printf("\n ATTENTION: The file must not contain empty lines!");
    printf("\n Give name of input file..................................: ");
    scanf("%s",filename);
    fin = fopen(filename,"r");
    if (fin) fin_exists = 1;
    }

  int nlines=0, n, itest;
  char *pch;
  while ((itest=fgetc(fin))!=EOF)
    {
      fgets(line,200,fin);
      if (itest=='#') { // print comment lines
	printf("# %s",line);
      }
      else {
	nlines++;
      }
    }
  fclose(fin);
  //nlines --;

  printf("\n Give energy unit in eV ..................................: ");
  scanf("%lf",&e_unit);
  printf("\n Give cross section unit in cm^2 .........................: ");
  scanf("%lf",&cs_unit);

  printf("\n\n If the x-axis spacing is too coarse the convolution is bound to fail.");
  printf("\n In this case the convoluted cross section is substitued by the original.");
  printf("\n cross section. The associated quality of the convolution is monitored by");
  printf("\n checking whether the normalization integral of the distribution");
  printf("\n function evaluates to 1. The user can allow for tiny deviations from 1 by");
  printf("\n making an appropriate choice for the 'tolerance'. A reasonable value would");
  printf("\n be 1e-4. The check can be switched off by setting 'tolerance <= 0.0'.\n");

  printf("\n Give tolerance (<=0 means no check) .....................: ");
  scanf("%lf",&tolerance);

  double elo = 1E99, ehi = -1E99;
  double* ecm = new double[nlines];
  double* xsec = new double[nlines];

  fin = fopen(filename,"r");
  for(n=0;n<nlines;n++)
    {
      fgets(line,200,fin);
      if (!strncmp(line,"##",1)) { n--; continue; } // skip comment lines
      istringstream buf(line,istringstream::in);
      buf >> ecm[n] >> xsec[n];
      //cout << ecm[n] << "   " << xsec[n]<< endl;
      if (ecm[n]<0) {
	nlines--;
	n--;
	continue;
      }
      ecm[n] *= e_unit;    // energies are now in eV
      if (datatype==2) {
	xsec[n] /= (2.99792458e10*sqrt(2.0*ecm[n]/510.99906e3));
      }
      xsec[n] *= cs_unit;  // cross sections are now in cm^2
      if (ecm[n]<elo) elo = ecm[n];
      if (ecm[n]>ehi) ehi = ecm[n];
    }
  fclose(fin);
  printf("\n Unconvoluted cross sections read at %d energies",nlines+1);
  printf("\n                             ranging from %G to %G eV.\n",elo,ehi);

  // calculate second derivative for spline interpolation below
  double* dderiv = new double[nlines];
  spline(nlines,ecm,xsec,dderiv);

  
  FELE fele=NULL;
  double fwhm=0.0, ktpar, ktperp;
  switch (choice)
  {
  case 1:
    printf("\n Give parallel electron beam temperature in meV ..........: ");
    scanf("%lf",&ktpar);ktpar*=0.001;
    printf("\n Give perpendicular electron beam temperature in meV .....: ");
    scanf("%lf",&ktperp);ktperp*=0.001;
    fele = &fecool;
    break;
  case 2:
    printf("\n Give fwhm in eV .........................................: ");
    scanf("%lf",&ktpar);
    fele = &felegauss;
    break;
  case 3:
    printf("\n Give base width in eV ...................................: ");
    scanf("%lf",&ktpar);
    printf("\n  Give top width in eV ...................................: ");
    scanf("%lf",&ktperp);
    fele = &trapezoid;
    break;
    }

  double eV,emin,emax,edelta;
  int steps;
  // in order to generate theoretical rate coefficients at exactly the
  // same energies as given in an experimental data file, the energies
  // from this file can be read in. It is assumed that energies are
  // listed in the first column of the experimental data file.
  printf("\n Read energies from file? Give filename (0=no file) ......: ");
  strcpy(line,"0");
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
       while ((itest=fgetc(fenergy))!=EOF) {
         fgets(line,200,fenergy);
         steps++;
       }
       steps--;
       fclose(fenergy);
      }
    }
  if (!read_energy)
    {
    printf("\n Give energy range (emin,emax,delta) .....................: ");
    scanf("%lf %lf %lf",&emin,&emax,&edelta);
    steps = int((emax-emin)/edelta);
    }
  printf("\n Convoluted cross sections will be calculated at %d energies.\n",steps);
  printf("\n Give name of output file ................................: ");
  scanf("%s",filenameout);
  printf("\n Progress: ");
  fout = fopen(filenameout,"w");
  fprintf(fout,"      Erel           sigma           alpha                  FWHM\n");
  fprintf(fout,"      (eV)           (cm^2)          (cm^3/s)               (eV)\n");

  if (read_energy) {
    fenergy = fopen(filename,"r");
     if (fenergy == NULL) perror ("Error opening file");
  }

  for (int i=0; i<=steps; i++) { // loop over energies
    if ( (i % 100) == 0) printf(".");
    if (read_energy)
      {
      fgets(line,200,fenergy);
      sscanf(line,"%lf",&eV);
      }
    else
      {
      eV = emin +i*edelta;
      if (emin<0) eV = exp(log(10.0)*eV);
      }

    double fpar, fperp;
    switch (choice)
      {
      case 1:
	fpar = 4*sqrt(ktpar*eV*log(2.0));
        fperp = log(2.0)*ktperp;
        fwhm = sqrt(fpar*fpar+fperp*fperp);
        break;
      case 2:
	fwhm = ktpar;
	break;
      case 3:
	fwhm = 0.5*(ktpar+ktperp);
        break;
      }

    double delta = ecm[1]-ecm[0];
    double df = delta*fele(eV,ecm[0],ktpar,ktperp);
    double sumdf  = ecm[0]<elo ? 0.0 : df;
    double sumdfx = ecm[0]<elo ? 0.0 : df*xsec[0];

    for (n=1;n<nlines-1;n++)
      {
      if ((elo<ecm[n]) && (ecm[n]<ehi))
        {
        delta = 0.5*(ecm[n+1]-ecm[n-1]);
	df = delta*fele(eV,ecm[n],ktpar,ktperp);
        sumdf  += df;
        sumdfx += df*xsec[n];
        }
      }
    
    if (ecm[nlines-1]<=ehi)
      {
      delta = ecm[nlines-1]-ecm[nlines-2];
      df = delta*fele(eV,ecm[n],ktpar,ktperp);
      sumdf  += df;
      sumdfx += df*xsec[nlines-1];
      }

    // Here we check wether the integration yields a normalized distribution function.
    // This is not the case if the cross-section grid is too coarse, i.e. when the
    // width of the distribution function is smaller than the difference between
    // adjacent energy points of the cross-section grid. In this case we substitute
    // the convoluted cross section by the cross section itself.
    if ((tolerance > 0.0) && (fabs(sumdf-1.0) > tolerance))
      {
	// cubic spline interpolation
	
	sumdfx = splint(eV,nlines,ecm,xsec,dderiv);
	// linear interpolation
	/*
	for (n=1; n<nlines; n++)
	  {
	    if ((ecm[n-1]<= eV) && (eV < ecm[n])) break;
	  }
	sumdfx = xsec[n-1] + (eV-ecm[n-1]) *(xsec[n]-xsec[n-1])/(ecm[n]-ecm[n-1]);
	*/
      }

    
    if (fabs(sumdfx)<1E-99) sumdfx = 0.0;
    double rate_coeff = sumdfx*sqrt(2*eV/510999.0)*2.99792458e10;
    fprintf(fout,"%15.8G   %15.8G   %15.8G   %15.8G\n",eV,sumdfx,rate_coeff,fwhm);
  } // end for(i...)
  fclose(fout);
  if (read_energy) fclose(fenergy);
  printf("\n\n");
  delete[] xsec;
  delete[] ecm;
}

