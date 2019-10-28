/**
 * @file drmodel.cxx
 *
 * @brief Calculation of model DR spectra from parameterized atomic rates
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: drmodel.cxx 446 2017-08-28 16:08:49Z iamp $
   @endverbatim
 *
 */
#include <stdio.h>
#include <math.h>
#include <cstring>
#include "clebsch.h"
#include "fele.h"
#include "fraction.h"
#include "lifetime.h"

using namespace std;

void calcdr(void)
  {
  const double ryd=13.6058;  // eV
  const double ln2=0.6931472;
  const double mec2 = 510999.6; // eV
  const double clight = 2.99792458e10; // cm/s
  const int maxn = 1000; // highest n+1 

  char answer[2];
  double qdef[maxn],q,slim,arI,arII,aa0,ld1,ld2,ld3,ld4,factor,q2=1.0,q4=1.0;
  int model,hydrogenic=0,field=0,l,ncore,lcore,nstart,lstart=0,nmin,nmax,defmax;

  printf("\n Which model to use?");
  printf("\n       nl-dependent Auger rates (with DRF): 1");
  printf("\n   only n-dependent Auger rates   (no DRF): 2");
  printf("\n                            make a choice : ");
  scanf("%d",&model);

  
  printf("\n Give ion charge ..........................: ");
  scanf("%lf",&q);
  if (q<1) { q = -q; q2 = q*q; q4 = q2*q2; }

  printf("\n Give series limit in eV ..................: ");
  scanf("%lf",&slim);
  slim *= q2;
  if (model==2)
    {
    printf("\n Give quantum defect ......................: ");
    scanf("%lf",&qdef[0]);
    defmax = 0;
    }
  else
    {
    l=0;
    printf("\n Different quantum defects can be spefied");
    printf("\n for each l separately. From some l onwards");
    printf("\n all quantum defects will be the same. The");
    printf("\n input can then be terminated by giving -1.");
    printf("\n If the input for the l=0 quantum defect");
    printf("\n is -1, then all quantum defects will be 0.\n");
    do
      {
      printf("\n Give quantum defect for l =%2d (-1 quits) .: ",l);
      scanf("%lf",&qdef[l]);
      }
    while (qdef[l++]>=0);
    defmax=l-2;
    }
  if (defmax<0) {qdef[0]=0.0; defmax=0;}

  double emin,emax,edelta,ktpar,ktperp,norm; 
  printf("\n Give energy range (min, max, delta in eV) : ");
  scanf("%lf %lf %lf",&emin,&emax,&edelta);
  emax *= q2; // q2 is 1 if on input q>0 or q^2 if on input q<0

  printf("\n Give ktpar and ktperp in meV .............: ");
  scanf("%lf %lf",&ktpar,&ktperp);
  ktperp *= 0.001;
  ktpar *= 0.001;
  int i,epts = 1+int((emax-emin)/edelta);
  double *energy  = new double[epts];
  double *esigma  = new double[epts];
  double *esigmaF = new double[epts];

  new_model:

  printf("\n Give constant factor S0 in eV2 cm2 s .....: ");
  scanf("%lf",&factor);
    
  printf("\n Give core transition rate ArI in 1/s .....: ");
  scanf("%lf",&arI);
  arI *= q4; // q4 is 1 if on input q>0 or q^4 if on input q<0

  printf("\n Radiative rates of type II transitions are");
  printf("\n modelled as ArII/n^3");
  if (model==1) 
    { 
    printf("\n If a negative value is given for ArII then");
    printf("\n hydrogenic rates will be used instead.");
    printf("\n  -ArII is then interpreted as the number");
    printf("\n of fully occupied subshells in the ionic"); 
    printf("\n core, e.g. core = 1s2 2s: ArII = -1.");
    printf("\n Give ArII = -0.5 if no subshell is fully"); 
    printf("\n occupied (bare or hydrogenlike core).\n");
    printf("\n Give ArII in 1/s (if <0 hydrogenic rates).: ");
    }  
  else
    {
    printf("\n Give ArII in 1/s .........................: ");
    }

  scanf("%lf",&arII);
  arII *= q4; // q4 is 1 if on input q>0 or q^4 if on input q<0
  if ((arII<0)&&(model==1))  
    {
    hydrogenic = 1;
    if (arII>-0.9) 
      { ncore = 0; lcore = 0; }
    else
      {
      // -arII = 0.5*ncore*(ncore-1)+lcore
      // calculate ncore and lcore from -arII
      double shell = sqrt(-2*arII+0.25)-0.5;
      ncore = int(shell+0.1);
      lcore = int(-arII - ncore*(ncore+1)/2-1+0.1);
      if (lcore >= ncore-1) { nstart = ncore+1; lstart=0; }
                       else { nstart = ncore; lstart=lcore+1;}
      }
    }

  if (model==1)
    {
    printf("\n Auger rates are modeled as");
    printf(" Aa(n,l) = Aa0*exp(d1*l+d2*l^2+d3*l^3+d4*l^4)/n^3");
    printf("\n Give Aa0 (in 1/s), d1, d2, d3, d4 .......: ");
    scanf("%lf %lf %lf %lf %lf",&aa0,&ld1,&ld2,&ld3,&ld4);
    }
  else    
    {
    printf("\n Auger rates are modelled as Aa0/n^3 (up to lmax)");
    printf("\n Give Aa0 (in 1/s) ........................: ");
    scanf("%lf",&aa0);
    ld1=0.0; ld2=0.0; ld3=0.0;
    }

  printf("\n Give nmin ................................: ");
  scanf("%d",&nmin);

  new_nmax:
  printf("\n Give nmax ................................: ");
  scanf("%d",&nmax);


  int nAa=0;
  if ((model==1)&&(defmax==0))
    {
    printf("\n Include electric field mixing ? (y/n) ....: ");
    scanf("%s",&answer);
    field = (!strcmp(answer,"y"));
    }
  if (field)
    {
    printf("\n Write field mixed Auger rates Aa(n,k,m) to file ?");
    printf("\n Give one n value larger than 0 ..........:  ");
    scanf("%d",&nAa); 
    }


  int read_fractions, fracdim = nmax*(nmax+1)/2;
  double *fraction = new double[fracdim];
  char header[200];
  for(i=0; i<fracdim; i++) fraction[i] = 1.0;
  printf("\n Read surviving fractions from file ? (y/n): ");
  scanf("%s",&answer);
  if (!strcmp(answer,"y")) readfraction(fraction,header,nmax);  
 
  FILE *fout, *fAnkm;
  char filename[200],fn[200],fnAnkm[200];
  printf("\n Give filename for output (*.spe) .........: ");
  scanf("%s",&fn);
  strcpy(filename,fn);
  strcat(filename,".str");
  fout = fopen(filename,"w");

  fprintf(fout,"kTpar                : %10.5f meV\n",ktpar*1000);
  fprintf(fout,"kTperp               : %10.5f meV\n",ktperp*1000);
  fprintf(fout,"ion charge state     : %10.5f\n",q);
  fprintf(fout,"core radiative rate  : %10.4g /s \n",arI);
  if (hydrogenic)
    {
    fprintf(fout,"type II rad. rates   : hydrogenic\n");
    fprintf(fout,"ion core filled up to: n = %2d, l = %2d\n",ncore,lcore);
    fprintf(fout,"highest final shell  : n = %2d\n",nmin-1);
    }
  else
    { 
    fprintf(fout,"type II rad. rate    : %10.4g /s \n",arII);
    }
  fprintf(fout,"Aa0                  : %10.4g /s \n",aa0);
  if (model == 1)
    {
    fprintf(fout,"l-decay coefficient 1: %10.4g\n",ld1);
    fprintf(fout,"l-decay coefficient 2: %10.4g\n",ld2);
    fprintf(fout,"l-decay coefficient 3: %10.4g\n",ld3);
    fprintf(fout,"l-decay coefficient 4: %10.4g\n",ld4);
    }
  fprintf(fout,"factor               : %10.4g eV2 cm2 s\n",factor);
  fprintf(fout,"series limit         : %10.5f eV\n",slim);
  for(l=0; l<defmax; l++)
     fprintf(fout,"q-defect for l  = %2d : %10.5f\n",l,qdef[l]);  
  fprintf(fout,"q-defect for l >= %2d : %10.5f\n",defmax,qdef[defmax]);
  fprintf(fout,"\n\n In the follwing listed are:");
  fprintf(fout,"\n       n: main quantum number");
  fprintf(fout,"\n   E(eV): resonance energy in eV");
  fprintf(fout,"\n      sn: DR resonance strength in eVcm^2");
  fprintf(fout,"\n  Sum sn: sum of sn up to actual n");
  fprintf(fout,"\n     nDR: number of states participating in DR");
  if (field)
    {
    fprintf(fout,"\n     sFn: DR resonance strength at maximum field mixing");
    fprintf(fout,"\n Sum sFn: sum of sFn up to actual n");
    fprintf(fout,"\n    nDRF: number of mixed states participating in DR");
    fprintf(fout,"\n    enh1: nDRF/nDR");
    fprintf(fout,"\n    enh2: sn/sFn");
    }
  fprintf(fout,"\n");
    
  fprintf(fout,"\n    n      E(eV)         sn     Sum sn    nDR");
  if (field) { fprintf(fout,"        sFn    Sum sFn   nDRF   enh1   enh2\n");}
        else { fprintf(fout,"\n"); }
  fclose(fout);


// l-dependent part of the autoionization rates

  double *lexpdecay = new double[nmax];

  int ldecayflag=1; double ldecay=0.0, ldecayold;
  for (l=0; l<nmax;l++)
     { 
     if (l>defmax) qdef[l] = qdef[defmax];
     ldecayold = ldecay;
     ldecay = (((ld4*l+ld3)*l+ld2)*l+ld1)*l;
     if ((ldecayflag) && 
          ( (ldecay<-30) || ( (l>10) && (ldecay>ldecayold)) )) ldecayflag = 0;
     lexpdecay[l] = (ldecayflag) ? exp(ldecay) : -1.0;
     }


  if (nAa)
    {  // ouput of Stark mixed autonionization rates for n = nAa 
    strcpy(fnAnkm,fn);
    strcat(fnAnkm,".Ankm");
    fAnkm = fopen(fnAnkm,"w");

    double aa, j1 = 0.5*(nAa-1);
    for (int m = 0; m<nAa; m++)
      {
      int k;
      for (k= -nAa+1;k<-nAa+m+1;k++) fprintf(fAnkm," --");
      int kk=1;
      for (k= -nAa+m+1; k<= nAa-m-1; k++) 
        { // the quantum number k is counted in steps of 2
        if (kk)
          {
          aa=0;
          for (l=m;l<nAa;l++)  
            { 
            if (lexpdecay[l]<=0.0) break;
            double m1 = 0.5*(m-k);
            double m2 = 0.5*(m+k);
            double cg = CG(j1,j1,m1,m2,l,m);
            aa += cg*cg*lexpdecay[l];
            }
          aa *= aa0/pow(double(nAa),3);
          fprintf(fAnkm," %10.4g",aa);
          kk = 0;
          }
        else
          { // every other k does not occur
          fprintf(fAnkm," --");
          kk = 1;
          }
        } // end for (k...)
        for (k= nAa-m;k<nAa;k++) fprintf(fAnkm," --");
        fprintf(fAnkm,"\n");
      } // end for (m...)
    fclose(fAnkm);
    } // end if(nAa)


// initialize spectra
  for(i=0; i<epts; i++)
    {
    energy[i] = emin+i*edelta;
    esigma[i] = 0.0;
    esigmaF[i] = 0.0;
    }  

// calculate spectra
  double aa,ar=0.0,e,es,statistic,lsum=0.0,lsumF=0.0,lsum2=0.0,nsum=0.0,nsumF=0.0,fsum=0.0;
  
  int arIIflag = hydrogenic;
  for (int n=nmin; n<=nmax; n++)
    {
    printf(" n = %4d",n);

    double n3 = n*n*n;

    double *arIIhydro = new double[n];
    double arIImax = 0; // maximum hydrogenic rate with one n manifold
    if (hydrogenic)
      {
      for (l=0;l<n;l++)
        {
        arIIhydro[l] = 0.0;
        if (arIIflag)
          {
          for(int n2=nstart;n2<nmin;n2++)
            {
            for (int l2=l-1;l2<=l+1;l2+=2)
              {
              if ((n2==nstart) && (l2<lstart)) break;
              if ((l2>=0) && (l2<n2)) arIIhydro[l] += hydrotrans(q,n,l,n2,l2);
              }
            }
          if (arIIhydro[l] > arIImax) arIImax = arIIhydro[l];
          } // end if (arIIflag) 
        } // end for (l...)
      if (arIImax < 0.01*arI) arIIflag = 0;
      }
    else //not hydrogenic
      {
      ar = arI+arII/n3;
      }
    if (arIIflag) printf("   arIImax = %12.6g /s",arIImax);
    printf("\n");

    int countF = 0;
    int count0 = 0;
    if (model==1)
      { 
      lsum = 0;
      lsumF = 0;
      lsum2 = 0;
      for (l=0; l<n; l++)
        {
        e=slim-ryd*q*q/(n-qdef[l])/(n-qdef[l]);
        if (e<0) continue;

        if (field) 
          { // attention: l has the meaning of m now
          double j1 = 0.5*(n-1);
          for (int k=n-l-1; k>=0; k-=2) 
            {
            aa=0; double arIIF = 0;
            for (int j=l;j<n;j++)  
              { // attention: j has the meaning of l now
              if (lexpdecay[j]<=0.0) break;
              double l1 = 0.5*(l-k);
              double l2 = 0.5*(l+k);
              double cg = CG(j1,j1,l1,l2,j,l);
              aa += cg*cg*lexpdecay[j];
              if (arIIflag) arIIF += cg*cg*arIIhydro[j];
              }
            ar = arI + arIIF;
            aa *= aa0/n3;
            int mult1 = l==0 ? 1 : 2;
            int mult2 = k==0 ? 1 : 2;
            lsumF += mult1*mult2*factor*ar*aa/(ar+aa);
            if (aa>ar) countF += mult1*mult2;
            } // end for (k...)
          } // end if (field)
        if (lexpdecay[l]<=0.0) break;
        aa = aa0*lexpdecay[l]/n3;
        if (hydrogenic) 
           { ar = arI + arIIhydro[l]; }
        else
           { ar = arI + arII/n3; }
        if (aa>ar) count0 += 4*l+2;
        es = factor*(2*l+1)*ar*aa/(ar+aa)*fraction[(n-1)*n/2+l];
        if (l<defmax)
          { // quantum defects are mutually different for l <= defmax
          for (i=0; i<epts; i++) 
            {
            esigma[i] += es*fecool(energy[i],e,ktpar,ktperp);
            }
          lsum2 += es;
          } // end if (l<defmax) 
        lsum += es;
	}  // end for(l...)
      for (i=0; i<epts; i++) 
        {
        double fe = fecool(energy[i],e,ktpar,ktperp);
        esigma[i] += (lsum-lsum2)*fe;
        if (field) esigmaF[i] += lsumF*fe;
        }
      } // end (model==1)
    else
      { // (model==2)
      e=slim-ryd*q*q/(n-qdef[0])/(n-qdef[0]);
      aa = aa0/n3;
      es = factor*aa*ar/(aa+ar);
      for (i=0; i<epts; i++) 
        {
        esigma[i] += es*fecool(energy[i],e,ktpar,ktperp);
        }
      lsum=es;    
      delete[] arIIhydro;
      } // end for (n...) 
    nsum += lsum;
    nsumF += lsumF;    
    if (e>0) 
      {
      fout = fopen(filename,"a");
      fprintf(fout,"%5d %10.4f %10.4g %10.4g %6d",n,e,lsum,nsum,count0);
      if (field)
        {
        double enh1 = count0> 0 ? double(countF)/(double(count0)) : 0.0;
        double enh2 = lsum > 0 ? lsumF/lsum : 0.0;
        fprintf(fout," %10.4g %10.4g %6d %6.2f %6.2f",
                lsumF,nsumF,countF,enh1,enh2);
        }
      fprintf(fout,"\n");
      fclose(fout);
      }
    }
  printf("\n List of line strengths written to file %s\n",filename);  

  strcpy(filename,fn);
  strcat(filename,".spe");
  fout = fopen(filename,"w");
  for (i=0; i<epts; i++)
    {
    if (energy[i]==0) continue;
    double alphafac = sqrt(2.0/energy[i]/mec2)*clight;
    double alpha = esigma[i] > 1e-99 ? alphafac*esigma[i] : 0.0;
    fprintf(fout,"%12.6g %12.6g",energy[i],alpha);
    if (field)
      {
      double alphaF = esigma[i] > 1e-99 ? alphafac*esigmaF[i] : 0.0;
      double enh = alpha > 0 ? alphaF/alpha : 0.0;
      fprintf(fout," %12.6g %6.2f",alphaF,enh);
      }
    fprintf(fout,"\n");
    }
  fclose(fout);
  printf("\n DR rate coefficient written to file %s\n",filename);  

  delete[] lexpdecay;
  delete[] fraction;

  printf("\n New nmax ? (y/n) .........................: ");
  scanf("%s",&answer);
  if (!strcmp(answer,"y")) goto new_nmax;

  printf("\n New model parameters ? (y/n) .............: ");
  scanf("%s",&answer);
  if (!strcmp(answer,"y")) goto new_model;

  delete[] esigmaF;
  delete[] esigma;
  delete[] energy;
  }


double HFAugerfactor(double I,  double j1core, double l1Rydberg, 
		     double F1, double j2core, double F2)
{
    // see Eq.9 of Pindzola et al., Phys. Rev. A 45 (1992) R7659
    double sum = 0.0;
    for(double j1 = fabs(l1Rydberg-0.5); j1<= l1Rydberg+0.6; j1 += 1.0)
    for(double J = fabs(j1-F1); J<= F1+j1; J += 1.0)
    for(double j2 = fabs(F2-J); j2<= F2+J+0.1; j2 += 1.0)
    {
//	printf("%4.1f %4.1f %4.1f\n",j1,J,j2);
	double ninej = NineJ(I,j1core,F1,j1core,l1Rydberg,j1,F2,j2,J);
	sum += (2*j1+1)*(2*j2+1)*ninej*ninej;
    }
    return 4*sum*(2*F1+1)*(2*F2+1)*pow(2*l1Rydberg+1,-2);
}

void HyperfineDR(void)
{
    double I, j1core, l1Rydberg, j2core;

    printf("\n Give I, j1core, l1Rydberg, j2core: ");
    scanf(" %lf %lf %lf %lf",&I,&j1core,&l1Rydberg,&j2core);

    double testsum=0;
    for (double F1=fabs(I-j1core); F1 <= I+j1core+0.1; F1 += 1.0)
    for (double F2=fabs(I-j2core); F2 <= I+j2core+0.1; F2 += 1.0)
    {
	testsum += HFAugerfactor(I, j1core, l1Rydberg, F1, j2core, F2);
    }

    printf("\n testsum = %10.4g\n",testsum);
}









