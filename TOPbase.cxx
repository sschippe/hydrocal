// $Id: TOPbase.cxx 301 2014-12-01 20:10:11Z iamp $
/** \file TOPbase.cxx
 *  \brief Handling of data from the TOPbase atomic data base
 */
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include "fele.h"
#include "hydromath.h"

using namespace std;

double FWHM_ElHassan_Fe2(double eV);
double FWHM_ElHassan_Fe3(double eV);
double FWHM_ElHassan_Fe4(double eV);

void convolute_TOPbase_PI(void)
{
  char answer[2], line[200], filename[200];
  //istringstream buf(line,istringstream::in);
  FILE *fin, *fout;
  int Znucl, Nele;
  string asymb[]={"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si",
		  "P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe"}; //26
  string lsymb[]={"S","P","D","F","G","H","I","K","L","M","N","O","Q"}; //13
  string psymb[]={"e","o"};
  
  printf("\n Convolution of photoionization (PI) cross sections from the TOPbase.");
  printf("\n Cross section files can be retrieved from the internet, see");
  printf("\n       http://cdsweb.u-strasbg.fr/topbase/ftp.html.\n");
  printf("\n The cross sections can be selected by term.");
  printf("\n Cross sections can be written to a file before or after convolution.");
  printf("\n A gaussian width can be specified for the convolution.");
  printf("\n During the convolution the TOP cross sections are linearly interpolated.\n\n");
  printf("\n Give name of TOBbase data file containing PI cross sections: ");
  scanf("%s",filename);
  fin = fopen(filename,"r");
  if (!fin) {
    printf("\n ERROR::File %s not found.\n\n",filename);
    return;
  }
  
  if (fgets(line,200,fin)==NULL)
    {
    printf("\n ERROR::File %s is empty.\n\n",filename);
    return;
    }
  sscanf(line,"%d %d", &Znucl, &Nele);
  int Qcharge = Znucl-Nele;
  
  printf("\n nuclear charge .....: %2d",Znucl);
  printf("\n number of electrons.: %2d",Nele);
  if (Znucl > 26)
    {
      printf("\n ERROR: Nuclear charge > 26.\n\n");
      return;
    }
  printf("\n ion ................: %s%d+\n",asymb[Znucl-1].c_str(),Qcharge);
    
  unsigned int Ssearch, Lsearch, Psearch, Nsearch;
  printf("\n Give 2S+1, L, parity (0 or 1), and term number (1,2,...) ..: ");
  scanf("%u %u %u %u", &Ssearch, &Lsearch, &Psearch, &Nsearch);
  if (Psearch>1)
    {
      printf("\n ERROR: invalid parity value, please specifiy 0 or 1.");
      return;
    }
  if (Lsearch>12)
    {
      printf("\n ERROR: invalid L value, please specifiy one of 0,1,...,12.");
      return;
    }
       
  unsigned int Sread, Lread, Pread, Nread;
  int nlines, ndummy;
  double *ecm=NULL, *xsec=NULL;
  double elo = 1E99, ehi = -1E99;
  int search_flag = 1;
  while(search_flag)
    {
      if (fgets(line,200,fin)==NULL) break;
      sscanf(line,"%d %d %d %d", &Sread, &Lread, &Pread, &Nread);
      printf("\n Reading cross section for term %d %d %d %d ...",Sread,Lread,Pread,Nread);
      if (fgets(line,200,fin)==NULL) break;
      sscanf(line,"%d %d", &ndummy, &nlines);
      if (fgets(line,200,fin)==NULL) break;

      ecm = new double[nlines];
      xsec = new double[nlines];
      elo = 1E99, ehi = -1E99;
      for(int n=0; n<nlines; n++)
	{
	  fgets(line,200,fin);
	  sscanf(line,"%lf %lf", &(ecm[n]), &(xsec[n]));
	  //printf("%20G %20G\n",ecm[n],xsec[n]);
	  ecm[n] = ecm[n]*13.6057;
	  if (ecm[n]<elo) elo = ecm[n];
	  if (ecm[n]>ehi) ehi = ecm[n];
	}
      if ( (Sread==Ssearch) && (Lread==Lsearch) && (Pread==Psearch) && (Nread==Nsearch) )
	{
	  printf("\n\n %d cross section values found.\n",nlines);
	  search_flag = 0;
	}
      else
	{
	  delete[] xsec;
	  delete[] ecm;
	}
    }
  fclose(fin);
  
  if (search_flag) {
    printf("\n ERROR::Term not found.\n\n");
    return;
  }
  char term_string[32];
  sprintf(term_string,"%d%s%s%d",Ssearch,lsymb[Lsearch].c_str(),psymb[Psearch].c_str(),Nsearch);
  
  FELE fele=&felegauss;
  double fwhm=0.0;
  int fwhm_mode = 0;
  printf("\n Specify photon energy resolution mode");
  printf("\n  0: unconvoluted cross section");
  printf("\n  1: constant resolution");
  printf("\n  2: relative resolution");
  printf("\n  3: specific photon energy dependence");
  printf("\n Make your choice ........................................: ");
  scanf("%d",&fwhm_mode);
  
  if (fwhm_mode==0)
    { // output of unconvoluted cross section data
      printf("\n Output of raw data or interpolated data (r/i) ...........: ");
      scanf("%s",&answer);
      int raw_flag = ((!strcmp(answer,"r")) || (!strcmp(answer,"R")));
      
      printf("\n Give name of output file ................................: ");
      scanf("%s",filename);
      fout = fopen(filename,"w");
      if (raw_flag) {
	fprintf(fout,"# TOPbase raw PI cross sections\n");
      }
      else {
	fprintf(fout,"# TOPbase interpolated PI cross sections\n");
      }
      fprintf(fout,"# %s%d+ %s (statistical weight: %d)\n",asymb[Znucl-1].c_str(),Qcharge,term_string,Ssearch*(2*Lsearch+1));
      fprintf(fout,"      Erel           s%s\n",term_string);
      fprintf(fout,"      (eV)           (Mb)\n");
      if (raw_flag)
	{
	  for(int n=0; n<nlines; n++)
	    {
	      fprintf(fout,"%15.8G   %15.8G\n",ecm[n],xsec[n]);
	    }
	  printf("\n Raw cross sections written to %s",filename);	  
	}
      else
	{
	  double emin, emax, edelta;
	  printf("\n Input energy range in eV: [%f, %f]\n",elo,ehi);
	  printf("\n Give output energy range (emin,emax,delta) ..............: ");
	  scanf("%lf %lf %lf",&emin,&emax,&edelta);
	  for (double eV = emin ; eV <= emax ; eV += edelta)
	    {
	      double xsec_interp = 0.0;
	      if ( (eV>elo) && (eV<ehi) )
		{
		  int n=0;
		  for (int nn=1; nn<nlines; nn++) {
		    if (ecm[nn]>eV) {
		      n = nn;
		      break;
		    }
		  }
		  if (n>0) {
		    xsec_interp = xsec[n-1]+(xsec[n]-xsec[n-1])/(ecm[n]-ecm[n-1])*(eV-ecm[n-1]);
		  }
		}  
	      fprintf(fout,"%15.8G   %15.8G\n",eV,xsec_interp);
	    }
	  printf("\n Linearly interpolated cross sections written to %s",filename);
	}
      fclose(fout);
      delete[] xsec;
      delete[] ecm;
      return;
    }
  else if (fwhm_mode==1)
    {
      printf("\n Give constant fwhm in eV ................................: ");
      scanf("%lf",&fwhm);
    }
  else if (fwhm_mode==2)
    {
      printf("\n Give relative fwhm (proportional to photon energy) ......: ");
      scanf("%lf",&fwhm);
    }
  else
    {
      printf("\n Which fwhm energy dependence to use?");
      printf("\n 2: El Hassan et al., PRA 79 (2009) 033415: Fe2+");
      printf("\n 3: El Hassan et al., PRA 79 (2009) 033415: Fe3+");
      printf("\n 4: El Hassan et al., PRA 79 (2009) 033415: Fe4+");
      int choice;
      printf("\n Make your choice ........................................: ");
      scanf("%d",&choice);
      if ((choice<2) || (choice>4)) {
	printf("\n ERROR::Illeagal choice.\n\n");
	return;
      }    
      fwhm_mode = 300+choice;
    }
  
  printf("\n Give name of output file ................................: ");
  scanf("%s",filename);

  double emin, emax, edelta, eV;
  int steps;
  printf("\n Input energy range in eV: [%f, %f]\n",elo,ehi);
  printf("\n Give output energy range (emin,emax,delta) ..............: ");
  scanf("%lf %lf %lf",&emin,&emax,&edelta);
  steps = int((emax-emin)/edelta)+1;     

  printf("\n Convoluted cross sections will be calculated at %d energies.\n",steps);
  printf("\n Progress: ");
  fout = fopen(filename,"w");
  fprintf(fout,"# TOPbase convoluted PI cross sections\n");
  fprintf(fout,"# %s%d+ %s (statistical weight: %d)\n",asymb[Znucl-1].c_str(),Qcharge,term_string,Ssearch*(2*Lsearch+1));
  fprintf(fout,"      Erel           fwhm        s%s\n",term_string);
  fprintf(fout,"      (eV)           (eV)        (Mb)\n");
  
  for (int i=0; i<steps; i++)
    { // loop over energies
      if ( (i % 100) == 0) {
	printf(".");
	fflush(stdout);
      }
      eV = emin +i*edelta;
      double conv_fwhm = fwhm;
      if (fwhm_mode==2) {
        conv_fwhm = fwhm*eV;
      }
      else if (fwhm_mode==302) {
	conv_fwhm = FWHM_ElHassan_Fe2(eV);
      }
      else if (fwhm_mode==303) {
	conv_fwhm = FWHM_ElHassan_Fe3(eV);
      }
      else if (fwhm_mode==304) {
	conv_fwhm = FWHM_ElHassan_Fe4(eV);
      }
      double conv_e; // integration variable     
      double conv_width = 10.0*conv_fwhm; // left and right from the gaussian
      double conv_delta = 0.05*conv_fwhm; // step width for stepping through the gaussian
      double conv_emin = eV-conv_width;
      double conv_emax = eV+conv_width;

      if ((conv_emax<elo) || (conv_emin>ehi))
	{
	  fprintf(fout,"%15.8G   %15.8G   %15.8g\n",eV,conv_fwhm,0.0);
	  continue;
	}
	
      // make conv_delta smaller than the minimum raw energy difference 
      for (int n = 0; n<(nlines-1); n++)
	{
	  if ( (ecm[n] > conv_emin) && (ecm[n] < conv_emax) ) {
	    if (0.1*(ecm[n+1]-ecm[n]) < conv_delta) conv_delta = 0.1*(ecm[n+1]-ecm[n]);
	  }
	}
      
      double conv_sum = 0.0;
      for (conv_e = conv_emin; conv_e <= conv_emax; conv_e += conv_delta)
	{
	  double xsec_interp = 0.0;
	  if ( (conv_e>elo) && (conv_e<ehi) )
	    {
	      int n=0;
	      for (int nn=1; nn<nlines; nn++) {
		if (ecm[nn]>conv_e)
		  {
		    n = nn;
		    break;
		  }
	      }
	      if (n>0) {
		xsec_interp = xsec[n-1]+(xsec[n]-xsec[n-1])/(ecm[n]-ecm[n-1])*(conv_e-ecm[n-1]);
	      }
	    }  
	  conv_sum += conv_delta*fele(eV,conv_e,conv_fwhm,conv_fwhm)*xsec_interp;
	}
      
      if (fabs(conv_sum)<1E-99) conv_sum = 0.0;
      fprintf(fout,"%15.8G   %15.8G   %15.8G\n",eV,conv_fwhm,conv_sum);
    }   
  fclose(fout);
  printf("\n Convoluted cross sections written to %s",filename);
  printf("\n\n");
  delete[] xsec;
  delete[] ecm;
}

double FWHM_ElHassan_Fe2(double eV)
{
  if (eV<30.0) {
    return 0.1;
    }
  else if ((eV >= 30.0) && (eV<=60.0)) {
    return 0.1+0.77/30.0*(eV-30.0);
  }
  else if ((eV > 60.0) && (eV<=100)) {
    return 0.24+0.63/40.0*(eV-60.0);
  }
  else if ((eV > 100.0) && (eV<=160.0)) {
    return 1.67+3.301/60.0*(eV-100.0);
  }
  else {
    return 4.98;
  }
}

double FWHM_ElHassan_Fe3(double eV)
{
  if (eV<30.0) {
    return 0.06;
    }
  else if ((eV >= 30.0) && (eV<=45.0)) {
    return 0.06+0.19/15.0*(eV-30.0);
  }
  else if ((eV > 45.0) && (eV<=65.0)) {
    return 0.13+0.30/20.0*(eV-45.0);
  }
  else if ((eV > 65.0) && (eV<=100.0)) {
    return 0.30+0.60/35.0*(eV-65.0);
  }
  else if ((eV > 100.0) && (eV<=160.0)) {
    return 1.67+3.301/60.0*(eV-100.0);
  }
  else {
    return 4.98;
  }
}

double FWHM_ElHassan_Fe4(double eV)
{
  if (eV<60.0) {
    return 0.13;
    }
  else if ((eV >= 60.0) && (eV<=80.0)) {
    return 0.13+0.14/20.0*(eV-60.0);
  }
  else if ((eV > 80.0) && (eV<=140.0)) {
    return 0.49+1.38/60.0*(eV-80.0);
  }
  else {
    return 1.87;
  }
}

