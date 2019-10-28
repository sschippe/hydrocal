/**
 * @file populations.cxx
 *
 * @brief populations of hydrogenic energy levels by cascades (unfinished)
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: populations.cxx 375 2016-02-05 16:40:12Z iamp $
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

using namespace std;

//#define printcasc 1

const double pi = 3.1415926536;
const double clight=2.99792458e10;
const double melectron=510999.6;


double cascade(int counter, int ncascstep, int n1, int l1, int n_zero,
               int *n_list, int *l_list, double *pd_list, double *lt_list,
               const RADRATE& hydro, double *pcascdecay)
  {
  counter++;

  n_list[counter] = n1;
  l_list[counter] = l1;
  pd_list[counter]  = hydro.pdecay(n1,l1);
  lt_list[counter]  = hydro.life(n1,l1);

  static int i, k, n2max, docascades;
  static double prod, fcascade;
  
  double pdecay = 0.0;
  for (i=0; i<=counter; i++)
    {
    prod = 1.0;
    for (int k=0; k<=counter; k++) 
      { 
        if (k!=i) prod *= 1.0-lt_list[k]/lt_list[i];
      } 
    pdecay += pd_list[i]/prod;
    }
 
  *pcascdecay = pdecay; // fraction of flux going into cascades

  double fnl = 0.0;
  if (counter<ncascstep) { n2max = n1-1;} // full calc. for remaining casc.
             else { n2max = (n1<=n_zero) ? n1-1 : n_zero;}
  for(int n2=n2max; n2>0; n2--)
    {
    for(int l2=l1-1; l2<=l1+1; l2+=2)
      {
      if ( ( l2<0) || (l2>=n2) ) continue;
      double pdcasc = 0.0;
      fcascade = 0.0;
      if (counter<ncascstep)
         {
         fcascade = cascade(counter, ncascstep, n2, l2, n_zero, n_list, l_list,
                            pd_list, lt_list, hydro, &pdcasc);
         }
      fnl += hydro.branch(n1,l1,n2,l2)*((pdecay-pdcasc)+fcascade);
      } // end for (l2...)
    } // end (for (n2...)

#ifdef printcasc
  for(i=0;i<counter;i++)
    {
    printf("%3d %3d -> ",n_list[i],l_list[i]);
    }
  printf(" %3d %3d : %8.5f\n",n1,l1,fnl);
#endif

  return fnl;
  }

//////////////////////////////////////////////////////////////////////////

void Population(void)
{
  const double clight=2.99792458e10;
  const double melectron=510999.6;
  const double pi = 3.1415926536;
  double z,x1,x2,xd,ecool,rho,phi,mass;
  int nselective, ncut, nmax, softcut, ncascstep, selection; 
  char answer[2], filenameroot[200], filename[200], pfn[200];
  FILE *fin, *fout, *fout2, *fmatrix1, *fmatrix2, *fmatrix3;

  printf("\n Now calculating hydrogenic transition rates and");
  printf(" decay probabilities  \n");

  RADRATE hydro(nmax,z,vel,x1/gamma,x2/gamma);
 
  double* pd_list = new double[nmax];
  double* lt_list = new double[nmax];
  int* n_list = new int[nmax];
  int* l_list = new int[nmax];
  for(int i=0; i<nmax; i++)  {
    pd_list[i]=1.0; 
    lt_list[i]=1.0; 
    n_list[i]=1; 
    l_list[i]=0;}

  printf("\n\n Give number of cascade steps ( max if <0 ) ......: ");
  scanf("%d",&ncascstep);
  if (ncascstep<0) {ncascstep = nmax;}

  int read_previous = 0;
  printf("\n Multiply with values from a previous calculation?: ");
  scanf("%s",&answer);
  if ( (!strcmp(answer,"y")) || (!strcmp(answer,"Y")) )
    {
    read_previous = 1;

    printf("\n Give filename of previous calculation (*.fnl) ...: ");
    scanf("%s",&pfn);
    strcat(pfn,".fnl");

    char header[200];

    fin = fopen(pfn,"r");
    fgets(header,200,fin);
    printf("\n Header of previous calculation:\n %s \n",header);
 // overread next two lines
    fgets(header,200,fin);
    fgets(header,200,fin);
    }
      
  printf("\n Give filename for output (*.fnl, *.fm<i> i=1-3) .: ");
  scanf("%s",&filenameroot);

  strcpy(filename,filenameroot);
  strcat(filename,".fnl");
  fout = fopen(filename,"w");
  fprintf(fout,"z=%5.1f, ecool=%7.2f, x1=%7.2f cm, ",z,ecool,x1); 
  fprintf(fout,"x2=%7.2f cm, phi=%7.2f deg, ncascstep=%2d\n",
          x2,phi,ncascstep); 
  if (read_previous) {fprintf(fout," previous calculation: %s",pfn);}
  fprintf(fout,"\n   n    l            f  sigma(0 eV)      f*sigma\n");

  strcpy(filename,filenameroot);
  strcat(filename,".fn");
  fout2 = fopen(filename,"w");

  strcpy(filename,filenameroot);
  strcat(filename,".fm1");
  fmatrix1 = fopen(filename,"w");
  strcpy(filename,filenameroot);
  strcat(filename,".fm2");
  fmatrix2 = fopen(filename,"w");
  strcpy(filename,filenameroot);
  strcat(filename,".fm3");
  fmatrix3 = fopen(filename,"w");

  printf("\n Now calculating detection probabilities \n");

  // calculate detection probabilities

  int n1, l1;
     
  for(n1=1; n1<=nmax; n1++)
    {
    double fn = 0.0;
    int n12 = (n1-1)*n1/2;
    // find maximum RR cross section per n and store cross sections
    // for latter use:

    for(l1=0; l1<n1; l1++)
      {
      double tau = hydro.life(n1,l1);
      double pdecay = hydro.pdecay(n1,l1);
      double fnl = (1.0-pdecay);
      
      if (((n1==1)||(n1==2))&&(l1==0))
	{
	  fnl = 1.0; // 1s and 2s states are assumed to be always detected
	}
      else
	{
	  pd_list[0] = pdecay; lt_list[0]  = tau; n_list[0] = n1; l_list[0] = l1;
      
	  int n2max;
	  if (ncascstep) {n2max = n1-1;}        // full calculation for cascades
	  else { n2max = (n1<=n_zero) ? n1-1 : n_zero;}
	  for(int n2=n2max; n2>0; n2--)
	    {
	      int n22 = (n2-1)*n2/2;
	      for(int l2=l1-1; l2<=l1+1; l2+=2)
		{
		  if ( ( l2<0) || (l2>=n2) ) continue;
		  double pcascdecay = 0.0, fcascade = 0.0;
		  if (ncascstep)
		    {
		      fcascade = cascade(0, ncascstep, n2, l2, n_zero, n_list, l_list,
                                pd_list, lt_list, hydro, &pcascdecay);
		    }
		  double br = hydro.branch(n1,l1,n2,l2);
		  fnl += br*((pdecay-pcascdecay)+fcascade);
		} // end for (l2...)
	    } // end for (n2...)
	} // end else
      // read values from previous calculation

      int np,lp;
      double fnlp, sigmap, fsp;

      if (read_previous)
        {
        if (fscanf(fin,"%d %d %lf %lf %lf",&np,&lp,&fnlp,&sigmap,&fsp)==EOF) 
          {
          printf("\n !!! ATTENTION !!!");
          printf("\n Actual calculation not compatible with previous one!");
          printf("\n Program terminated at n = %4d.",n1-1);
          fclose(fin);
          fclose(fout);
          fclose(fmatrix1);
          fclose(fmatrix2);
          fclose(fmatrix3);
          break;
          }
        fnl*=fnlp;
        } // end if (read_previous)


      // output

      fprintf(fout,"%4d %4d %12.4G %12.4G %12.4G\n", 
                    n1,l1,fnl,sigma[l1],fnl*sigma[l1]/sigma_max);
      fprintf(fmatrix1,"%12.5g",fnl);
      fprintf(fmatrix2,"%12.5g",fnl*sigma[l1]/sigma_max);
      fprintf(fmatrix3,"%12.5g",pdecay);
      fn += (2.0*l1+1.0)*fnl;
      } // end for (l1...)

    fn /= n1*n1;
    fprintf(fout2,"%5d %12.4G\n",n1,fn);
  
    for (int l=n1; l<nmax; l++)
      {
      fprintf(fmatrix1,"%12.5g",0.0);
      fprintf(fmatrix2,"%12.5g",0.0);
      fprintf(fmatrix3,"%12.5g",0.0);
      }
    fprintf(fmatrix1,"\n");fprintf(fmatrix2,"\n");fprintf(fmatrix3,"\n");
     
    if ( (n1 % 10) == 0) { cout << "|"; } else { cout << "." ; } cout.flush();
                
    } // end for (n1...)
  if (read_previous) fclose(fin);
  fclose(fout);
  fclose(fout2);
  fclose(fmatrix1);
  fclose(fmatrix2);
  fclose(fmatrix3);
 
  printf("\n\n Another number of cascade steps? (y/n) : ");
  scanf("%s",&answer);
  if ( (!strcmp(answer,"y")) || (!strcmp(answer,"Y")) ) goto new_cascade;

  delete[] pd_list;
  delete[] lt_list;
  delete[] n_list;
  delete[] l_list;
  delete[] sigma;
  }

//////////////////////////////////////////////////////////////////////////

