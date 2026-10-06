// SPDX-License-Identifier: MIT
/**
 * @file populations.cxx
 *
 * @brief populations of hydrogenic energy levels by cascades (unfinished)
 *
 * @author Stefan Schippers
 *
 */
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include "radrate.h"
#include "sigmarr.h"
#include "hydromath.h"
#include "hydroconst.h"

using namespace std;

//#define printcasc 1

const double pi = hydroconst::pi;
const double clight=hydroconst::clight_cm_s;
const double melectron=hydroconst::mec2_eV;


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
    cout << setw(3) << n_list[i] << " " << setw(3) << l_list[i] << " -> ";
    }
  cout << " " << setw(3) << n1 << " " << setw(3) << l1 << " : " << setw(8) << setprecision(5) << fixed << fnl << "\n";
#endif

  return fnl;
  }

//////////////////////////////////////////////////////////////////////////

void Population(void)
{
  const double clight=hydroconst::clight_cm_s;
  const double melectron=hydroconst::mec2_eV;
  const double pi = hydroconst::pi;
  double z,x1,x2,xd,ecool,rho,phi,mass;
  int nselective, ncut, nmax, softcut, ncascstep, selection; 
  char answer;
  string filenameroot, filename, pfn;
  ifstream fin;
  ofstream fout, fout2, fmatrix1, fmatrix2, fmatrix3;

  cout << "\n Now calculating hydrogenic transition rates and";
  cout << " decay probabilities  \n";

  RADRATE hydro(nmax,z,vel,x1/gamma,x2/gamma);
  
  std::vector<double> pd_list(nmax);
  std::vector<double> lt_list(nmax);
  std::vector<int> n_list(nmax);
  std::vector<int> l_list(nmax);
  for(int i=0; i<nmax; i++)  {
    pd_list[i]=1.0; 
    lt_list[i]=1.0; 
    n_list[i]=1; 
    l_list[i]=0;}

  cout << "\n\n Give number of cascade steps ( max if <0 ) ......: ";
  cin >> ncascstep;
  if (ncascstep<0) {ncascstep = nmax;}

  int read_previous = 0;
  cout << "\n Multiply with values from a previous calculation?: ";
  cin >> answer;
  if ( (answer=='y') || (answer=='Y') )
    {
    read_previous = 1;

    cout << "\n Give filename of previous calculation (*.fnl) ...: ";
    cin >> pfn;
    pfn += ".fnl";

    string header;

    fin.open(pfn);
    getline(fin,header);
    cout << "\n Header of previous calculation:\n " << header << "\n";
 // overread next two lines
    getline(fin,header);
    getline(fin,header);
    }
      
  cout << "\n Give filename for output (*.fnl, *.fm<i> i=1-3) .: ";
  cin >> filenameroot;

  filename = filenameroot;
  filename += ".fnl";
  fout.open(filename);
  fout << "z=" << fixed << setprecision(1) << z << ", ecool=" << setw(7) << setprecision(2) << ecool << ", x1=" << setw(7) << setprecision(2) << x1 << " cm, "; 
  fout << "x2=" << setw(7) << setprecision(2) << x2 << " cm, phi=" << setw(7) << setprecision(2) << phi << " deg, ncascstep=" << setw(2) << ncascstep << "\n"; 
  if (read_previous) {fout << " previous calculation: " << pfn;}
  fout << "\n   n    l            f  sigma(0 eV)      f*sigma\n";

  filename = filenameroot;
  filename += ".fn";
  fout2.open(filename);

  filename = filenameroot;
  filename += ".fm1";
  fmatrix1.open(filename);
  filename = filenameroot;
  filename += ".fm2";
  fmatrix2.open(filename);
  filename = filenameroot;
  filename += ".fm3";
  fmatrix3.open(filename);

  cout << "\n Now calculating detection probabilities \n";

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
		      fcascade = cascade(0, ncascstep, n2, l2, n_zero, n_list.data(), l_list.data(),
                                pd_list.data(), lt_list.data(), hydro, &pcascdecay);
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
        if (!(fin >> np >> lp >> fnlp >> sigmap >> fsp))
          {
          cout << "\n !!! ATTENTION !!!";
          cout << "\n Actual calculation not compatible with previous one!";
          cout << "\n Program terminated at n = " << setw(4) << n1-1 << ".";
          fin.close();
          fout.close();
          fmatrix1.close();
          fmatrix2.close();
          fmatrix3.close();
          break;
          }
        fnl*=fnlp;
        } // end if (read_previous)


      // output

      fout << setw(4) << n1 << " " << setw(4) << l1 << " " << setw(12) << uppercase << defaultfloat << setprecision(4) << fnl << " " << setw(12) << uppercase << defaultfloat << setprecision(4) << sigma[l1] << " " << setw(12) << uppercase << defaultfloat << setprecision(4) << fnl*sigma[l1]/sigma_max << "\n";
      fmatrix1 << setw(12) << setprecision(5) << fnl;
      fmatrix2 << setw(12) << setprecision(5) << fnl*sigma[l1]/sigma_max;
      fmatrix3 << setw(12) << setprecision(5) << pdecay;
      fn += (2.0*l1+1.0)*fnl;
      } // end for (l1...)

    fn /= n1*n1;
    fout2 << setw(5) << n1 << " " << setw(12) << uppercase << defaultfloat << setprecision(4) << fn << "\n";
  
    for (int l=n1; l<nmax; l++)
      {
      fmatrix1 << setw(12) << setprecision(5) << 0.0;
      fmatrix2 << setw(12) << setprecision(5) << 0.0;
      fmatrix3 << setw(12) << setprecision(5) << 0.0;
      }
    fmatrix1 << "\n"; fmatrix2 << "\n"; fmatrix3 << "\n";
     
    if ( (n1 % 10) == 0) { cout << "|"; } else { cout << "." ; } cout << flush;
                
    } // end for (n1...)
  if (read_previous) fin.close();
  fout.close();
  fout2.close();
  fmatrix1.close();
  fmatrix2.close();
  fmatrix3.close();
 
  cout << "\n\n Another number of cascade steps? (y/n) : ";
  cin >> answer;
  if ( (answer=='y') || (answer=='Y') ) goto new_cascade;


  }

//////////////////////////////////////////////////////////////////////////

