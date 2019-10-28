/**
 * @file radrate.cxx
 *
 * @brief Class for holding and manipulating hydrogenic dipole transition rates 
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: radrate.cxx 375 2016-02-05 16:40:12Z iamp $
   @endverbatim
 *
 */
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include "osci.h"
#include "lifetime.h"
#include "radrate.h"

using namespace std;

RADRATE::~RADRATE()
  { 
  delete[] transrate_array; 
  delete[] branch_array; 
  delete[] lifetime_array; 
  delete[] pdecay_array;
  } 

RADRATE::RADRATE(int nmax, double z, double vel, double x1, double x2)
  {
  maxn = nmax; // stored internally
  int tdim = (nmax-1)*nmax*(nmax+1)+3;
  int ldim = (nmax+1)*nmax/2;
  fstream f1;

  // try to open file for reading
  f1.open("RADRATES",fstream::in);
  int readflag = f1.is_open();

  int nmaxf; 
  double zf;
  if (readflag) 
  { 
      f1 >> zf >> nmaxf ;
      if ( (!f1.good()) || (nmaxf < nmax) || (fabs(zf-z)>1.0E-6) )
      { // do not read from file 
	  f1.close(); 
	  readflag = 0; 
      }
  }

  int writeflag=1;
  if (!readflag)
  { // open file for writing
      f1.clear();
      f1.open("RADRATES",ios::out);
      if (f1.is_open()) 
      {
	  f1 << z << " " << nmax << endl;
          f1.setf(ios::scientific);
	  f1.precision(35);
      }
      else
      {
	  cout << "\nERROR:: file RADRATES not open for writing" << endl;
          writeflag = 0;
      }
  }
  else
    {
      cout << " reading previously calculated hydrogenic rates from file RADRATES" << endl;
    }
 

  transrate_array = new double[tdim];
  branch_array = new double[tdim];
  lifetime_array = new double[ldim];
  pdecay_array = new double[ldim];
  for (int i=0;i<ldim;i++) 
    {
    transrate_array[i] = 0.0;
    lifetime_array[i] = 1.0;
    pdecay_array[i] = 0.0;
    }
  for (int j=ldim; j<tdim; j++) transrate_array[j] = 0.0;

  transrate_array[0] = 0.0;  // 1s state
  lifetime_array[0] = 1.0;
  pdecay_array[0] = 0.0;
  cout << "." ; cout.flush();
  for (int n1=2; n1<=nmax; n1++)
  { 
      int n12 = (n1-1)*n1/2;
      for (int l1=0; l1<n1; l1++)
      {
	  int n2min = l1>0 ? l1 : 1;
	  double ht, sum = 0.0;
	  for (int n2=n2min; n2<n1; n2++)
	  {
	      int index = 1+(n1-2)*(n1-1)*n1+2*l1*(n1-1)+2*(n2-1);
	      if (l1+1<n2) 
	      {
		  if (readflag) { f1 >> ht; }
		  else
		  {
		      ht = hydrotrans(z,n1,l1,n2,l1+1); 
		      if (writeflag) f1 << ht <<endl; 
		  }
		  transrate_array[index]=ht; 
		  sum += ht;
	      }
	      if ((l1-1<n2)&&(l1>0))
	      { 
		  if (readflag) { f1 >> ht; }
		  else
		  {
		      ht = hydrotrans(z,n1,l1,n2,l1-1); 
		      if (writeflag) f1 << ht <<endl; 
		  }
		  transrate_array[index+1]=ht; 
		  sum += ht;
	      }
	  }
	  if (sum>0)
	  {
	      double tau = 1.0/sum;
	      lifetime_array[n12+l1] = tau;
	      for (int n2=n2min; n2<n1; n2++)
		{
		  int index = 1+(n1-2)*(n1-1)*n1+2*l1*(n1-1)+2*(n2-1);
		  branch_array[index]   =  tau*transrate_array[index]; 
		  branch_array[index+1] =  tau*transrate_array[index+1]; 
		}
	      double d = vel*tau;
	      if ( (x1-x2) > 1.e6)
	      { 
		  pdecay_array[n12+l1] = 1.0-d*(exp(-x1/d)-exp(-x2/d))/(x2-x1);
	      }
	      else
	      { 
		  pdecay_array[n12+l1] = 1-exp(-x1/d);
	      }
	  }
	  else
	  {
	      lifetime_array[n12+l1] = 1.0;
	      pdecay_array[n12+l1] = 0.0;
	  }
      } // end for (l1...)
      if ( (n1 % 10) == 0) { cout << "|";  }  else  { cout << "." ;  } 
      cout.flush();
  }// end for (n1...)
  f1.flush();
  f1.close();
}
             
void RADRATE::print() const
  {
  for (int n1=1; n1<=maxn; n1++)
    {
    int n12 = (n1-1)*n1/2;
    for (int l1=0; l1<n1; l1++)
      {
          cout.precision(10);
	  cout << n1 << " " << l1 << " " << lifetime_array[n12+l1] 
               << " " << pdecay_array[n12+l1] << endl;
      }
    }
  }

double RADRATE::trans(int n1, int l1, int n2, int l2) const
  { 
  if ( (n1<2) || (n1>maxn) ) return 0.0;
  l2 = (1+(l1-l2))/2; 
  if ( (l2<0) || (l2>1) ) return 0.0;
  int index = 1+(n1-2)*(n1-1)*n1+2*l1*(n1-1)+2*(n2-1);
  return transrate_array[index+l2];
  }

double RADRATE::branch(int n1, int l1, int n2, int l2) const
  { 
  if ( (n1<2) || (n1>maxn) ) return 0.0;
  l2 = (1+(l1-l2))/2; 
  if ( (l2<0) || (l2>1) ) return 0.0;
  int index = 1+(n1-2)*(n1-1)*n1+2*l1*(n1-1)+2*(n2-1);
  return branch_array[index+l2];
  }


double RADRATE::life(int n1, int l1) const
  { 
  return lifetime_array[(n1-1)*n1/2+l1]; 
  }

double RADRATE::pdecay(int n1, int l1) const
  { 
  return pdecay_array[(n1-1)*n1/2+l1]; 
  }

