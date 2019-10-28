/**
 * @file readxsec.cxx
 *
 * @brief Reading theoretical cross sections from files 
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: readxsec.cxx 653 2019-10-27 13:59:15Z iamp $
   @endverbatim
 *
 */
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include "fraction.h"
#include "peakfunctions.h"

using namespace std;

double efac = 1.0, wfac = 1.0, sfac = 1.0;
int nmax = 0;

///////////////////////////////////////////////////////////////////////////////

void read_vector(string line, int linewidth, double *v, int pos0, int no)
  {  
  string zahl;
  for (int k=0; k<no; k++)
    {
    zahl=line.substr(k*12,12);
    istringstream buf(zahl,istringstream::in);
    buf>>v[pos0+k];
    }
  }

///////////////////////////////////////////////////////////////////////////////

int open_peak(char *filename)
  {
  int npts=0;
  char c,line[200];  
  ifstream fin(filename); 
  if (!fin) {return -1;}
  
  while (!fin.eof())
    {
    fin.get(line,199,'\n');
    fin.get(c);
    npts++;
    }
  fin.close();
  return npts-1;
  } 

///////////////////////////////////////////////////////////////////////////////

void read_peak(int mode, char *filename, double *energy, double *strength, double *qFano,
		double *wLorentz, double *wGauss, int npts, double &emin, double &emax)
  {
  ifstream fin(filename); 
  
  emin =  1e99; emax = -1e99;
  for(int i=0; i<npts; i++)
    {
      qFano[i] = 0.0;
      wGauss[i] = 0.0;
      if (mode == PEAK_Lorentz) 
	{
	  fin >> energy[i] >> strength[i]  >> wLorentz[i];
	}
      else if (mode == PEAK_Lorentz_Steih) 
	{
	  fin >> energy[i] >> wLorentz[i] >> strength[i];
	}
      else if (mode == PEAK_Fano) 
	{// Fano peaks
	  fin >> energy[i] >> strength[i]  >> qFano[i] >> wLorentz[i];
	}
      else if (mode == PEAK_Voigt) 
	{
	  fin >> energy[i] >> strength[i]  >> wLorentz[i] >> wGauss[i];
	}
      else if (mode == PEAK_FanoVoigt) 
	{// FanoVoigt peaks
	  fin >> energy[i] >> strength[i]  >> qFano[i] >> wLorentz[i] >> wGauss[i];
	}
      else
	{ //PEAK_delta
	  fin >> energy[i] >> strength[i];
	}
      energy[i] *= efac;
      if (energy[i]<(emin)) emin = energy[i];
      if (energy[i]>(emax)) emax = energy[i];
      wLorentz[i] *= wfac;
      strength[i] *= sfac;
    }
  fin.close();
  cout.flush();
  }

///////////////////////////////////////////////////////////////////////////////

int open_Griffin_autostructure(char *filename)
  {
  int npts=0, incr = 0, i, n, l;
  char line[200], teststr1[6], teststr2[4];
  istringstream buf(line,istringstream::in);  

  ifstream fin(filename); 
  if (!fin) {return -1;}

  nmax = 0; // nmax is a global variable in this module
  while (!fin.eof())
    {
    buf.clear();
    fin.getline(line,200); 
    for (i=0;i<5;i++) { teststr1[i]=line[i+1]; } teststr1[5] = '\0';
    for (i=0;i<3;i++) { teststr2[i]=line[i+5]; } teststr2[3] = '\0';
 
    if (!strcmp(teststr1,"total"))
       {
       fin.getline(line,200);
       fin.getline(line,200);
       incr = 0;
       }
    else if (!strcmp(teststr2,"n=:"))
       {
       buf.seekg(0); buf.ignore(8); buf >> n; buf.ignore(6); buf >> l;
       if (n>nmax) nmax=n; // nmax is a global variable in this module
       fin.getline(line,200);
       incr = 1;
       }
    else if (incr) 
       {
       npts++;
       }
    }
  fin.close();

  return npts;
  } 

///////////////////////////////////////////////////////////////////////////////

void read_Griffin_autostructure(char *filename, double *energy,
                 double *strength, double *width, double &emin, double &emax)
  {
  int i, n, l, fdim = nmax*(nmax+1)/2, use_fractions=0;
  double *fraction = new double[fdim];
  for (i=0;i<fdim;i++) {fraction[i] = 1.0;}
  char header[200], answer;
  cout << "\n Read surviving fractions from file? (y/n) : ";
  cin >> answer;
  if ( (answer=='y') || (answer=='Y') ) 
     {
     use_fractions = readfraction(fraction, header, nmax);
     cout << " surviving fractions " << header << '\n';
     cout.flush();
     cout << "\n Perform l-averaging of fractions? (y/n) ..: ";
     cin >> answer; 
     if ( (answer=='y') || (answer=='Y') ) 
       {
       for(n=1;n<=nmax;n++)
         {
         double sum=0.0, n2=n*n;
         for(l=0;l<n;l++)
           {
           sum += (2.0*l+1.0)*fraction[n*(n-1)/2+l];
           }
         double laverage = sum/n2;
         for(l=0;l<=n;l++)
           {
           fraction[n*(n-1)/2+l] = laverage;
           }
         cout << " n = " << n << " : " << laverage << "\n"; cout.flush();
         }
       }
     }  

  int npts=-1, incr = 0, badcount=0;
  char line[200], teststr1[6], teststr2[4];
  istringstream buf(line,istringstream::in);  
  emin = 1E99; emax=-1E99;
  ifstream fin(filename); 

  while (!fin.eof())
    {
    buf.clear();
    fin.getline(line,200); 
    for (i=0;i<5;i++) { teststr1[i]=line[i+1]; } teststr1[5] = '\0';
    for (i=0;i<3;i++) { teststr2[i]=line[i+5]; } teststr2[3] = '\0';
 
    if (!strcmp(teststr1,"total"))
       {
       fin.getline(line,200);
       fin.getline(line,200);
       buf.clear();
       incr = 0;
       }
    else if (!strcmp(teststr2,"n=:"))
       {
       buf.seekg(0); buf.ignore(8); buf >> n; buf.ignore(6); buf >> l;
       fin.getline(line,200);
       buf.clear();
       incr = 1;
       }
    else if (incr) 
       {
       double e, s;
       buf.seekg(0);
       if (buf >> e >> s)
          {
          npts++;
          e *= efac;
          s *= sfac;
          if (e>(emax)) emax=e;
          if (e<(emin)) emin=e;
          energy[npts] = e;
          strength[npts] = s*fraction[n*(n-1)/2+l];
          width[npts] = 0.0;
          }
       else
          {
          badcount++;
          }
       }
    }
  if (badcount) 
    {
    cout << "\n " << badcount << " reading errors\n";
    cout.flush();
    }
  fin.close();
  if (use_fractions) delete[] fraction;
  }  

///////////////////////////////////////////////////////////////////////////////

int open_Pindzola(char *filename)
  {
  int npts;  
  double dummy;
  ifstream fin(filename); 
  if (!fin) {return -1;}
  fin >> dummy >> dummy >> dummy >> dummy;
  fin >> dummy >> dummy;
  fin >> dummy >> dummy;
  fin >> npts >> dummy >> dummy >> dummy;
  fin.close();

  return npts;
  } 

///////////////////////////////////////////////////////////////////////////////

void read_Pindzola(char *filename, double *energy, double *sigma, 
                   double *width, double &emin, double &emax)
  {

  int npts;  
  double dummy, field;
  ifstream fin(filename); 
  fin >> dummy >> dummy >> dummy >> dummy;
  fin >> dummy >> dummy;
  fin >> dummy >> dummy;
  fin >> npts >> field >> dummy >> dummy;

  emin =  1e99; emax = -1e99;
  int i;
  for (i=0; i<npts; i++) 
    { 
    fin >> energy[i]; 
    energy[i] *= efac; 
    if (energy[i]<emin) emin = energy[i];
    if (energy[i]>emax) emax = energy[i];
    }
  for (i=0; i<npts; i++) 
    { 
    fin >> sigma[i]; 
    sigma[i] *= sfac*(energy[1]-energy[0]); 
    width[i] = 0.0; 
    }
  fin.close();

  char header[200], answer[2];
  cout << " Read surviving fractions from file? (y/n) : ";
  cin >> answer;
  if ( (!strcmp(answer,"y")) || (!strcmp(answer,"Y")) ) 
     {
     int nmax, n, l;
     double slim,z,elo,ehi,fn;

     cout << "\n Give effective ionic charge ..............: ";
     cin >> z; 
     cout << "\n Give series limit in eV ..................: ";
     cin >> slim;
     cout << "\n Give maximum n ...........................: ";
     cin >> nmax; n=nmax+1;
     int fdim = nmax*(nmax+1)/2;
     double *fraction = new double[fdim];
     for (i=0;i<fdim;i++) {fraction[i] = 1.0;}
     nmax = readfraction(fraction, header, nmax);
     cout << " surviving fractions " << header << '\n';
     cout.flush();

     elo = slim-13.606*z*z/(n-0.5)/(n-0.5);
     ehi = slim;
     fn  = 0.0;
     for (i=npts-1;i>=0;i--)
       {
       if (energy[i]<elo)
         { 
         n--; 
         elo = slim-13.606*z*z/(n-0.5)/(n-0.5);
         ehi = slim-13.606*z*z/(n+0.5)/(n+0.5);
         fn = 0.0;
         for (l=0; l<n; l++) { fn += (2*l+1)*fraction[n*(n-1)/2+l]; }
         fn /= (n*n);
         cout << "n=" << n << ": fraction = " << fn <<"\n";
         cout.flush();
         }
       if (n>nmax) {fn = 0.0;} // else {fn = 1.0;}
       if (fn==1.0) break; 
       if (sigma[i]>0)
          {
          sigma[i] *= fn;
          cout << n << " : " << elo << " " << energy[i] << " " << ehi << "\n"; 
          }
       } // end for(i...)
     delete[] fraction;
     } // end if

  cout << "\n Calculate integrated rate coeff? (y/n) ...: ";
  cin >> answer;
  if ( (!strcmp(answer,"y")) || (!strcmp(answer,"Y")) ) 
    {
    double elo, ehi, strength = 0.0;
    cout << "\n Give boundaries of integration range in eV: ";
    cin >> elo >> ehi; 
    for(i=0;i<npts;i++)
      {
      if ((energy[i]>=elo) && (energy[i]<=ehi)) 
         strength += sigma[i]*2.998e10*sqrt(energy[i]/511000);
      }
    cout << "\n" << field  << " V/cm : " << strength << " eVcm^2\n";       
    ofstream fout("integrals.dat",ios::app);
    fout << elo << " " << ehi << " " << field  << " " << strength << "\n";
    fout.close();
    }  
  cout.flush();
  }

///////////////////////////////////////////////////////////////////////////////

int open_Cowan(char *filename)
  {
  int npts;  
  ifstream fin(filename); 
  if (!fin) {return -1;}
  fin >> npts;
  fin.close();
  return npts;
  } 

///////////////////////////////////////////////////////////////////////////////

void read_Cowan(char *filename, double *energy, double *sigma, 
                double *width, double &emin, double &emax)

  {
  int i, n, l, npts;
  char line[200];
  istringstream buf(line,istringstream::in);

  ifstream fin(filename); 
  fin.getline(line,200);
  buf >> npts;
  buf.seekg(0);

  int *n_store = new int[npts];
  int *l_store = new int[npts];

  nmax = 0; // nmax is a global variable in this module
  for (i=0; i<npts; i++)
    {
    fin.getline(line,200);
    buf >> n_store[i] >> l_store[i] >> energy[i] >> sigma[i] >> width[i];
    buf.seekg(0);
//    cout <<  n_store[i] << ' ' << l_store[i] << ' ' << energy[i] << ' ';
//    cout <<  sigma[i] << ' ' << width[i] << '\n';
    if (n_store[i]>nmax) { nmax = n_store[i]; }
    }
  fin.close();

  int fdim = nmax*(nmax+1)/2;
  double *fraction = new double[fdim];
  for (i=0;i<fdim;i++) {fraction[i] = 1.0;}
  char header[200], answer[2];
  cout << " Read surviving fractions from file? (y/n) : ";
  cin >> answer;
  if ( (!strcmp(answer,"y")) || (!strcmp(answer,"Y")) ) 
     {
     nmax = readfraction(fraction, header, nmax);
     cout << " surviving fractions " << header << '\n';
     cout.flush();
     }

  double *psum = new double[nmax+1];
  double *wsum = new double[nmax+1];
  for (i=0; i<=nmax; i++) {psum[i] = 0.0; wsum[i] = 0.0;}
  emin =  1e99; emax = -1e99;
  for (i=0; i<npts; i++)
    {
    n = n_store[i]; l = l_store[i];
    energy[i] *= efac;
    if (energy[i]<emin) emin = energy[i];
    if (energy[i]>emax) emax = energy[i];
    sigma[i] *= sfac;
    wsum[n] += sigma[i];
    sigma[i] *= fraction[n*(n-1)/2+l];
    psum[n] += sigma[i];
    }

  char fn[200];  
  strcpy(fn,filename);
  strcat(fn,"_n");  
  ofstream fout(fn); 
  for (n=1;n<=nmax;n++)
     {
     double pmean = wsum[n]>0 ? psum[n]/wsum[n] : 0.0;
     fout << n << ' ' << pmean << '\n';
     }
  fout.close();
  cout << " n-specific detection probabilities written to file " << fn << '\n';
  cout.flush();
  delete[] wsum;
  delete[] psum;
  delete[] l_store;
  delete[] n_store;
  delete[] fraction; 
  }

///////////////////////////////////////////////////////////////////////////////

int open_autostructure(char *filename, int &number_of_levels)
  {
  int npts;  
  ifstream fin(filename); 
  if (!fin) {return -1;}
  fin >> npts;

  int line_count = 0;
  char line[200]; 
  while (!fin.eof())
    {
      fin >> line;
      line_count++;
    }
  fin.close();
  number_of_levels = 6*line_count/npts-1;
  return npts-1;
  } 

///////////////////////////////////////////////////////////////////////////////

void read_autostructure(char *filename, int level_number, double *energy, double *strength, 
                        double *width, double &emin, double &emax)
  {
  int npts;  
  ifstream fin(filename); 
  fin >> npts;
  char line[200]; 
  int i,j,k,pos0 = 0, nospec=1;
  int line_count = 0;

  for (i=0; i<npts/6; i++)
    {
    fin >> line;
    line_count++;
    read_vector(line,200,energy,pos0,6);
    pos0 += 6;
    }
  if (npts % 6 > 0) 
    {
    fin >> line;
    read_vector(line,200,energy,pos0,(npts % 6));
    if (npts % 6 > 1) line_count++;
    }
  npts = npts-1; // energy bins are given by min and max value, 
  // from here on we use the center of the bin as energy value 
  // therefore the number of bins in decreased by one 

  // overread information on lower levels that is not requested
  // 'level number' is a global variable set in open_autostructure()
  for (j=1; j < level_number; j++)
    {
      for (k=0; k<line_count; k++) fin >> line;
    }


  pos0 = 0;
  for (i=0; i<npts/6; i++)
    {
    fin >> line;
    read_vector(line,200,strength,pos0,6);
    pos0 += 6;
    }
  if (npts % 6 > 0) 
    {
    fin >> line;
    read_vector(line,200,strength,pos0,npts % 6);
    }
  fin.close();

  //strcat(filename,".col");
  //ofstream fout(filename);
  for (i=0; i<npts; i++)
    {
      width[i] = 0.0;
      double binwidth = (energy[i+1]-energy[i])*efac; // not that we have one more energy point than cross section point
      energy[i] = energy[i]*efac+0.5*binwidth; // use the center of a bin as energy
      strength[i] *=sfac*binwidth;
      //fout << energy[i] << " " << strength[i] << endl;
    }
  //fout.close();
  emin =  energy[0];
  emax =  energy[npts-1];
}

///////////////////////////////////////////////////////////////////////////////

int open_autosrr(char *filename)
  {
  int npts;  
  ifstream fin(filename); 
  if (!fin) {return -1;}
  fin >> npts;
  fin.close();
  return npts;
  } 

///////////////////////////////////////////////////////////////////////////////

void read_autosrr(char *filename, double *energy, double *strength, 
                        double *width, double emin, double emax)
  {
  int npts;  
  ifstream fin(filename); 
  fin >> npts;
  for (int j=0; j<npts; j++)
    {
      fin >> energy[j] >> strength[j];
      energy[j] *= efac;
      strength[j] *= sfac;
      width[j] = 0.0;
      cout << energy[j] << " " << strength[j] << endl;
    }
  fin.close();

  emin =  energy[0];
  emax =  energy[npts-1];
}

///////////////////////////////////////////////////////////////////////////////

int open_Badnell_nl(char *filename)
  {
  int npts=0, n, l;
  double e, s;
  char line[200];
  istringstream buf(line,istringstream::in);

  ifstream fin(filename); 
  if (!fin) {return -1;}
  for (int i=1;i<=3;i++)
    {
    fin.getline(line,200);
    buf.seekg(0);
    }
  while (!fin.eof())
    {
    fin.getline(line,200);
    buf >> n >> l >> e >> s;
    if (buf.fail()) 
      {
      buf.clear();
      }
    else
      {
      npts++;
      }
    buf.seekg(0);
    }
  fin.close();
  return npts;
  } 

///////////////////////////////////////////////////////////////////////////////

void read_Badnell_nl(char *filename, double *energy, double *sigma, 
                     double *width, int npts, double emin, double emax)

  {
  int i=0, n=0, l=0, nn=0, ll=0;
  double e, s, binwidth;
  int *n_store = new int[npts];  // storage for n and l
  int *l_store = new int[npts];  // of energy bins
  char line[200];
  istringstream buf(line,istringstream::in);
  
  ifstream fin(filename); 

  fin.getline(line,200);
  buf.seekg(11);
  buf >> binwidth;
  fin.getline(line,200);
  fin.getline(line,200);
  buf.seekg(0);
  i = 0;
  nmax = 0; // namx is a global variable in this module
  while (!fin.eof())
    {
    fin.getline(line,200);
    buf >> nn >> ll >> e >> s; // read 4 numbers
    if (buf.fail())   // in case of error only 2 numbers are there
      {               // which are to be interpreted as n and l
      buf.clear();   // clear error condition
      n = nn;
      l = ll;
      if (n>nmax) { nmax=n; }
      }
    else
      {
      energy[i] = e;
      sigma[i] = s;
      n_store[i] = n;
      l_store[i] = l;
      i++;
      if (i==npts) break;
      }
    buf.seekg(0);
    }
  fin.close();

  int fdim = nmax*(nmax+1)/2;
  double *fraction = new double[fdim];
  for (i=0;i<fdim;i++) {fraction[i] = 1.0;}
  char header[200], answer[2];
  cout << " Read surviving fractions from file? (y/n) : ";
  cin >> answer;
  if ( (!strcmp(answer,"y")) || (!strcmp(answer,"Y")) ) 
     {
     nmax = readfraction(fraction, header, nmax);
     cout << " surviving fractions " << header << '\n';
     cout.flush();
     }

  double *psum = new double[nmax+1];
  double *wsum = new double[nmax+1];
  for (i=0; i<=nmax; i++) {psum[i] = 0.0; wsum[i] = 0.0;}
  emin =  1e99; emax = -1e99;
  for (i=0; i<npts; i++)
    {
    n = n_store[i]; l = l_store[i];
//    if (n<10) 
//      {cout << n << ' ' << l << ' ' << fraction[(n-1)*(n-1)+l] << '\n';}
    energy[i] *= efac;
    if (energy[i]<emin) emin = energy[i];
    if (energy[i]>emax) emax = energy[i];
    sigma[i] *= sfac*binwidth;
    wsum[n] += sigma[i];
    sigma[i] *= fraction[(n-1)*n/2+l];
    psum[n] += sigma[i];
    width[i] = 0.0;
    }

  char fn[200];  
  strcpy(fn,filename);
  strcat(fn,"_n");  
  ofstream fout(fn); 
  for (n=1;n<=nmax;n++)
     {
     double pmean = wsum[n]>0 ? psum[n]/wsum[n] : 0.0;
     fout << n << ' ' << pmean << '\n';
     }
  fout.close();
  cout << " n-specific detection probabilities written to file " << fn << '\n';
  cout.flush();
  delete[] wsum;
  delete[] psum;
  delete[] l_store;
  delete[] n_store;
  delete[] fraction; 
  }

///////////////////////////////////////////////////////////////////////////////

int open_Griffin_n(char *filename)
  {
  int n, npts=0;
  double d1,d2,d3,d4,d5,d6;
  
  char line[200];  
  istringstream buf(line,istringstream::in);
  ifstream fin(filename); 
  if (!fin) {return -1;}
  
  fin.getline(line,200);
  fin.getline(line,200);
  fin.getline(line,200);
  fin.getline(line,200);
  fin.getline(line,200);
  while (!fin.eof())
    {
    buf.clear();
    fin.getline(line,200);
    buf.seekg(0);
    if (buf >> n >> d1 >> d2 >> d3 >> d4 >> d5 >> d6) 
      {
      npts++;
//    cout << n << "\n";cout.flush();
      }
    }
  fin.close();

  nmax = n; // nmax is a global variable within this module
  return 2*npts;
  } 

///////////////////////////////////////////////////////////////////////////////

void read_Griffin_n(char *filename, double *energy,
                 double *strength, double *width, double emin, double emax)
  {
  int i, n, l, fdim = nmax*(nmax+1)/2, use_fractions = 0;
  double *fraction = new double[fdim];
  for (i=0;i<fdim;i++) {fraction[i] = 1.0;}
  char header[200], answer;
  cout << "\n Read surviving fractions from file? (y/n) : ";
  cin >> answer;
  if ( (answer=='y') || (answer=='Y') ) 
     {
     use_fractions = readfraction(fraction, header, nmax);
     cout << " surviving fractions " << header << '\n';
     cout.flush();

     for(n=1;n<=nmax;n++)
       {
       double sum=0.0, n2 =n*n;
       for(l=0;l<n;l++)
         {
         sum += (2.0*l+1.0)*fraction[n*(n-1)/2+l];
         }
       double laverage = sum/n2;
       for(l=0;l<n;l++)
         {
         fraction[n*(n-1)/2+l] = laverage;
         }
       cout << " n = " << n << " : " << laverage << '\n'; cout.flush(); 
       }
     }  

  char line[200];
  istringstream buf(line,istringstream::in);

  emin = 1E99; emax=-1E99;

  double deltae,e1,s1,e2,s2,et,st;
  ifstream fin(filename); 

  fin.getline(line,200);
  fin.getline(line,200);
  buf.clear();
  fin.getline(line,200);
  buf.seekg(47); buf >> deltae;
  fin.getline(line,200);
  fin.getline(line,200);

  i=0;
  while((!fin.eof()) && (2*i <nmax))
    {
    buf.clear();
    fin.getline(line,200);
    buf.seekg(0);
    if (buf >> n >> e1 >> s1 >> e2 >> s2 >> et >> st)
       {
       e1 *= efac;
       s1 *= sfac*deltae;
       e2 *= efac;
       s2 *= sfac*deltae;
       if (e1>(emax)) emax=e1;
       if (e1<(emin)) emin=e1;
       if (e2>(emax)) emax=e2;
       if (e2<(emin)) emin=e2;
       energy[2*i] = e1;
       strength[2*i] = s1*fraction[n*(n-1)/2];
       width[2*i] = 0.0;
       energy[2*i+1] = e2;
       strength[2*i+1] = s2*fraction[n*(n-1)/2];
       width[2*i+1] = 0.0;
       //cout <<energy[2*i]<<" "<<strength[2*i]<<" "<< width[2*i]<<"\n";
       //cout <<energy[2*i+1]<<" "<<strength[2*i+1]<<" "<< width[2*i+1]<<"\n";
       }
    i++;
    }
  if (use_fractions) delete[] fraction;
  cout.flush();
  fin.close();
  }  

///////////////////////////////////////////////////////////////////////////////

int open_Griffin_nl(char *filename)
  {
  int n, l, npts=0;
  double d1,d2;
  
  char line[200];  
  istringstream buf(line,istringstream::in);
  ifstream fin(filename); 
  if (!fin) {return -1;}
  
  fin.getline(line,200);
  fin.getline(line,200);
  while (!fin.eof())
    {
    buf.clear();
    fin.getline(line,200);
    buf.seekg(0);
    if (buf >> n >> l >> d1 >> d2)
      {
      npts++;
//    cout << n << "\n";cout.flush();
      }
    }
  fin.close();

  nmax = n; // nmax is a global variable within this module
  return npts;
  } 

///////////////////////////////////////////////////////////////////////////////

void read_Griffin_nl(char *filename, double *energy,
                 double *strength, double *width, double emin, double emax)
  {
  int i, n, l, fdim = nmax*(nmax+1)/2, use_fractions = 0;
  double *fraction = new double[fdim];
  for (i=0;i<fdim;i++) {fraction[i] = 1.0;}
  char header[200], answer;
  cout << "\n Read surviving fractions from file? (y/n) : ";
  cin >> answer;
  if ( (answer=='y') || (answer=='Y') ) 
     {
     use_fractions = readfraction(fraction, header, nmax);
     cout << " surviving fractions " << header << '\n';
     cout.flush();
     cout << "\n Perform l-averaging of fractions? (y/n) ..: ";
     cin >> answer; 
     if ( (answer=='y') || (answer=='Y') ) 
       { 
       for(n=1;n<=nmax;n++)
         {
         double sum=0.0, n2 =n*n;
         for(l=0;l<n;l++)
           {
           sum += (2.0*l+1.0)*fraction[n*(n-1)/2+l];
           }
         double laverage = sum/n2;
         for(l=0;l<n;l++)
           {
           fraction[n*(n-1)/2+l] = laverage;
           }
         cout << " n = " << n << " : " << laverage << '\n'; cout.flush(); 
         }
       }
     }  

  char line[200];
  istringstream buf(line,istringstream::in);

  emin = 1E99; emax=-1E99;

  double e,s;
  ifstream fin(filename); 

  fin.getline(line,200);
  fin.getline(line,200);

  i=0; n=0; l=0;
  while( (!fin.eof()) && (n<=nmax) && (l<nmax-1))
    {
    buf.clear();
    fin.getline(line,200);
    buf.seekg(0);
    if (buf >> n >> l >> e >> s)
       {
       e *= efac;
       s *= sfac;
       if (e>(emax)) emax=e;
       if (e<(emin)) emin=e;
       energy[i] = e;
       strength[i] = s*fraction[n*(n-1)/2+l];
       width[i] = 0.0;
     //cout << n << " " << l << " " << energy[i] << " " << strength[i] << "\n";
       }
    i++;
    }
  if (use_fractions) delete[] fraction;
  cout.flush();
  fin.close();
  }  


///////////////////////////////////////////////////////////////////////////////

int open_Fritzsche(char *filename)
  {
  int npts;  
  char line[200];
  ifstream fin(filename); 
  if (!fin) {return -1;}
  npts=0;
  while (!fin.eof())
    {
    fin.getline(line,200);
    npts++;
    }
  npts--;
  fin.close();
  return npts;
  } 

///////////////////////////////////////////////////////////////////////////////

void read_Fritzsche(char *filename, double *energy, double *sigma, 
		    double *width, int npts, double &emin, double &emax)

  {
  int i, Babuskin=0, widthfactor=1;
  ifstream fin(filename); 

  double *Aami  = new double[npts];
  double *Aa    = new double[npts];
  double *ArBab  = new double[npts];
  double *ArCoul = new double[npts];
  double *SBab  = new double[npts];
  double *SCoul  = new double[npts];

  cout << "\n Coulomb or Babuskin gauge? (0=Coul/1=Bab/2=negCoul/3=negBab) : ";
  cin >> Babuskin;
  if (Babuskin==2) {
    widthfactor=-1;
    Babuskin=0;
    // cout<<endl<<"Widht are multiplied by -1"<<endl;
  } else if (Babuskin==3) {
    widthfactor=-1;
    Babuskin=1;
    // endl cout<<"Widht are multiplied by -1"<<endl;
  } 

  emin = 1E99; emax=-1E99;
  for (i=0; i<npts; i++) {
    fin.ignore(40);
    fin >> energy[i] >> Aami[i] >> Aa[i] >> ArBab[i] >> ArCoul[i] >> SBab[i] >> SCoul[i];
    // cout <<  energy[i] << ' ' << Aami[i] << ' ' << Aa[i] << ' ';
    // cout <<  ArBab[i] << ' ' << ArCoul[i] << ' ' <<  SBab[i] << ' ' << SCoul[i] << '\n';
    if (Babuskin) { 
        // cout<<"Babuskin gauge used"<<endl;
	width[i] = widthfactor*6.58212e-16*(Aa[i]+ArBab[i]);
	sigma[i] = SBab[i]*sfac;
    } else {
        // cout<<"Coulomb gauge used"<<endl;
	width[i] = widthfactor*6.58212e-16*(Aa[i]+ArCoul[i]);
	sigma[i] = SCoul[i]*sfac;
    }
    energy[i] *= efac;
    // cout<<energy[i]<<" "<<width[i]<<" "<<sigma[i]<<endl;
    if (energy[i]>emax) emax=energy[i];
    if (energy[i]<emin) emin=energy[i];
  }
  fin.close();

  delete[] SCoul;
  delete[] SBab;
  delete[] ArCoul;
  delete[] ArBab;
  delete[] Aa; 
  delete[] Aami; 
  }

///////////////////////////////////////////////////////////////////////////////

int open_Grasp2K(char *filename)
  {
  int npts;  
  char line[200];
  ifstream fin(filename); 
  if (!fin) {return -1;}
  npts=0;
  while (!fin.eof())
    {
    fin.getline(line,200);
    npts++;
    }
  npts--;
  fin.close();
  return npts;
  } 

///////////////////////////////////////////////////////////////////////////////
void read_Grasp2K(char *filename, double *energy, double *strength, int npts, double &emin, double &emax)
{
  const double invcm2eV = 1.239842e-4; // conversion from 1/cm to eV
  const double f2MbeV= 109.761;     // conversion from oscillator strength to integrated cross section
  const string Lsymb="SPDFGHIKLMNOPQRSTUVWXYZ";

  double wavenumber,wavelength,gf,A,S,dT;

  int table_format=0;
  cout << "\n Which table format has been used in Grasp2K rtransitiontable?";
  cout << "\n 1: Lower & Upper & Energy diff. & wavelength & S  & gf & A  & dT"; 
  cout << "\n 2: Lower & Upper & Energy diff. & wavelength & gf & A  & dT"; 
  cout << "\n 3: Lower & Upper & Energy diff. & wavelength & gf & A"; 
  cout << "\n 4: Lower & Upper & Energy diff. & S          & gf & A  & dT"; 
  cout << "\n 5: Lower & Upper & Energy diff. & gf         & A  & dT"; 
  cout << "\n 6: Lower & Upper & Energy diff. & gf         & A"; 
  while ((table_format<1)||(table_format>6))
    {
      cout << "\n Specify table format (1, 2, 3, 4, 5, or 6)..........: ";
      cin >> table_format;
    }
  const unsigned int gfpos_in_table[6] = {6,5,5,5,4,4};
  unsigned int gfpos=gfpos_in_table[table_format-1];
  
  emin = 1E99; emax=-1E99;
  string line, helpstr;
  string *label = new string[npts];
  string *lower = new string[npts];
  double *gi = new double[npts];
  double *gilower = new double[npts];
  ifstream fin(filename); 
  unsigned int posL, posD, pos[8];
  double value, L;
  cout << endl;

  int n, nlower = 0;
  double gisum =0.0;
  for (int i=0; i<npts; i++) 
    {
      getline(fin,line);

      // search for '&' field separators
      pos[0]=0;
      for(int n=1; n<8; n++)
	{
	  pos[n] = line.find('&',pos[n-1]+1);
	  if (pos[n]==string::npos) break;
	}
      //label and weight of lower term
      label[i] = line.substr(0,pos[1]-1);
      L = Lsymb.find(label[i].at(label[i].length()-1));
      posL=label[i].find_last_of('_'); 
      helpstr=line.substr(posL+1,label[i].length()-posL-1);       
      sscanf(helpstr.data(),"%lf",&value); // value = 2*S+1      
      gi[i] = value*(2*L+1);

      // list of initial terms
      for (n=0; n < nlower; n++)
	{
	  if (label[i]==lower[n]) break;
	}
      if (n>=nlower) 
	{
	  lower[n]=label[i];
	  gilower[n] = gi[i];
	  gisum += gi[i];
	  nlower++;
	}

      //transition energy
      helpstr = line.substr(pos[2]+2,pos[3]-pos[2]-3);
      sscanf(helpstr.data(),"%lf",&value);
      energy[i] = value*invcm2eV;

      //weighted oscillator strength
      helpstr = line.substr(pos[gfpos-1]+2,pos[gfpos]-pos[gfpos-1]-3);
      posD=helpstr.find('D');
      helpstr.replace(posD,1,1,'E');
      sscanf(helpstr.data(),"%lf",&value);
      strength[i] = value*f2MbeV;

      //cout << energy[i] << " " << strength[i] <<endl;
      if (energy[i]>emax) emax=energy[i];
      if (energy[i]<emin) emin=energy[i];
    }
  fin.close();

  cout << "  " << nlower << " different inital terms found.";
  cout << " Which term is to be used? " << endl;
  cout << "    0: statistical average (Sum gi = " << gisum << ")" << endl;
  for (n=0; n<nlower; n++)
    {
      cout << "    " << n+1 << ": " << lower[n] << " (gi = " << gilower[n]  << ")" << endl;
    }
  cout << " Make a choice ........................................: ";
  int nini;
  cin >> nini;
   
  if ((nini==0)||(nini>nlower))
    { 
      for (int i=0; i<npts; i++) if (gisum>0.0) strength[i] /= gisum; 
    }
  else
    { 
      for (int i=0; i<npts; i++)
	{ 
	  if (label[i]==lower[nini])
	    {
	      if (gi[i]>0.0) strength[i] /= gi[i];
	    }
	  else
	    {
	      strength[i] = 0.0;
	    }
	}
    }
  delete[] gilower;
  delete[] gi;
  delete[] lower;
  delete[] label;
}

///////////////////////////////////////////////////////////////////////////////

int open_JAC(char *filename)
  {
  int npts;  
  ifstream fin(filename); 
  if (!fin) {return -1;}
  cout << endl << " Opening JAC sum file assuming 4 header lines and" << endl;
  cout << " alternating Coulomb and Babushkin data ..." << endl;
  npts=0;
  while (!fin.eof())
    {
      fin.ignore(1024,'\n');
    npts++;
    }
  npts--;
  fin.close();
													
  return (npts-5)/2; // 5 lines do not contain data, factor 2 because of Coulomb and Babushkin data
  } 

///////////////////////////////////////////////////////////////////////////////

void read_JAC(char *filename, double *energy, double *sigma, 
		    double *width, int npts, double &emin, double &emax)

  {
  int i, widthfactor=1;
  char answer;
  bool Babushkin_flag = false;
  ifstream fin(filename); 

  double *EinsteinA  = new double[npts];
  double *EinsteinB  = new double[npts];
  double *fosci      = new double[npts];

  cout << "\n Coulomb or Babuskin gauge? (C/B) : ";
  cin >> answer;
  if ((answer=='B') || (answer=='b'))
    {
      Babushkin_flag = true;
    }
  else
    {
      Babushkin_flag = false;
    }
  
  // skip file header
  for (i=0;i<4;i++) fin.ignore(1024,'\n');

  cout << endl;
  emin = 1E99; emax=-1E99;
  for (i=0; i<npts; i++) {
    if (Babushkin_flag) fin.ignore(1024,'\n'); // discard Coubomb data

    fin.ignore(40);
    fin >> energy[i];
    fin.ignore(27);
    fin >> EinsteinA[i] >> EinsteinB[i] >> fosci[i] >> width[i];
    fin.ignore(1024,'\n'); // read until end of line
    
    //cout << energy[i] << " " << EinsteinA[i] << " " << EinsteinB[i] << " " <<  fosci[i] << " " <<  width[i] << endl;
													       
    energy[i] *= efac;
    sigma[i] = 109.761*fosci[i]; // conversion to line strength in eV Mb
    cout << " " << energy[i] << " " << width[i] << " " << sigma[i] <<endl;
    if (energy[i]>emax) emax=energy[i];
    if (energy[i]<emin) emin=energy[i];

    if (!Babushkin_flag) fin.ignore(1024,'\n'); // discard Babushkin data
  }
  fin.close();

  delete[] fosci;
  delete[] EinsteinB;
  delete[] EinsteinA;
  }

///////////////////////////////////////////////////////////////////////////////

void info_DRtheory(void)
  {
  cout << "\n Predefined theory file formats are recognized by the extensions";
  cout << "\n\n *.lorentz, *.steih, *.bdr, *.brr, *.autos, *.rrautos, *.gautos, *.cowan,";
  cout << "\n *.lindroth, *.bnl, *.pindzola, *.gn, *.gnl, *.fritzsche, *.grasp2K\n";
  cout << "\n\n A theory data file with any other extension is expected to";
  cout << "\n contain lines with three numbers characterizing a single";
  cout << "\n DR or PI resonance each by its peak position, strength, and width.";
  cout << "\n A width of 0 eV is allowed. A negative width indicates that";
  cout << "\n a 1/E times Lorentzian line profile is to be used for the";
  cout << "\n specific DR resonance instead of a pure Lorentzian."; 
  cout << "\n After having read in the theory file the program will ask you";
  cout << "\n to supply the energy and cross section units you use.\n";
  cout.flush();
  }

///////////////////////////////////////////////////////////////////////////////

int open_DRtheory(char *filename, int &theo_mode, int &number_of_levels)
  {
    int npts = 0;
    theo_mode = 0;  
    number_of_levels = 1;

    efac = 1.0; wfac = 1.0; sfac = 1.0;

    if (strstr(filename,".lindroth")) 
      { theo_mode = 1; sfac = 1e-20; }
    else if (strstr(filename,".lorentz")) 
      { theo_mode = 1; efac=1; sfac = 1; }
    else if (strstr(filename,".bdr"))  
      { theo_mode = 2; efac = 13.606; sfac = 1e-18; }
    else if (strstr(filename,".pindzola")) 
      { theo_mode = 3; efac = 27.2116; }
    else if (strstr(filename,".bnl")) 
      { theo_mode = 4; efac = 1; sfac = 1e-18; }
    else if (strstr(filename,".cowan")) 
      { theo_mode = 5; efac = 1; sfac = 1; }
    else if (strstr(filename,".gautos")) 
      { theo_mode = 6; efac = 1; sfac = 1; }
    else if (strstr(filename,".gn")) 
      { theo_mode = 7; efac = 1; sfac = 1e-18; }
    else if (strstr(filename,".gnl")) 
      { theo_mode = 8; efac = 1; sfac = 1e-18; }
    else if (strstr(filename,".fritzsche")) 
      { theo_mode = 9; efac = 1; sfac = 1; }
    else if (strstr(filename,".brr"))  
      { theo_mode =10; efac = 13.606; sfac = 1e-18*efac; }
    else if (strstr(filename,".steih")) 
      { theo_mode = 11; efac=1; sfac = 1E-24; }
    else { theo_mode = 0; } 
    
    if (theo_mode==0)
      {
	cout << "\n The following theory data formats are defined:";
	cout << "\n      Lorentzian(e,S,wL)  or Voigt(Er,S,wl,wG) : 1";
	cout << "\n         wL > 0 : sigma(E) ~ Lorentzian";
	cout << "\n         wL < 0 : sigma(E) ~ Lorentzian/E";
	cout << "\n          S > 0 : norm. factor = Er";
	cout << "\n          S < 0 : norm. factor = Er+G^2/4Er";
	cout << '\n';
	cout << "\n            binned Xsec (autostructure) ......: 2";
	cout << '\n';
	cout << "\n            binned Xsec (Pindzola) ...........: 3";
	cout << '\n';
	cout << "\n                                  make a choice : ";
	cin >> theo_mode;
	switch (theo_mode)
	  {
	  case 1 : 
	    cout << "\n                         Give energy unit in eV : "; 
	    cin >> efac; wfac = efac;
	    cout << "\n                  Give strength unit in eV cm^2 : ";
	    cin >> sfac;
	    break;
	  case 2 : sfac = 1e-18; efac = 13.6056; break;
	  case 3 : efac = 27.2116; break;
	  default: break;
	  }
      }
    
    switch (theo_mode)
      {
      case 1 : npts = open_peak(filename); break;
      case 2 : npts = open_autostructure(filename, number_of_levels); break;
      case 3 : npts = open_Pindzola(filename); break;
      case 4 : npts = open_Badnell_nl(filename); break;
      case 5 : npts = open_Cowan(filename); break;
      case 6 : npts = open_Griffin_autostructure(filename); break;
      case 7 : npts = open_Griffin_n(filename); break;
      case 8 : npts = open_Griffin_nl(filename); break;
      case 9 : npts = open_Fritzsche(filename); break;
      case 10 : npts = open_autosrr(filename); break;
      case 11 : npts = open_peak(filename); break;
      default: theo_mode = 0;
      }
    return npts;
  }

///////////////////////////////////////////////////////////////////////////////

void read_DRtheory(char *filename, int theo_mode, int level_number, double *energy, double *strength,
                 double *wLorentz, int npts, double &emin, double &emax)
{
  double *qFano = new double [npts];
  double *wGauss = new double [npts];
  switch (theo_mode)
    {
    case 1 : read_peak(PEAK_Lorentz,filename,energy,strength,qFano,wLorentz,wGauss,npts,emin,emax); 
             break;
    case 2 : read_autostructure(filename,level_number,energy,strength,wLorentz,emin,emax); 
             break;
    case 3 : read_Pindzola(filename,energy,strength,wLorentz,emin,emax); 
             break;
    case 4 : read_Badnell_nl(filename,energy,strength,wLorentz,npts,emin,emax); 
             break;
    case 5 : read_Cowan(filename,energy,strength,wLorentz,emin,emax); 
             break;
    case 6 : read_Griffin_autostructure(filename,energy,strength,wLorentz,emin,emax);
             break;
    case 7 : read_Griffin_n(filename,energy,strength,wLorentz,emin,emax); 
             break; 
    case 8 : read_Griffin_nl(filename,energy,strength,wLorentz,emin,emax); 
             break; 
    case 9 : read_Fritzsche(filename,energy,strength,wLorentz,npts,emin,emax); 
             break; 
    case 10: read_autosrr(filename,energy,strength,wLorentz,emin,emax); 
             break;
    case 11: read_peak(PEAK_Lorentz_Steih,filename,energy,strength,qFano,wLorentz,wGauss,npts,emin,emax); 
             break;
    }
  delete[] wGauss;
  delete[] qFano;
}
                       

///////////////////////////////////////////////////////////////////////////////

void info_peak_files(void)
  {
  cout << "\n Predefined theory file formats are recognized by the extensions";
  cout << "\n *.grasp2K *.jac\n";
  cout << "\n A theory data file with any other extension is expected to have";
  cout << "\n at least two columns of resonance positions and strengths.";
  cout.flush();
  }

///////////////////////////////////////////////////////////////////////////////

int open_peak_file(char *filename, int &read_mode, int &number_of_levels)
{
  int npts = 0;
  read_mode = 0;  
  number_of_levels = 1;
  
  efac = 1.0; wfac = 1.0; sfac = 1.0;
  
  cout << "\n Cross sections will be calculated as a sum peaks.";
  cout << "\n The following peak functions are defined:";
  cout << "\n           Delta(Er,S) .......................: 1";
  cout << '\n';
  cout << "\n           Lorentz(Er,S,wL) ..................: 2";
  cout << "\n";
  cout << "\n           Fano(Er,S,q,wL) ...................: 3";
  cout << '\n';
  cout << "\n           Gauss(Er,S,wG) ....................: 4";
  cout << '\n';
  cout << "\n           Voigt(Er,S,wL,wG) .................: 5";
  cout << "\n";
  cout << "\n           FanoVoigt(Er,S,q,wL,wG) ...........: 6";
  cout << "\n";       
  cout << "\n The peak parameters are expected to be on one";
  cout << "\n line per peak and to be separeted by spaces.";
  cout << "\n Zero widths are allowed.";
  cout << "\n";       
  cout << "\n In addition, differently formatted output from";
  cout << "\n theoretical calculations is handled by the";
  cout << "\n subsequent menu items:";
  cout << '\n';
  cout << "\n   Grasp2K ascii output from rtransitiontable : 7";
  cout << '\n';
  cout << "\n   JAC sum output from excitation calculation : 8";
  cout << '\n';
  cout << "\n                                  make a choice : ";
  cin >> read_mode;
  cout << "\n                    Give name of peak data file : ";
  cin >> filename;
  
  if (read_mode==7)
    {
      npts = open_Grasp2K(filename); 
    }
  else if (read_mode==8)
    {
      npts = open_JAC(filename); 
    }
  else if (read_mode>0)
    {
      npts = open_peak(filename); 
    }
  return npts;
}


///////////////////////////////////////////////////////////////////////////////

int read_peak_file(char *filename, int read_mode, int level_number, double *energy, double *strength,
		    double *qFano, double *wLorentz, double *wGauss, int npts, double &emin, double &emax)
{
  int peak_shape=PEAK_undefined;  // for peak shape see enumeration in peakfunctions.h
  for (int i=0; i<npts; i++)
    {
      energy[i] = 0.0;
      strength[i] = 0.0;
      qFano[i] = 0.0;
      wLorentz[i] = 0.0;
      wGauss[i] = 0.0;
    }
  if  ((read_mode>=1) && (read_mode<=6))
    {
      peak_shape = read_mode;
      read_peak(read_mode,filename,energy,strength,qFano,wLorentz,wGauss,npts,emin,emax);
    }
  else if (read_mode==7)
    {
      peak_shape = PEAK_delta;
      read_Grasp2K(filename,energy,strength,npts,emin,emax); 
    }
  else if (read_mode==8)
    {
      peak_shape = PEAK_Lorentz;
      read_JAC(filename,energy,strength,wLorentz,npts,emin,emax); 
    }
  else
    {
      peak_shape = PEAK_undefined;
    }
  return peak_shape;
}
                       













