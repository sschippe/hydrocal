/**
 * @file readxsec.cxx
 *
 * @brief Reading theoretical cross sections from files
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: readxsec.cxx 2038 2026-07-17 16:41:52Z iamp $
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */
#include "readxsec.h"
#include "CoolerFractions.h"
#include "hydroconst.h"
#include "peakfunctions.h"
#include "readadas.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

using namespace std;

double efac = 1.0, wfac = 1.0, sfac = 1.0, dummy;
int nmax = 0;

ADASadf09 *adf09;

///////////////////////////////////////////////////////////////////////////////

void read_vector(string line, int, vector<double> &v, int pos0, int no) {
  string zahl;
  for (int k = 0; k < no; k++) {
    zahl = line.substr(k * 12, 12);
    istringstream buf(zahl, istringstream::in);
    buf >> v[pos0 + k];
  }
}

///////////////////////////////////////////////////////////////////////////////

static void skip_comment_lines(ifstream &fin, string &line) {
  while (!fin.eof()) {
    streampos pos = fin.tellg();
    getline(fin, line);
    if (line.empty() || line[0] == '#')
      continue;
    fin.seekg(pos);
    break;
  }
}

///////////////////////////////////////////////////////////////////////////////

int open_peak(string filename) {
  int npts = 0;
  string line;
  ifstream fin(filename);
  if (!fin) {
    return -1;
  }

  while (!fin.eof()) {
    streampos pos = fin.tellg();
    getline(fin, line);
    if (line.empty() || line[0] == '#')
      continue;
    npts++;
  }
  fin.close();
  return npts;
}

///////////////////////////////////////////////////////////////////////////////

void read_peak(int mode, string filename, vector<double> &energy,
               vector<double> &strength, vector<double> &qFano,
               vector<double> &wLorentz, vector<double> &wGauss, int npts,
               double &emin, double &emax) {
  ifstream fin(filename);
  string line;

  emin = 1e99;
  emax = -1e99;
  for (int i = 0; i < npts; i++) {
    skip_comment_lines(fin, line);
    qFano[i] = 0.0;
    wGauss[i] = 0.0;
    if (mode == PEAK_Lorentz) {
      fin >> energy[i] >> strength[i] >> wLorentz[i];
      cout << energy[i] << "  " << strength[i] << "  " << wLorentz[i] << endl;
    } else if (mode == PEAK_Lorentz_Steih) {
      fin >> energy[i] >> wLorentz[i] >> strength[i];
      cout << energy[i] << "  " << strength[i] << "  " << wLorentz[i] << endl;
    } else if (mode == PEAK_Fano) { // Fano peaks
      fin >> energy[i] >> strength[i] >> qFano[i] >> wLorentz[i];
      cout << energy[i] << "  " << strength[i] << "  " << qFano[i] << "  "
           << wLorentz[i] << endl;
    } else if (mode == PEAK_Voigt) {
      fin >> energy[i] >> strength[i] >> wLorentz[i] >> wGauss[i];
      cout << energy[i] << "  " << strength[i] << "  " << wLorentz[i] << "  "
           << wGauss[i] << endl;
    } else if (mode == PEAK_FanoVoigt) { // FanoVoigt peaks
      fin >> energy[i] >> strength[i] >> qFano[i] >> wLorentz[i] >> wGauss[i];
      cout << energy[i] << "  " << strength[i] << "  " << qFano[i] << "  "
           << wLorentz[i] << "  " << wGauss[i] << endl;
    } else if (mode == PEAK_Lorentz_Amaro) {
      // fin >> energy[i] >> dummy >> dummy >> dummy >> dummy >> dummy >>
      // wLorentz[i] >> strength[i];
      fin >> energy[i] >> dummy >> dummy >> strength[i] >> dummy >> dummy >>
          wLorentz[i];
      fin.ignore(numeric_limits<streamsize>::max(),
                 '\n'); // skip the reaminder of the line
      cout << energy[i] << "  " << strength[i] << "  " << wLorentz[i] << endl;
    } else { // PEAK_delta
      fin >> energy[i] >> strength[i];
      cout << energy[i] << "  " << strength[i] << endl;
    }
    // energy[i] *= efac;
    if (energy[i] < (emin))
      emin = energy[i];
    if (energy[i] > (emax))
      emax = energy[i];
    // wLorentz[i] *= wfac;
    // strength[i] *= sfac;
  }
  fin.close();
  cout.flush();
}

///////////////////////////////////////////////////////////////////////////////

int open_Griffin_autostructure(string filename) {
  int npts = 0, incr = 0, n, l;
  string line, teststr1, teststr2;
  istringstream buf;

  ifstream fin(filename);
  if (!fin) {
    return -1;
  }

  nmax = 0; // nmax is a global variable in this module
  while (!fin.eof()) {
    getline(fin, line);
    teststr1 = line.substr(1, 5);
    teststr2 = line.substr(5, 3);

    if (teststr1 == "total") {
      getline(fin, line);
      getline(fin, line);
      incr = 0;
    } else if (teststr2 == "n=:") {
      buf.str(line);
      buf.clear();
      buf.ignore(8);
      buf >> n;
      buf.ignore(6);
      buf >> l;
      if (n > nmax)
        nmax = n; // nmax is a global variable in this module
      getline(fin, line);
      incr = 1;
    } else if (incr) {
      npts++;
    }
  }
  fin.close();

  return npts;
}

///////////////////////////////////////////////////////////////////////////////

void read_Griffin_autostructure(string filename, vector<double> &energy,
                                vector<double> &strength, vector<double> &width,
                                double &emin, double &emax) {
  int n, l, fdim = nmax * (nmax + 1) / 2;
  vector<double> fraction(fdim, 1.0);
  string header;
  char answer;
  cout << " Read surviving fractions from file? (y/n) : ";
  cin >> answer;
  if ((answer == 'y') || (answer == 'Y')) {
    nmax = readfraction(fraction, header, nmax);
    cout << " surviving fractions " << header << '\n';
    cout.flush();
    cout << "\n Perform l-averaging of fractions? (y/n) ..: ";
    cin >> answer;
    if ((answer == 'y') || (answer == 'Y')) {
      for (n = 1; n <= nmax; n++) {
        double sum = 0.0, n2 = n * n;
        for (l = 0; l < n; l++) {
          sum += (2.0 * l + 1.0) * fraction[n * (n - 1) / 2 + l];
        }
        double laverage = sum / n2;
        for (l = 0; l < n; l++) {
          fraction[n * (n - 1) / 2 + l] = laverage;
        }
        cout << " n = " << n << " : " << laverage << "\n";
        cout.flush();
      }
    }
  }

  int npts = -1, incr = 0, badcount = 0;
  string line, teststr1, teststr2;
  istringstream buf;
  emin = 1E99;
  emax = -1E99;
  ifstream fin(filename);

  while (!fin.eof()) {
    getline(fin, line);
    teststr1 = line.substr(1, 5);
    teststr2 = line.substr(5, 3);

    if (teststr1 == "total") {
      getline(fin, line);
      getline(fin, line);
      incr = 0;
    } else if (teststr2 == "n=:") {
      buf.str(line);
      buf.clear();
      buf.ignore(8);
      buf >> n;
      buf.ignore(6);
      buf >> l;
      getline(fin, line);
      incr = 1;
    } else if (incr) {
      double e, s;
      buf.str(line);
      buf.clear();
      if (buf >> e >> s) {
        npts++;
        e *= efac;
        s *= sfac;
        if (e > (emax))
          emax = e;
        if (e < (emin))
          emin = e;
        energy[npts] = e;
        strength[npts] = s * fraction[n * (n - 1) / 2 + l];
        width[npts] = 0.0;
      } else {
        badcount++;
      }
    }
  }
  if (badcount) {
    cout << "\n " << badcount << " reading errors\n";
    cout.flush();
  }
  fin.close();
}

///////////////////////////////////////////////////////////////////////////////

int open_Pindzola(string filename) {
  int npts;
  double dmy;
  ifstream fin(filename);
  if (!fin) {
    return -1;
  }
  fin >> dmy >> dmy >> dmy >> dmy;
  fin >> dmy >> dmy;
  fin >> dmy >> dmy;
  fin >> npts >> dmy >> dmy >> dmy;
  fin.close();

  return npts;
}

///////////////////////////////////////////////////////////////////////////////

void read_Pindzola(string filename, vector<double> &energy,
                   vector<double> &sigma, vector<double> &width, double &emin,
                   double &emax) {
  const double Rydberg = hydroconst::Ryd_eV;
  int npts;
  double dmy, field;
  ifstream fin(filename);
  fin >> dmy >> dmy >> dmy >> dmy;
  fin >> dmy >> dmy;
  fin >> dmy >> dmy;
  fin >> npts >> field >> dummy >> dummy;

  emin = 1e99;
  emax = -1e99;
  int i;
  for (i = 0; i < npts; i++) {
    fin >> energy[i];
    energy[i] *= efac;
    if (energy[i] < emin)
      emin = energy[i];
    if (energy[i] > emax)
      emax = energy[i];
  }
  for (i = 0; i < npts; i++) {
    fin >> sigma[i];
    sigma[i] *= sfac * (energy[1] - energy[0]);
    width[i] = 0.0;
  }
  fin.close();

  string header;
  char answer{};
  cout << " Read surviving fractions from file? (y/n) : ";
  cin >> answer;
  if ((answer == 'y') || (answer == 'Y')) {
    int n, l;
    double slim, z, elo, ehi, fn;

    cout << "\n Give effective ionic charge ..............: ";
    cin >> z;
    cout << "\n Give series limit in eV ..................: ";
    cin >> slim;
    cout << "\n Give maximum n ...........................: ";
    cin >> n;
    n = n + 1;
    int fdim = n * (n + 1) / 2;
    vector<double> fraction(fdim, 1.0);
    int nmax_loc = readfraction(fraction, header, n - 1);
    cout << " surviving fractions " << header << '\n';
    cout.flush();

    elo = slim - Rydberg * z * z / (n - 0.5) / (n - 0.5);
    ehi = slim;
    fn = 0.0;
    for (i = npts - 1; i >= 0; i--) {
      if (energy[i] < elo) {
        n--;
        elo = slim - Rydberg * z * z / (n - 0.5) / (n - 0.5);
        ehi = slim - Rydberg * z * z / (n + 0.5) / (n + 0.5);
        fn = 0.0;
        for (l = 0; l < n; l++) {
          fn += (2 * l + 1) * fraction[n * (n - 1) / 2 + l];
        }
        fn /= (n * n);
        cout << "n=" << n << ": fraction = " << fn << "\n";
        cout.flush();
      }
      if (n > nmax_loc) {
        fn = 0.0;
      } // else {fn = 1.0;}
      if (fn == 1.0)
        break;
      if (sigma[i] > 0) {
        sigma[i] *= fn;
        cout << n << " : " << elo << " " << energy[i] << " " << ehi << "\n";
      }
    } // end for(i...)
  } // end if

  cout << "\n Calculate integrated rate coeff? (y/n) ...: ";
  cin >> answer;
  if ((answer == 'y') || (answer == 'Y')) {
    double elo, ehi, strength = 0.0;
    cout << "\n Give boundaries of integration range in eV: ";
    cin >> elo >> ehi;
    for (i = 0; i < npts; i++) {
      if ((energy[i] >= elo) && (energy[i] <= ehi))
        strength += sigma[i] * hydroconst::clight_cm_s *
                    sqrt(energy[i] / hydroconst::mec2_eV);
    }
    cout << "\n" << field << " V/cm : " << strength << " eVcm^2\n";
    ofstream fout("integrals.dat", ios::app);
    fout << elo << " " << ehi << " " << field << " " << strength << "\n";
    fout.close();
  }
  cout.flush();
}

///////////////////////////////////////////////////////////////////////////////

int open_Cowan(string filename) {
  int npts;
  ifstream fin(filename);
  if (!fin)
    return -1;
  fin >> npts;
  fin.close();
  return npts;
}

///////////////////////////////////////////////////////////////////////////////

void read_Cowan(string filename, vector<double> &energy, vector<double> &sigma,
                vector<double> &width, double &emin, double &emax) {
  int i, n, l, npts;
  string line;
  istringstream buf;

  ifstream fin(filename);
  getline(fin, line);
  buf.str(line);
  buf >> npts;

  vector<int> n_store(npts, 0);
  vector<int> l_store(npts, 0);

  nmax = 0; // nmax is a global variable in this module
  for (i = 0; i < npts; i++) {
    getline(fin, line);
    buf.str(line);
    buf >> n_store[i] >> l_store[i] >> energy[i] >> sigma[i] >> width[i];
    //    cout <<  n_store[i] << ' ' << l_store[i] << ' ' << energy[i] << ' ';
    //    cout <<  sigma[i] << ' ' << width[i] << '\n';
    if (n_store[i] > nmax) {
      nmax = n_store[i];
    }
  }
  fin.close();

  int fdim = nmax * (nmax + 1) / 2;
  vector<double> fraction(fdim, 1.0);
  string header;
  char answer;
  cout << " Read surviving fractions from file? (y/n) : ";
  cin >> answer;
  if ((answer == 'y') || (answer == 'Y')) {
    nmax = readfraction(fraction, header, nmax);
    cout << " surviving fractions " << header << '\n';
    cout.flush();
  }

  vector<double> psum(nmax + 1, 0.0);
  vector<double> wsum(nmax + 1, 0.0);
  emin = 1e99;
  emax = -1e99;
  for (i = 0; i < npts; i++) {
    n = n_store[i];
    l = l_store[i];
    energy[i] *= efac;
    if (energy[i] < emin)
      emin = energy[i];
    if (energy[i] > emax)
      emax = energy[i];
    sigma[i] *= sfac;
    wsum[n] += sigma[i];
    sigma[i] *= fraction[n * (n - 1) / 2 + l];
    psum[n] += sigma[i];
  }

  string fn = filename + "_n";
  ofstream fout(fn);
  for (n = 1; n <= nmax; n++) {
    double pmean = wsum[n] > 0 ? psum[n] / wsum[n] : 0.0;
    fout << n << ' ' << pmean << '\n';
  }
  fout.close();
  cout << " n-specific detection probabilities written to file " << fn << '\n';
  cout.flush();
}

///////////////////////////////////////////////////////////////////////////////

int open_autosdr(string filename, int &number_of_levels) {
  int npts;
  ifstream fin(filename);
  if (!fin) {
    return -1;
  }
  fin >> npts;

  int line_count = 0;
  string line;
  while (!fin.eof()) {
    fin >> line;
    line_count++;
  }
  fin.close();
  number_of_levels = 6 * line_count / npts - 1;
  return npts - 1;
}

///////////////////////////////////////////////////////////////////////////////

void read_autosdr(string filename, int level_number, vector<double> &energy,
                  vector<double> &strength, vector<double> &width, double &emin,
                  double &emax) {
  int npts;
  ifstream fin(filename);
  fin >> npts;
  string line;
  int i, j, k, pos0 = 0;
  int line_count = 0;

  for (i = 0; i < npts / 6; i++) {
    fin >> line;
    line_count++;
    read_vector(line, 200, energy, pos0, 6);
    pos0 += 6;
  }
  if (npts % 6 > 0) {
    fin >> line;
    read_vector(line, 200, energy, pos0, (npts % 6));
    if (npts % 6 > 1)
      line_count++;
  }
  npts = npts - 1; // energy bins are given by min and max value,
  // from here on we use the center of the bin as energy value
  // therefore the number of bins in decreased by one

  // overread information on lower levels that is not requested
  // 'level number' is a global variable set in open_autostructure()
  for (j = 1; j < level_number; j++) {
    for (k = 0; k < line_count; k++)
      fin >> line;
  }

  pos0 = 0;
  for (i = 0; i < npts / 6; i++) {
    fin >> line;
    read_vector(line, 200, strength, pos0, 6);
    pos0 += 6;
  }
  if (npts % 6 > 0) {
    fin >> line;
    read_vector(line, 200, strength, pos0, npts % 6);
  }
  fin.close();

  // strcat(filename,".col");
  // ofstream fout(filename);
  for (i = 0; i < npts; i++) {
    width[i] = 0.0;
    double binwidth =
        (energy[i + 1] - energy[i]) *
        efac; // not that we have one more energy point than cross section point
    energy[i] =
        energy[i] * efac + 0.5 * binwidth; // use the center of a bin as energy
    strength[i] *= sfac * binwidth;
    // fout << energy[i] << " " << strength[i] << endl;
  }
  // fout.close();
  emin = energy[0];
  emax = energy[npts - 1];
}

///////////////////////////////////////////////////////////////////////////////

int open_autosrr(string filename) {
  int npts;
  ifstream fin(filename);
  if (!fin)
    return -1;
  fin >> npts;
  fin.close();
  return npts;
}

///////////////////////////////////////////////////////////////////////////////

void read_autosrr(string filename, vector<double> &energy,
                  vector<double> &strength, vector<double> &width, double &emin,
                  double &emax) {
  int npts;
  ifstream fin(filename);
  fin >> npts;
  for (int j = 0; j < npts; j++) {
    fin >> energy[j] >> strength[j];
    energy[j] *= efac;
    strength[j] *= sfac;
    width[j] = 0.0;
    cout << energy[j] << " " << strength[j] << endl;
  }
  fin.close();

  emin = energy[0];
  emax = energy[npts - 1];
}

///////////////////////////////////////////////////////////////////////////////

int open_Badnell_nl(string filename) {
  int npts = 0, n, l;
  double e, s;
  string line;
  istringstream buf;

  ifstream fin(filename);
  if (!fin)
    return -1;
  for (int i = 1; i <= 3; i++) {
    getline(fin, line);
  }
  while (!fin.eof()) {
    getline(fin, line);
    buf.str(line);
    if (buf >> n >> l >> e >> s) {
      npts++;
    }
  }
  fin.close();
  return npts;
}

///////////////////////////////////////////////////////////////////////////////

void read_Badnell_nl(string filename, vector<double> &energy,
                     vector<double> &sigma, vector<double> &width, int npts,
                     double emin, double emax) {
  int i = 0, n = 0, l = 0, nn = 0, ll = 0;
  double e, s, binwidth;
  vector<int> n_store(npts, 0); // storage for n of energy bins
  vector<int> l_store(npts, 0); // storage for l of energy bins
  string line;
  istringstream buf;

  ifstream fin(filename);

  getline(fin, line);
  buf.str(line);
  buf.seekg(11);
  buf >> binwidth;
  getline(fin, line);
  getline(fin, line);
  i = 0;
  nmax = 0; // namx is a global variable in this module
  while (!fin.eof()) {
    getline(fin, line);
    buf.str(line);
    buf >> nn >> ll >> e >> s; // read 4 numbers
    if (buf.fail()) { // in case of error only 2 numbers are there, which are to
                      // be interpreted as n and l
      buf.clear();    // clear error condition
      n = nn;
      l = ll;
      if (n > nmax)
        nmax = n;
    } else {
      energy[i] = e;
      sigma[i] = s;
      n_store[i] = n;
      l_store[i] = l;
      i++;
      if (i == npts)
        break;
    }
  }
  fin.close();

  int fdim = nmax * (nmax + 1) / 2;
  vector<double> fraction(fdim, 1.0);
  string header;
  char answer;
  cout << " Read surviving fractions from file? (y/n) : ";
  cin >> answer;
  if ((answer == 'y') || (answer == 'Y')) {
    nmax = readfraction(fraction, header, nmax);
    cout << " surviving fractions " << header << '\n';
    cout.flush();
  }

  vector<double> psum(nmax + 1, 0.0);
  vector<double> wsum(nmax + 1, 0.0);
  emin = 1e99;
  emax = -1e99;
  for (i = 0; i < npts; i++) {
    n = n_store[i];
    l = l_store[i];
    //    if (n<10)
    //      {cout << n << ' ' << l << ' ' << fraction[(n-1)*(n-1)+l] << '\n';}
    energy[i] *= efac;
    if (energy[i] < emin)
      emin = energy[i];
    if (energy[i] > emax)
      emax = energy[i];
    sigma[i] *= sfac * binwidth;
    wsum[n] += sigma[i];
    sigma[i] *= fraction[(n - 1) * n / 2 + l];
    psum[n] += sigma[i];
    width[i] = 0.0;
  }

  string fn = filename + "_n";
  ofstream fout(fn);
  for (n = 1; n <= nmax; n++) {
    double pmean = wsum[n] > 0 ? psum[n] / wsum[n] : 0.0;
    fout << n << ' ' << pmean << '\n';
  }
  fout.close();
  cout << " n-specific detection probabilities written to file " << fn << '\n';
  cout.flush();
}

///////////////////////////////////////////////////////////////////////////////

int open_Griffin_n(string filename) {
  int n, npts = 0;
  double d1, d2, d3, d4, d5, d6;

  string line;
  istringstream buf;
  ifstream fin(filename);
  if (!fin)
    return -1;

  for (int i = 0; i < 5; i++)
    getline(fin, line);
  while (!fin.eof()) {
    getline(fin, line);
    buf.str(line);
    if (buf >> n >> d1 >> d2 >> d3 >> d4 >> d5 >> d6) {
      npts++;
      //    cout << n << "\n";cout.flush();
    }
  }
  fin.close();

  nmax = n; // nmax is a global variable within this module
  return 2 * npts;
}

///////////////////////////////////////////////////////////////////////////////

void read_Griffin_n(string filename, vector<double> &energy,
                    vector<double> &strength, vector<double> &width,
                    double emin, double emax) {
  int i, n, l, fdim = nmax * (nmax + 1) / 2;
  vector<double> fraction(fdim, 1.0);
  string header;
  char answer;
  cout << " Read surviving fractions from file? (y/n) : ";
  cin >> answer;
  if ((answer == 'y') || (answer == 'Y')) {
    nmax = readfraction(fraction, header, nmax);
    cout << " surviving fractions " << header << '\n';
    cout.flush();

    for (n = 1; n <= nmax; n++) {
      double sum = 0.0, n2 = n * n;
      for (l = 0; l < n; l++) {
        sum += (2.0 * l + 1.0) * fraction[n * (n - 1) / 2 + l];
      }
      double laverage = sum / n2;
      for (l = 0; l < n; l++) {
        fraction[n * (n - 1) / 2 + l] = laverage;
      }
      cout << " n = " << n << " : " << laverage << '\n';
      cout.flush();
    }
  }

  string line;
  istringstream buf;

  emin = 1E99;
  emax = -1E99;

  double deltae, e1, s1, e2, s2, et, st;
  ifstream fin(filename);

  getline(fin, line);
  getline(fin, line);
  getline(fin, line);
  buf.str(line);
  buf.seekg(47);
  buf >> deltae;
  getline(fin, line);
  getline(fin, line);

  i = 0;
  while ((!fin.eof()) && (2 * i < nmax)) {
    getline(fin, line);
    buf.str(line);
    if (buf >> n >> e1 >> s1 >> e2 >> s2 >> et >> st) {
      e1 *= efac;
      s1 *= sfac * deltae;
      e2 *= efac;
      s2 *= sfac * deltae;
      if (e1 > (emax))
        emax = e1;
      if (e1 < (emin))
        emin = e1;
      if (e2 > (emax))
        emax = e2;
      if (e2 < (emin))
        emin = e2;
      energy[2 * i] = e1;
      strength[2 * i] = s1 * fraction[n * (n - 1) / 2];
      width[2 * i] = 0.0;
      energy[2 * i + 1] = e2;
      strength[2 * i + 1] = s2 * fraction[n * (n - 1) / 2];
      width[2 * i + 1] = 0.0;
      // cout <<energy[2*i]<<" "<<strength[2*i]<<" "<< width[2*i]<<"\n";
      // cout <<energy[2*i+1]<<" "<<strength[2*i+1]<<" "<< width[2*i+1]<<"\n";
    }
    i++;
  }
  cout.flush();
  fin.close();
}

///////////////////////////////////////////////////////////////////////////////

int open_Griffin_nl(string filename) {
  int n, l, npts = 0;
  double d1, d2;

  string line;
  istringstream buf;
  ifstream fin(filename);
  if (!fin) {
    return -1;
  }

  getline(fin, line);
  getline(fin, line);
  while (!fin.eof()) {
    getline(fin, line);
    buf.str(line);
    if (buf >> n >> l >> d1 >> d2) {
      npts++;
      //    cout << n << "\n";cout.flush();
    }
  }
  fin.close();

  nmax = n; // nmax is a global variable within this module
  return npts;
}

///////////////////////////////////////////////////////////////////////////////

void read_Griffin_nl(string filename, vector<double> &energy,
                     vector<double> &strength, vector<double> &width,
                     double emin, double emax) {
  int i, n, l, fdim = nmax * (nmax + 1) / 2;
  vector<double> fraction(fdim, 1.0);
  string header;
  char answer;
  cout << "\n Read surviving fractions from file? (y/n) : ";
  cin >> answer;
  if ((answer == 'y') || (answer == 'Y')) {
    readfraction(fraction, header, nmax);
    cout << " surviving fractions " << header << '\n';
    cout.flush();
    cout << "\n Perform l-averaging of fractions? (y/n) ..: ";
    cin >> answer;
    if ((answer == 'y') || (answer == 'Y')) {
      for (n = 1; n <= nmax; n++) {
        double sum = 0.0, n2 = n * n;
        for (l = 0; l < n; l++) {
          sum += (2.0 * l + 1.0) * fraction[n * (n - 1) / 2 + l];
        }
        double laverage = sum / n2;
        for (l = 0; l < n; l++) {
          fraction[n * (n - 1) / 2 + l] = laverage;
        }
        cout << " n = " << n << " : " << laverage << '\n';
        cout.flush();
      }
    }
  }

  string line;
  istringstream buf;

  emin = 1E99;
  emax = -1E99;

  double e, s;
  ifstream fin(filename);

  getline(fin, line);
  getline(fin, line);

  i = 0;
  n = 0;
  l = 0;
  while ((!fin.eof()) && (n <= nmax) && (l < nmax - 1)) {
    getline(fin, line);
    buf.str(line);
    if (buf >> n >> l >> e >> s) {
      e *= efac;
      s *= sfac;
      if (e > (emax))
        emax = e;
      if (e < (emin))
        emin = e;
      energy[i] = e;
      strength[i] = s * fraction[n * (n - 1) / 2 + l];
      width[i] = 0.0;
      // cout << n << " " << l << " " << energy[i] << " " << strength[i] <<
      // "\n";
    }
    i++;
  }
  cout.flush();
  fin.close();
}

///////////////////////////////////////////////////////////////////////////////

int open_Fritzsche(string filename) {
  int npts;
  string line;
  ifstream fin(filename);
  if (!fin)
    return -1;
  npts = 0;
  while (!fin.eof()) {
    getline(fin, line);
    npts++;
  }
  npts--;
  fin.close();
  return npts;
}

///////////////////////////////////////////////////////////////////////////////

void read_Fritzsche(string filename, vector<double> &energy,
                    vector<double> &sigma, vector<double> &width, int npts,
                    double &emin, double &emax) {
  int i, Babuskin = 0, widthfactor = 1;
  ifstream fin(filename);

  vector<double> Aami(npts, 0.0);
  vector<double> Aa(npts, 0.0);
  vector<double> ArBab(npts, 0.0);
  vector<double> ArCoul(npts, 0.0);
  vector<double> SBab(npts, 0.0);
  vector<double> SCoul(npts, 0.0);

  cout << "\n Coulomb or Babuskin gauge? (0=Coul/1=Bab/2=negCoul/3=negBab) : ";
  cin >> Babuskin;
  if (Babuskin == 2) {
    widthfactor = -1;
    Babuskin = 0;
    // cout<<endl<<"Widht are multiplied by -1"<<endl;
  } else if (Babuskin == 3) {
    widthfactor = -1;
    Babuskin = 1;
    // endl cout<<"Widht are multiplied by -1"<<endl;
  }

  emin = 1E99;
  emax = -1E99;
  for (i = 0; i < npts; i++) {
    fin.ignore(40);
    fin >> energy[i] >> Aami[i] >> Aa[i] >> ArBab[i] >> ArCoul[i] >> SBab[i] >>
        SCoul[i];
    // cout <<  energy[i] << ' ' << Aami[i] << ' ' << Aa[i] << ' ';
    // cout <<  ArBab[i] << ' ' << ArCoul[i] << ' ' <<  SBab[i] << ' ' <<
    // SCoul[i] << '\n';
    if (Babuskin) {
      // cout<<"Babuskin gauge used"<<endl;
      width[i] = widthfactor * 6.58212e-16 * (Aa[i] + ArBab[i]);
      sigma[i] = SBab[i] * sfac;
    } else {
      // cout<<"Coulomb gauge used"<<endl;
      width[i] = widthfactor * 6.58212e-16 * (Aa[i] + ArCoul[i]);
      sigma[i] = SCoul[i] * sfac;
    }
    energy[i] *= efac;
    // cout<<energy[i]<<" "<<width[i]<<" "<<sigma[i]<<endl;
    if (energy[i] > emax)
      emax = energy[i];
    if (energy[i] < emin)
      emin = energy[i];
  }
  fin.close();
}

///////////////////////////////////////////////////////////////////////////////

int open_Grasp2K(string filename) {
  int npts;
  string line;
  ifstream fin(filename);
  if (!fin)
    return -1;
  npts = 0;
  while (!fin.eof()) {
    getline(fin, line);
    npts++;
  }
  npts--;
  fin.close();
  return npts;
}

///////////////////////////////////////////////////////////////////////////////
void read_Grasp2K(string filename, vector<double> &energy,
                  vector<double> &strength, int npts, double &emin,
                  double &emax) {
  const double invcm2eV = 1.239842e-4; // conversion from 1/cm to eV
  const double f2MbeV = 109.761;       // conversion from oscillator strength to
                                       // integrated cross section
  const string Lsymb = "SPDFGHIKLMNOPQRSTUVWXYZ";

  int table_format = 0;
  cout << "\n Which table format has been used in Grasp2K rtransitiontable?";
  cout << "\n 1: Lower & Upper & Energy diff. & wavelength & S  & gf & A  & dT";
  cout << "\n 2: Lower & Upper & Energy diff. & wavelength & gf & A  & dT";
  cout << "\n 3: Lower & Upper & Energy diff. & wavelength & gf & A";
  cout << "\n 4: Lower & Upper & Energy diff. & S          & gf & A  & dT";
  cout << "\n 5: Lower & Upper & Energy diff. & gf         & A  & dT";
  cout << "\n 6: Lower & Upper & Energy diff. & gf         & A";
  while ((table_format < 1) || (table_format > 6)) {
    cout << "\n Specify table format (1, 2, 3, 4, 5, or 6)..........: ";
    cin >> table_format;
  }
  const unsigned int gfpos_in_table[6] = {6, 5, 5, 5, 4, 4};
  unsigned int gfpos = gfpos_in_table[table_format - 1];

  emin = 1E99;
  emax = -1E99;
  string line, helpstr;
  vector<string> label(npts, "");
  vector<string> lower(npts, "");
  vector<double> gi(npts, 0.0);
  vector<double> gilower(npts, 0.0);
  ifstream fin(filename);
  unsigned int posL, posD, pos[8];
  double value, L;
  cout << endl;

  int n, nlower = 0;
  double gisum = 0.0;
  for (int i = 0; i < npts; i++) {
    getline(fin, line);

    // search for '&' field separators
    pos[0] = 0;
    for (int m = 1; m < 8; m++) {
      pos[m] = static_cast<unsigned int>(line.find('&', pos[m - 1] + 1));
      if (pos[m] == string::npos)
        break;
    }
    // label and weight of lower term
    label[i] = line.substr(0, pos[1] - 1);
    L = static_cast<double>(Lsymb.find(label[i].at(label[i].length() - 1)));
    posL = static_cast<unsigned int>(label[i].find_last_of('_'));
    helpstr = line.substr(posL + 1, label[i].length() - posL - 1);
    value = stod(helpstr);
    gi[i] = value * (2 * L + 1);

    // list of initial terms
    for (n = 0; n < nlower; n++) {
      if (label[i] == lower[n])
        break;
    }
    if (n >= nlower) {
      lower[n] = label[i];
      gilower[n] = gi[i];
      gisum += gi[i];
      nlower++;
    }

    // transition energy
    helpstr = line.substr(pos[2] + 2, pos[3] - pos[2] - 3);
    value = stod(helpstr);
    energy[i] = value * invcm2eV;

    // weighted oscillator strength
    helpstr = line.substr(pos[gfpos - 1] + 2, pos[gfpos] - pos[gfpos - 1] - 3);
    posD = static_cast<unsigned int>(helpstr.find('D'));
    helpstr.replace(posD, 1, 1, 'E');
    value = stod(helpstr);
    strength[i] = value * f2MbeV;

    // cout << energy[i] << " " << strength[i] <<endl;
    if (energy[i] > emax)
      emax = energy[i];
    if (energy[i] < emin)
      emin = energy[i];
  }
  fin.close();

  cout << "  " << nlower << " different inital terms found.";
  cout << " Which term is to be used? " << endl;
  cout << "    0: statistical average (Sum gi = " << gisum << ")" << endl;
  for (n = 0; n < nlower; n++) {
    cout << "    " << n + 1 << ": " << lower[n] << " (gi = " << gilower[n]
         << ")" << endl;
  }
  cout << " Make a choice ........................................: ";
  int nini;
  cin >> nini;

  if ((nini == 0) || (nini > nlower)) {
    for (int i = 0; i < npts; i++)
      if (gisum > 0.0)
        strength[i] /= gisum;
  } else {
    for (int i = 0; i < npts; i++) {
      if (label[i] == lower[nini - 1]) {
        if (gi[i] > 0.0)
          strength[i] /= gi[i];
      } else {
        strength[i] = 0.0;
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
bool JAC_photon_flag = false;
const double JAC_cutoff =
    1E-40; // partial line strength below the cutoff will be ignored

int open_JAC(string filename) {
  ifstream fin(filename);
  if (!fin) {
    return -1;
  }

  string line;
  string teststringr1 = "Partial resonance strength";
  string teststringr2 =
      "Total Auger rates, radiative rates and resonance strengths";
  int npts = 0, counter1 = 0, counter2 = 0;
  double dvalue;
  vector<double> dvector1, dvector2;
  char answer;

  cout << endl << " Opening JAC summary file ... " << endl;

  while (!fin.eof()) {
    getline(fin, line);

    // count partial data
    if (line.find(teststringr1) != string::npos) {
      cout << "    counting partial resonance strengths ..." << endl;
      counter1 = 1;
    }
    if ((counter1 > 0) && (counter1 < 7)) {
      counter1++; // skip header
    } else if (line.find("----------") != string::npos) {
      counter1 = 0; // end of partial data reached
    }
    if (counter1 >= 7) {
      {
        istringstream iss(line);
        iss.ignore(127);
        iss >> dvalue;
      }
      // cout << line << ":::" << dvalue << endl;
      if (dvalue > JAC_cutoff) {
        dvector1.push_back(dvalue);
      }
    }

    // count total data
    if (line.find(teststringr2) != string::npos) {
      cout << "    counting total rates ..." << endl;
      counter2 = 1;
    }
    if ((counter2 > 0) && (counter2 < 26)) { // used to be 8
      counter2++;                            // skip header
    } else if (line.find("----------") != string::npos) {
      counter2 = 0; // end of total data
    }
    if (counter2 >= 26) { // used to be 8
      {
        istringstream iss(line);
        iss.ignore(132);
        iss >> dvalue;
      }
      // cout << line << ":::" << dvalue << endl;
      dvector2.push_back(dvalue);
    }
  } // end while

  cout << " Photon spectra or DR cross sections ? (p/D) : ";
  cin >> answer;
  if ((answer == 'p') || (answer == 'P')) {
    JAC_photon_flag = true;
    npts = static_cast<int>(dvector1.size());
  } else {
    JAC_photon_flag = false;
    npts = static_cast<int>(dvector2.size());
  }

  fin.close();
  return npts;
}

void read_JAC(string filename, vector<double> &energy, vector<double> &sigma,
              vector<double> &width, int npts, double &emin, double &emax) {
  ifstream fin(filename);
  vector<int> pindex_m, tindex_m;
  vector<double> Ephoton, SphotonC, SphotonB;
  vector<double> Eres, SresC, SresB, WresC, WresB;
  int counter1 = 0, counter2 = 0;
  string line;
  string teststringr1 = "Partial resonance strength";
  string teststringr2 =
      "Total Auger rates, radiative rates and resonance strengths";
  int ivalue, nptot = 0, nttot = 0;
  double dvalue1, dvalue2, dvalue3, dvalue4, dvalue5;
  char answer;

  cout << endl << " Opening JAC summary file ... " << endl;

  while (!fin.eof()) {
    getline(fin, line);

    // read partial data
    if (line.find(teststringr1) != string::npos) {
      cout << "    reading partial resonance strengths ..." << endl;
      counter1 = 1;
    }
    if ((counter1 > 0) && (counter1 < 7)) {
      counter1++; // skip header
    } else if (line.find("----------") != string::npos) {
      counter1 = 0; // end of partial data reached
    }
    if (counter1 >= 7) {
      {
        istringstream iss(line);
        iss.ignore(11);
        iss >> ivalue;
        iss.ignore(65);
        iss >> dvalue1;
        iss.ignore(8);
        iss >> dvalue2 >> dvalue3;
      }
      // cout << line << endl;
      // cout << ivalue << " " << dvalue1  << " " << dvalue2  << " " << dvalue3
      // << endl;
      if (dvalue3 > JAC_cutoff) {
        pindex_m.push_back(ivalue);
        Ephoton.push_back(dvalue1);
        SphotonC.push_back(dvalue2);
        SphotonB.push_back(dvalue3);
      }
    }
    nptot = static_cast<int>(pindex_m.size());

    // read total data
    if (line.find(teststringr2) != string::npos) {
      cout << "    reading total rates ..." << endl;
      counter2 = 1;
    }
    if ((counter2 > 0) && (counter2 < 26)) { // used to be 8
      counter2++;                            // skip header
    } else if (line.find("----------") != string::npos) {
      counter2 = 0; // end of total data
    }
    if (counter2 >= 26) { // used to be 8
      {
        istringstream iss(line);
        iss.ignore(11);
        iss >> ivalue;
        iss.ignore(28);
        iss >> dvalue1;
        iss.ignore(50);
        iss >> dvalue2 >> dvalue3 >> dvalue4;
        iss.ignore(4);
        iss >> dvalue5;
      }
      // cout << line << endl;
      // cout << ivalue << " " << dvalue1  << " " << dvalue2  << " " << dvalue3
      // << "  "  << dvalue4 << " " << dvalue5 << endl;
      tindex_m.push_back(ivalue);
      Eres.push_back(dvalue1);
      SresC.push_back(dvalue2);
      SresB.push_back(dvalue3);
      WresC.push_back(dvalue4);
      WresB.push_back(dvalue5);
    }
  } // end while
  nttot = static_cast<int>(tindex_m.size());

  cout << "\n Coulomb or Babuskin gauge? (C/B) : ";
  cin >> answer;
  cout << endl;
  bool Babushkin_flag = ((answer == 'B') || (answer == 'b'));
  emin = 1E99;
  if (JAC_photon_flag) {
    if (nptot != npts) {
      cout << "ERROR:: inconsistent number of partial data: " << pindex_m.size()
           << " vs. " << npts << endl;
      return;
    }
    for (int i = 0; i < npts; i++) {
      energy[i] = Ephoton[i];
      sigma[i] =
          Babushkin_flag
              ? SphotonB[i]
              : SphotonC[i]; // line strength times  partial width in eV2 cm2
      int n = 0;             // next we look for the appropriate total width
      for (n = 0; n < nttot; n++) {
        if (tindex_m[n] == pindex_m[i])
          break;
      }
      if (n < nttot)
        width[i] = Babushkin_flag ? WresB[n] : WresC[n];
      else
        width[i] = 0.0; // no matching total width found
      sigma[i] =
          (width[i] > 0.0)
              ? sigma[i] / width[i]
              : 0.0; // divide by total width to obtain line stength in eV cm2
      // cout  << "  " << energy[i] << "   " << width[i] << "   " << sigma[i] <<
      // "  " << pindex_m[i] << endl;
    }
  } else {
    if (nttot != npts) {
      cout << "ERROR:: inconsistent number of total data: " << tindex_m.size()
           << " vs. " << npts << endl;
      return;
    }
    for (int i = 0; i < npts; i++) {
      energy[i] = Eres[i];
      sigma[i] =
          Babushkin_flag ? SresB[i] : SresC[i]; // line strength in eV cm2
      width[i] = Babushkin_flag ? WresB[i] : WresC[i];
      // cout  << "  " << energy[i] << "   " << width[i] << "   " << sigma[i] <<
      // endl;
    }
  }

  emin = 1e99;
  emax = -1e99;
  efac = 1.0;
  sfac = 1.0;

  cout << endl << " First 10 and last 10 entries in list of ";
  if (JAC_photon_flag) {
    cout << " photon emission lines:" << endl;
  } else {
    cout << " dielectronic-recombination resonances:" << endl;
  }
  cout << "   energy (eV)   width (eV)   strength (cm2 eV) " << endl;

  for (int i = 0; i < npts; i++) {
    energy[i] *= efac;
    sigma[i] *= sfac;
    width[i] *= efac;
    if (energy[i] > emax)
      emax = energy[i];
    if (energy[i] < emin)
      emin = energy[i];
    if ((i < 10) || (i >= (npts - 10))) {
      cout << "    " << energy[i] << "         " << width[i] << "        "
           << sigma[i] << endl;
    }
    if (i == 10)
      cout << endl;
  }
}
///////////////////////////////////////////////////////////////////////////////

void info_DRtheory(void) {
  cout << "\n Predefined theory file formats are recognized by the extensions";
  cout << "\n *.lorentz (columns with Lorentzian peak parameters eres, "
          "strength, width),";
  cout
      << "\n *.bdr, *.brr (binned DR or RR cross sections from AUTOSTRUCTURE),";
  cout << "\n *.bnl, *cowan, *.fritzsche, *.gautos, *.gn, *gnl, *.jac, "
          "*.lindroth, *.pindzola,";
  cout << "\n *.steih *.amaro (specific file formats from various authors, see "
          "readxsec.cxx).";
  cout << "\n\n A theory data file with any other extension is expected to "
          "contain lines with three or four numbers "
          "characterizing";
  cout << "\n a single DR or PI resonance each by its peak position, strength, "
          "Lorentzian width (and Gaussian width).";
  cout << "\n Widths of 0 eV are allowed. A negative Lorentzian width "
          "indicates that a 1/E times Lorentzian line profile";
  cout << "\n is to be used for this specific DR resonance instead of a pure "
          "Lorentzian."
       << endl;
  cout.flush();
}

///////////////////////////////////////////////////////////////////////////////

int open_DRtheory(string filename, int &theo_format, int &number_of_levels) {
  int npts = 0;
  theo_format = 0;
  number_of_levels = 1;

  efac = 1.0;
  wfac = 1.0;
  sfac = 1.0;

  if (strstr(filename.data(), ".lindroth")) {
    theo_format = THEO_FORMAT_LORENTZ;
    sfac = 1e-20;
  } else if (strstr(filename.data(), ".lorentz")) {
    theo_format = THEO_FORMAT_LORENTZ;
    efac = 1;
    sfac = 1;
  } else if (strstr(filename.data(), ".bdr")) {
    theo_format = THEO_FORMAT_AUTOSDR;
    efac = hydroconst::Ryd_eV;
    sfac = 1e-18;
  } else if (strstr(filename.data(), ".pindzola")) {
    theo_format = THEO_FORMAT_PINDZOLA;
    efac = hydroconst::Ryd_eV * 2.0;
  } else if (strstr(filename.data(), ".bnl")) {
    theo_format = THEO_FORMAT_BADNELL;
    efac = 1;
    sfac = 1e-18;
  } else if (strstr(filename.data(), ".cowan")) {
    theo_format = THEO_FORMAT_COWAN;
    efac = 1;
    sfac = 1;
  } else if (strstr(filename.data(), ".gautos")) {
    theo_format = THEO_FORMAT_GRIFFIN_AS;
    efac = 1;
    sfac = 1;
  } else if (strstr(filename.data(), ".gn")) {
    theo_format = THEO_FORMAT_GRIFFIN_N;
    efac = 1;
    sfac = 1e-18;
  } else if (strstr(filename.data(), ".gnl")) {
    theo_format = THEO_FORMAT_GRIFFIN_NL;
    efac = 1;
    sfac = 1e-18;
  } else if (strstr(filename.data(), ".fritzsche")) {
    theo_format = THEO_FORMAT_FRITZSCHE;
    efac = 1;
    sfac = 1;
  } else if (strstr(filename.data(), ".jac")) {
    theo_format = THEO_FORMAT_JAC;
    efac = 1;
    sfac = 1;
  } else if (strstr(filename.data(), ".brr")) {
    theo_format = THEO_FORMAT_AUTOSRR;
    efac = hydroconst::Ryd_eV;
    sfac = 1e-18 * efac;
  } else if (strstr(filename.data(), ".steih")) {
    theo_format = THEO_FORMAT_STEIH;
    efac = 1;
    sfac = 1E-24;
  } else if (strstr(filename.data(), ".amaro")) {
    theo_format = THEO_FORMAT_AMARO;
    efac = 1;
    sfac = 1E-20;
    wfac = hydroconst::hbar_eV_s;
  } else if (strstr(filename.data(), ".adf09")) {
    theo_format = THEO_FORMAT_ADASADF09;
    efac = 1;
    sfac = 1E-20;
    wfac = hydroconst::hbar_eV_s;
  } else {
    cout << "\n Assuming  3 or 4 columns of peak parameters e, S, wL, (wG)";
    cout << "\n  for Lorentzian(e,S,wL) or Voigt(Er,S,wl,wG) peaks";
    cout << "\n         wL > 0 : sigma(E) ~ Lorentzian";
    cout << "\n         wL < 0 : sigma(E) ~ Lorentzian/E";
    cout << "\n          S > 0 : norm. factor = Er";
    cout << "\n          S < 0 : norm. factor = Er+G^2/4Er";
    cout << endl;
    cout << "\n                         Give energy unit in eV : ";
    cin >> efac;
    wfac = efac;
    cout << "\n                  Give strength unit in eV cm^2 : ";
    cin >> sfac;
    theo_format = THEO_FORMAT_LORENTZ;
  }
  cout << " theo_format: " << theo_format << endl;
  switch (theo_format) {
  case THEO_FORMAT_LORENTZ:
    npts = open_peak(filename);
    break;
  case THEO_FORMAT_AUTOSDR:
    npts = open_autosdr(filename, number_of_levels);
    break;
  case THEO_FORMAT_PINDZOLA:
    npts = open_Pindzola(filename);
    break;
  case THEO_FORMAT_BADNELL:
    npts = open_Badnell_nl(filename);
    break;
  case THEO_FORMAT_COWAN:
    npts = open_Cowan(filename);
    break;
  case THEO_FORMAT_GRIFFIN_AS:
    npts = open_Griffin_autostructure(filename);
    break;
  case THEO_FORMAT_GRIFFIN_N:
    npts = open_Griffin_n(filename);
    break;
  case THEO_FORMAT_GRIFFIN_NL:
    npts = open_Griffin_nl(filename);
    break;
  case THEO_FORMAT_FRITZSCHE:
    npts = open_Fritzsche(filename);
    break;
  case THEO_FORMAT_JAC:
    npts = open_JAC(filename);
    break;
  case THEO_FORMAT_AUTOSRR:
    npts = open_autosrr(filename);
    break;
  case THEO_FORMAT_STEIH:
    npts = open_peak(filename);
    break;
  case THEO_FORMAT_AMARO:
    npts = open_peak(filename);
    break;
  case THEO_FORMAT_ADASADF09:
    adf09 = new ADASadf09(filename);
    npts = adf09->get_temperatures().size();
    number_of_levels = adf09->get_nlvl();
    break;
  default:
    theo_format = 0;
  }
  return npts;
}

///////////////////////////////////////////////////////////////////////////////

void read_DRtheory(string filename, int theo_format, int level_number,
                   vector<double> &energy, vector<double> &strength,
                   vector<double> &wLorentz, int &npts, double &emin,
                   double &emax, double eshift) {
  vector<double> qFano(npts, 0.0);
  vector<double> wGauss(npts, 0.0);
  switch (theo_format) {
  case THEO_FORMAT_LORENTZ:
    read_peak(PEAK_Lorentz, filename, energy, strength, qFano, wLorentz, wGauss,
              npts, emin, emax);
    break;
  case THEO_FORMAT_AUTOSDR:
    read_autosdr(filename, level_number, energy, strength, wLorentz, emin,
                 emax);
    break;
  case THEO_FORMAT_PINDZOLA:
    read_Pindzola(filename, energy, strength, wLorentz, emin, emax);
    break;
  case THEO_FORMAT_BADNELL:
    read_Badnell_nl(filename, energy, strength, wLorentz, npts, emin, emax);
    break;
  case THEO_FORMAT_COWAN:
    read_Cowan(filename, energy, strength, wLorentz, emin, emax);
    break;
  case THEO_FORMAT_GRIFFIN_AS:
    read_Griffin_autostructure(filename, energy, strength, wLorentz, emin,
                               emax);
    break;
  case THEO_FORMAT_GRIFFIN_N:
    read_Griffin_n(filename, energy, strength, wLorentz, emin, emax);
    break;
  case THEO_FORMAT_GRIFFIN_NL:
    read_Griffin_nl(filename, energy, strength, wLorentz, emin, emax);
    break;
  case THEO_FORMAT_FRITZSCHE:
    read_Fritzsche(filename, energy, strength, wLorentz, npts, emin, emax);
    break;
  case THEO_FORMAT_JAC:
    read_JAC(filename, energy, strength, wLorentz, npts, emin, emax);
    break;
  case THEO_FORMAT_AUTOSRR:
    read_autosrr(filename, energy, strength, wLorentz, emin, emax);
    break;
  case THEO_FORMAT_STEIH:
    read_peak(PEAK_Lorentz_Steih, filename, energy, strength, qFano, wLorentz,
              wGauss, npts, emin, emax);
    break;
  case THEO_FORMAT_AMARO:
    read_peak(PEAK_Lorentz_Amaro, filename, energy, strength, qFano, wLorentz,
              wGauss, npts, emin, emax);
    break;
  case THEO_FORMAT_ADASADF09: {
    vector<double> const &te = adf09->get_temperatures();
    unsigned int nte = te.size();
    // Convert temperature to energy (kT in eV)
    for (unsigned int i = 0; i < nte && i < (unsigned int)npts; i++) {
      energy[i] = te[i] * hydroconst::kB_eV_K;
    }
    if (level_number >= 0 && level_number < (int)adf09->get_nlvl()) {
      vector<double> const &alf = adf09->get_alfi(level_number);
      for (unsigned int i = 0;
           i < nte && i < (unsigned int)npts && i < alf.size(); i++) {
        strength[i] = alf[i] * sfac;
        wLorentz[i] = 0.0;
      }
    }
    npts = nte;
    emin = energy[0];
    emax = energy[npts - 1];
    break;
  }
  }

  if (fabs(eshift) < 1.0E-9)
    return;

  // apply overall energy shift
  std::vector<int> mark_negative(npts);
  for (int n = 0; n < npts; n++) {
    energy[n] += eshift;
    mark_negative[n] = (energy[n] < 0.0);
  }

  int nn = 0;
  for (int n = 0; n < npts; n++) {
    if (mark_negative[n]) {
      nn++;
      continue;
    }
    energy[n - nn] = energy[n];
    strength[n - nn] = strength[n];
    wLorentz[n - nn] = wLorentz[n];
  }
  npts -= nn;
  emax = 0.0;
  emin = 999999.9;
  for (int n = 0; n < npts; n++) {
    if (emax < energy[n])
      emax = energy[n];
    if (emin > energy[n])
      emin = energy[n];
  }
}

///////////////////////////////////////////////////////////////////////////////

void info_peak_files(void) {
  cout << "\n Predefined theory file formats are recognized by the extensions";
  cout << "\n *.grasp2K *.jac\n";
  cout << "\n A theory data file with any other extension is expected to have";
  cout << "\n at least two columns of resonance positions and strengths.";
  cout.flush();
}

///////////////////////////////////////////////////////////////////////////////

int open_peak_file(string &filename, int &read_mode, int &number_of_levels) {
  int npts = 0;
  read_mode = 0;
  number_of_levels = 1;

  efac = 1.0;
  wfac = 1.0;
  sfac = 1.0;

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

  if (read_mode == 7) {
    npts = open_Grasp2K(filename);
  } else if (read_mode == 8) {
    npts = open_JAC(filename);
  } else if (read_mode > 0) {
    npts = open_peak(filename);
  }
  return npts;
}

///////////////////////////////////////////////////////////////////////////////

int read_peak_file(string filename, int read_mode, int, vector<double> &energy,
                   vector<double> &strength, vector<double> &qFano,
                   vector<double> &wLorentz, vector<double> &wGauss, int npts,
                   double &emin, double &emax) {
  int peak_shape =
      PEAK_undefined; // for peak shape see enumeration in peakfunctions.h

  energy.assign(npts, 0.0);
  strength.assign(npts, 0.0);
  qFano.assign(npts, 0.0);
  wLorentz.assign(npts, 0.0);
  wGauss.assign(npts, 0.0);

  if ((read_mode >= 1) && (read_mode <= 6)) {
    peak_shape = read_mode;
    read_peak(read_mode, filename, energy, strength, qFano, wLorentz, wGauss,
              npts, emin, emax);
  } else if (read_mode == 7) {
    peak_shape = PEAK_delta;
    read_Grasp2K(filename, energy, strength, npts, emin, emax);
  } else if (read_mode == 8) {
    peak_shape = PEAK_Lorentz;
    read_JAC(filename, energy, strength, wLorentz, npts, emin, emax);
  } else {
    peak_shape = PEAK_undefined;
  }
  return peak_shape;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief averages Lorentzian cross section over histogram bins
 *
 * @param nL number of Lorentzian peaks
 * @param eL array of Lorentzian resonance energies in eV
 * @param sL array of Lorentzian resonance strengths in cm2 eV
 * @param wL array of Lorentzian resonance widths in eV
 * @param nbin number of bins
 * @param emin energy of first bin
 * @param bin_width uniform width of all bins in eV
 * @param binned_energy on exit array of bin energies (center of bin) in eV
 * @param binned_xesc on exit array of binned cross section (averaged per bin)
 * in cm2
 *
 */
void bin_lorentzian_peaks(int nL, vector<double> eL, vector<double> sL,
                          vector<double> wL, int nbin, double emin,
                          double bin_width, vector<double> &binned_energy,
                          vector<double> &binned_xsec) {
  const double eps = 1E-9;
  double Pi = 4.0 * atan(1.0);

  for (int i = 0; i < nbin; i++) {
    double ebmin = emin + i * bin_width; // min energy of bin
    if (ebmin < 0)
      continue;
    double ebmax = ebmin + bin_width; // max energy of bin
    binned_energy[i] = ebmin + 0.5 * bin_width;
    binned_xsec[i] = 0.0;
    for (int n = 0; n < nL; n++) {
      double eres = eL[n];
      double wres = fabs(wL[n]);
      if (wres < eps) { // treat narrow resonances as delta-like resonances
        if ((eres >= ebmin) && (eres < ebmax))
          binned_xsec[i] += sL[n];
      } else { // integrate the Lorentzian resonance from lower to upper bin
               // boundary
        double xmin = 2.0 * (ebmin - eres) / wres;
        double xmax = 2.0 * (ebmax - eres) / wres;
        binned_xsec[i] += sL[n] * (atan(xmax) - atan(xmin)) / Pi;
      }
    } // end for(n...
    binned_xsec[i] /= bin_width;
  } // end for(i...
  return;
}
