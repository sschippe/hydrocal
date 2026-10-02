/**
 * @file TOPbase.cxx
 *
 * @brief Handling of data from the TOPbase atomic data base
 *
 * $Id: TOPbase.cxx 2037 2026-07-17 15:21:56Z iamp $
// SPDX-License-Identifier: MIT
 */
#include "fele.h"
#include "hydroconst.h"
#include "hydromath.h"
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

double FWHM_ElHassan_Fe2(double eV);
double FWHM_ElHassan_Fe3(double eV);
double FWHM_ElHassan_Fe4(double eV);

void convolve_TOPbase_PI(void) {
  string line, filename;
  char answer;
  ifstream fin;
  ofstream fout;
  int Znucl, Nele;
  string asymb[] = {"H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",
                    "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
                    "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe"}; // 26
  string lsymb[] = {"S", "P", "D", "F", "G", "H", "I",
                    "K", "L", "M", "N", "O", "Q"}; // 13
  string psymb[] = {"e", "o"};

  cout << "\n Convolution of photoionization (PI) cross sections from the "
          "TOPbase.";
  cout << "\n Cross section files can be retrieved from the internet, see";
  cout << "\n       http://cdsweb.u-strasbg.fr/topbase/ftp.html.\n";
  cout << "\n The cross sections can be selected by term.";
  cout << "\n Cross sections can be written to a file before or after "
          "convolution.";
  cout << "\n A gaussian width can be specified for the convolution.";
  cout << "\n During the convolution the TOP cross sections are linearly "
          "interpolated.\n\n";
  cout << "\n Give name of TOBbase data file containing PI cross sections: ";
  cin >> filename;
  fin.open(filename);
  if (!fin.is_open()) {
    cout << "\n ERROR::File " << filename << " not found.\n\n";
    return;
  }

  if (!getline(fin, line)) {
    cout << "\n ERROR::File " << filename << " is empty.\n\n";
    return;
  }
  {
    istringstream ss(line);
    ss >> Znucl >> Nele;
  }
  int Qcharge = Znucl - Nele;

  cout << "\n nuclear charge .....: " << setw(2) << Znucl;
  cout << "\n number of electrons.: " << setw(2) << Nele;
  if (Znucl > 26) {
    cout << "\n ERROR: Nuclear charge > 26.\n\n";
    return;
  }
  cout << "\n ion ................: " << asymb[Znucl - 1] << Qcharge << "+\n";

  unsigned int Ssearch, Lsearch, Psearch, Nsearch;
  cout << "\n Give 2S+1, L, parity (0 or 1), and term number (1,2,...) ..: ";
  cin >> Ssearch >> Lsearch >> Psearch >> Nsearch;
  if (Psearch > 1) {
    cout << "\n ERROR: invalid parity value, please specifiy 0 or 1.";
    return;
  }
  if (Lsearch > 12) {
    cout << "\n ERROR: invalid L value, please specifiy one of 0,1,...,12.";
    return;
  }

  unsigned int Sread, Lread, Pread, Nread;
  int nlines, ndummy;
  std::vector<double> ecm, xsec;
  double elo = 1E99, ehi = -1E99;
  int search_flag = 1;
  while (search_flag) {
    if (!getline(fin, line))
      break;
    {
      istringstream ss(line);
      ss >> Sread >> Lread >> Pread >> Nread;
    }
    cout << "\n Reading cross section for term " << Sread << " " << Lread << " "
         << Pread << " " << Nread << " ...";
    if (!getline(fin, line))
      break;
    {
      istringstream ss(line);
      ss >> ndummy >> nlines;
    }
    if (!getline(fin, line))
      break;

    ecm.resize(nlines);
    xsec.resize(nlines);
    elo = 1E99, ehi = -1E99;
    for (int n = 0; n < nlines; n++) {
      getline(fin, line);
      {
        istringstream ss(line);
        ss >> ecm[n] >> xsec[n];
      }
      ecm[n] = ecm[n] * hydroconst::Ryd_eV;
      if (ecm[n] < elo)
        elo = ecm[n];
      if (ecm[n] > ehi)
        ehi = ecm[n];
    }
    if ((Sread == Ssearch) && (Lread == Lsearch) && (Pread == Psearch) &&
        (Nread == Nsearch)) {
      cout << "\n\n " << nlines << " cross section values found.\n";
      search_flag = 0;
    } else {
      xsec.clear();
      ecm.clear();
    }
  }
  fin.close();

  if (search_flag) {
    cout << "\n ERROR::Term not found.\n\n";
    return;
  }
  string term_string =
      to_string(Ssearch) + lsymb[Lsearch] + psymb[Psearch] + to_string(Nsearch);

  FELE fele = &felegauss;
  double fwhm = 0.0;
  int fwhm_mode = 0;
  cout << "\n Specify photon energy resolution mode";
  cout << "\n  0: unconvoluted cross section";
  cout << "\n  1: constant resolution";
  cout << "\n  2: relative resolution";
  cout << "\n  3: specific photon energy dependence";
  cout << "\n Make your choice ........................................: ";
  cin >> fwhm_mode;

  if (fwhm_mode == 0) { // output of unconvoluted cross section data
    cout << "\n Output of raw data or interpolated data (r/i) ...........: ";
    cin >> answer;
    int raw_flag = ((answer == 'r') || (answer == 'R'));

    cout << "\n Give name of output file ................................: ";
    cin >> filename;
    fout.open(filename);
    if (raw_flag) {
      fout << "# TOPbase raw PI cross sections\n";
    } else {
      fout << "# TOPbase interpolated PI cross sections\n";
    }
    fout << "# " << asymb[Znucl - 1] << Qcharge << "+ " << term_string
         << " (statistical weight: " << Ssearch * (2 * Lsearch + 1) << ")\n";
    fout << "      Erel           s" << term_string << "\n";
    fout << "      (eV)           (Mb)\n";
    if (raw_flag) {
      for (int n = 0; n < nlines; n++) {
        fout << setw(15) << uppercase << defaultfloat << setprecision(8)
             << ecm[n] << "   " << setw(15) << uppercase << defaultfloat
             << setprecision(8) << xsec[n] << "\n";
      }
      cout << "\n Raw cross sections written to " << filename;
    } else {
      double emin, emax, edelta;
      cout << "\n Input energy range in eV: [" << elo << ", " << ehi << "]";
      cout << "\n Give output energy range (emin,emax,delta) ..............: ";
      cin >> emin >> emax >> edelta;
      for (double eV = emin; eV <= emax; eV += edelta) {
        double xsec_interp = 0.0;
        if ((eV > elo) && (eV < ehi)) {
          int n = 0;
          for (int nn = 1; nn < nlines; nn++) {
            if (ecm[nn] > eV) {
              n = nn;
              break;
            }
          }
          if (n > 0) {
            xsec_interp = xsec[n - 1] + (xsec[n] - xsec[n - 1]) /
                                            (ecm[n] - ecm[n - 1]) *
                                            (eV - ecm[n - 1]);
          }
        }
        fout << setw(15) << uppercase << defaultfloat << setprecision(8) << eV
             << "   " << setw(15) << uppercase << defaultfloat
             << setprecision(8) << xsec_interp << "\n";
      }
      cout << "\n Linearly interpolated cross sections written to " << filename;
    }
    fout.close();
    return;
  } else if (fwhm_mode == 1) {
    cout << "\n Give constant fwhm in eV ................................: ";
    cin >> fwhm;
  } else if (fwhm_mode == 2) {
    cout << "\n Give relative fwhm (proportional to photon energy) ......: ";
    cin >> fwhm;
  } else {
    cout << "\n Which fwhm energy dependence to use?";
    cout << "\n 2: El Hassan et al., PRA 79 (2009) 033415: Fe2+";
    cout << "\n 3: El Hassan et al., PRA 79 (2009) 033415: Fe3+";
    cout << "\n 4: El Hassan et al., PRA 79 (2009) 033415: Fe4+";
    int choice;
    cout << "\n Make your choice ........................................: ";
    cin >> choice;
    if ((choice < 2) || (choice > 4)) {
      cout << "\n ERROR::Illeagal choice.\n\n";
      return;
    }
    fwhm_mode = 300 + choice;
  }

  cout << "\n Give name of output file ................................: ";
  cin >> filename;

  double emin, emax, edelta, eV;
  int steps;
  cout << "\n Input energy range in eV: [" << elo << ", " << ehi << "]";
  cout << "\n Give output energy range (emin,emax,delta) ..............: ";
  cin >> emin >> emax >> edelta;
  steps = int((emax - emin) / edelta) + 1;

  cout << "\n Convoluted cross sections will be calculated at " << steps
       << " energies.\n";
  cout << "\n Progress: ";
  fout.open(filename);
  fout << "# TOPbase convoluted PI cross sections\n";
  fout << "# " << asymb[Znucl - 1] << Qcharge << "+ " << term_string
       << " (statistical weight: " << Ssearch * (2 * Lsearch + 1) << ")\n";
  fout << "      Erel           fwhm        s" << term_string << "\n";
  fout << "      (eV)           (eV)        (Mb)\n";

  for (int i = 0; i < steps; i++) { // loop over energies
    if ((i % 100) == 0) {
      cout << "." << flush;
    }
    eV = emin + i * edelta;
    double conv_fwhm = fwhm;
    if (fwhm_mode == 2) {
      conv_fwhm = fwhm * eV;
    } else if (fwhm_mode == 302) {
      conv_fwhm = FWHM_ElHassan_Fe2(eV);
    } else if (fwhm_mode == 303) {
      conv_fwhm = FWHM_ElHassan_Fe3(eV);
    } else if (fwhm_mode == 304) {
      conv_fwhm = FWHM_ElHassan_Fe4(eV);
    }
    double conv_e;                        // integration variable
    double conv_width = 10.0 * conv_fwhm; // left and right from the gaussian
    double conv_delta =
        0.05 * conv_fwhm; // step width for stepping through the gaussian
    double conv_emin = eV - conv_width;
    double conv_emax = eV + conv_width;

    if ((conv_emax < elo) || (conv_emin > ehi)) {
      fout << setw(15) << uppercase << defaultfloat << setprecision(8) << eV
           << "   " << setw(15) << uppercase << defaultfloat << setprecision(8)
           << conv_fwhm << "   " << setw(15) << uppercase << defaultfloat
           << setprecision(8) << 0.0 << "\n";
      continue;
    }

    // make conv_delta smaller than the minimum raw energy difference
    for (int n = 0; n < (nlines - 1); n++) {
      if ((ecm[n] > conv_emin) && (ecm[n] < conv_emax)) {
        if (0.1 * (ecm[n + 1] - ecm[n]) < conv_delta)
          conv_delta = 0.1 * (ecm[n + 1] - ecm[n]);
      }
    }

    double conv_sum = 0.0;
    for (conv_e = conv_emin; conv_e <= conv_emax; conv_e += conv_delta) {
      double xsec_interp = 0.0;
      if ((conv_e > elo) && (conv_e < ehi)) {
        int n = 0;
        for (int nn = 1; nn < nlines; nn++) {
          if (ecm[nn] > conv_e) {
            n = nn;
            break;
          }
        }
        if (n > 0) {
          xsec_interp = xsec[n - 1] + (xsec[n] - xsec[n - 1]) /
                                          (ecm[n] - ecm[n - 1]) *
                                          (conv_e - ecm[n - 1]);
        }
      }
      conv_sum +=
          conv_delta * fele(eV, conv_e, conv_fwhm, conv_fwhm) * xsec_interp;
    }

    if (fabs(conv_sum) < 1E-99)
      conv_sum = 0.0;
    fout << setw(15) << uppercase << defaultfloat << setprecision(8) << eV
         << "   " << setw(15) << uppercase << defaultfloat << setprecision(8)
         << conv_fwhm << "   " << setw(15) << uppercase << defaultfloat
         << setprecision(8) << conv_sum << "\n";
  }
  fout.close();
  cout << "\n Convoluted cross sections written to " << filename;
  cout << "\n\n";
}

double FWHM_ElHassan_Fe2(double eV) {
  if (eV < 30.0) {
    return 0.1;
  } else if ((eV >= 30.0) && (eV <= 60.0)) {
    return 0.1 + 0.77 / 30.0 * (eV - 30.0);
  } else if ((eV > 60.0) && (eV <= 100)) {
    return 0.24 + 0.63 / 40.0 * (eV - 60.0);
  } else if ((eV > 100.0) && (eV <= 160.0)) {
    return 1.67 + 3.301 / 60.0 * (eV - 100.0);
  } else {
    return 4.98;
  }
}

double FWHM_ElHassan_Fe3(double eV) {
  if (eV < 30.0) {
    return 0.06;
  } else if ((eV >= 30.0) && (eV <= 45.0)) {
    return 0.06 + 0.19 / 15.0 * (eV - 30.0);
  } else if ((eV > 45.0) && (eV <= 65.0)) {
    return 0.13 + 0.30 / 20.0 * (eV - 45.0);
  } else if ((eV > 65.0) && (eV <= 100.0)) {
    return 0.30 + 0.60 / 35.0 * (eV - 65.0);
  } else if ((eV > 100.0) && (eV <= 160.0)) {
    return 1.67 + 3.301 / 60.0 * (eV - 100.0);
  } else {
    return 4.98;
  }
}

double FWHM_ElHassan_Fe4(double eV) {
  if (eV < 60.0) {
    return 0.13;
  } else if ((eV >= 60.0) && (eV <= 80.0)) {
    return 0.13 + 0.14 / 20.0 * (eV - 60.0);
  } else if ((eV > 80.0) && (eV <= 140.0)) {
    return 0.49 + 1.38 / 60.0 * (eV - 80.0);
  } else {
    return 1.87;
  }
}
