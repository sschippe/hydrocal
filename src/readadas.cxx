/**
 * @file readadas.cxx
 *
 * @brief Reads and processes ADAS ADF09 files
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */
#include "readadas.h"
#include "hydroconst.h"
#include "hydromath.h"
#include "buildinfo.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include "stdin_guard.h"

using namespace std;

static int read_e10(string const &line, int pos, int nmax,
                    vector<double> &out) {
  int n = 0;
  for (int i = 0; i < nmax; i++) {
    int start = pos + i * 10;
    if (start + 10 > (int)line.size())
      break;
    string field = line.substr(start, 10);
    for (auto &c : field)
      if (c == 'D' || c == 'd')
        c = 'E';
    char *end;
    double v = strtod(field.c_str(), &end);
    if (end == field.c_str())
      break;
    out.push_back(v);
    n++;
  }
  return n;
}

static vector<double> interpolate(vector<double> const &tx,
                                  vector<double> const &ty,
                                  vector<double> const &new_t) {
  size_t n = tx.size();
  vector<double> out(new_t.size(), 0.0);
  if (n < 2) {
    out.assign(new_t.size(), n == 1 ? ty[0] : 0.0);
    return out;
  }

  // Transform to log-log space; clamp small/negative y to a floor
  vector<double> lx(n), ly(n);
  double ymin = 1e-40;
  for (size_t i = 0; i < n; i++) {
    lx[i] = log10(tx[i]);
    ly[i] = (ty[i] > ymin) ? log10(ty[i]) : log10(ymin);
  }

  // Natural cubic spline in log-log space
  vector<double> y2(n);
  spline((int)n, lx, ly, y2);

  // Evaluate spline at requested temperatures
  for (size_t i = 0; i < new_t.size(); i++) {
    double t = new_t[i];
    if (t <= tx[0]) {
      out[i] = ty[0];
    } else if (t >= tx[n - 1]) {
      out[i] = ty[n - 1];
    } else {
      double lt = log10(t);
      double ly_out = splint(lt, (int)n, lx, ly, y2);
      out[i] = pow(10.0, ly_out);
    }
  }
  return out;
}

ADASadf09::ADASadf09(string filename) {
  ifstream fin;

  s_filename = filename;

  ui_npts = 0;
  fin.open(filename);
  if (!fin.is_open()) {
    cout << "ERROR:: file " << filename << " does not exist!" << endl;
    return;
  }
  string line;
  stringstream ss;
  getline(fin, line);
  ss.str("");
  ss << line;
  ss.ignore(line.size(), '\'');
  ss.get(c3_isoseq, 3, '\'');
  ss.ignore(line.size(), '=');
  ss >> ui_nuccharge;
  ss.ignore(line.size(), '/');
  ss.get(c3_coupling, 3, '/');
  cout << " processing ADAS ADF09 file " << filename << ": " << c3_isoseq
       << " isoelectronic sequence, Z=" << ui_nuccharge << "  " << c3_coupling
       << endl;

  getline(fin, line);
  getline(fin, line);
  ss.clear();
  ss.str("");
  ss << line;
  ss.ignore(line.size(), '=');
  ss >> d_bwnp;
  ss.ignore(line.size(), '=');
  ss >> ui_nprnti;
  ss.ignore(line.size(), '=');
  ss >> ui_nprntf;

  getline(fin, line);
  getline(fin, line);
  getline(fin, line);

  vui_indpp.resize(ui_nprntf);
  vui_isp.resize(ui_nprntf);
  vui_ilp.resize(ui_nprntf);
  vd_xjp.resize(ui_nprntf);
  vd_wnpi.resize(ui_nprntf);
  vs_cfgp.clear();
  for (unsigned int i = 0; i < ui_nprntf; i++) {
    string helpstr(20, '\0');
    getline(fin, line);
    sscanf(line.data(), "%6u%*10c%20c%*1c%1u%*1c%1u%*1c%4lf%*1c%11lf",
           &vui_indpp[i], helpstr.data(), &vui_isp[i], &vui_ilp[i], &vd_xjp[i],
           &vd_wnpi[i]);
    vs_cfgp.push_back(helpstr);
  }

  getline(fin, line);
  getline(fin, line);
  ss.clear();
  ss.str("");
  ss << line;
  ss.ignore(line.size(), '=');
  ss >> d_bwnr;
  ss.ignore(line.size(), '=');
  ss >> ui_nlvl;
  cout << "   reading data for " << ui_nlvl << " levels ..." << endl;

  getline(fin, line);
  getline(fin, line);
  getline(fin, line);

  vui_indx.resize(ui_nlvl);
  vui_indp.resize(ui_nlvl);
  vui_is.resize(ui_nlvl);
  vui_il.resize(ui_nlvl);
  vd_xj.resize(ui_nlvl);
  vd_wnrl.resize(ui_nlvl);
  for (int n = 0; n < 5; n++) {
    vd_aalp[n].resize(ui_nlvl);
  }
  vs_cfgl.clear();

  for (unsigned int i = 0; i < ui_nlvl; i++) {
    string helpstr(20, '\0');
    char test;
    getline(fin, line);
    sscanf(line.data(), "%6u%6u%*4c%20c%*1c%1u%*1c%1u%*1c%4lf%*1c%11lf%1c",
           &vui_indx[i], &vui_indp[i], helpstr.data(), &vui_is[i], &vui_il[i],
           &vd_xj[i], &vd_wnrl[i], &test);

    vs_cfgl.push_back(helpstr);

    int nhelp = vui_indp[i] - 1;
    if ((nhelp > 0) && (test == '*')) {
      ui_npts++;
      ss.clear();
      ss.str("");
      ss << line;
      ss.seekg(60);
      // vd_aalp holds at most 5 parent levels; clamp to avoid overrunning it
      for (int n = 0; (n < nhelp) && (n < 5); n++) {
        ss >> vd_aalp[n][i];
      }
    }
  }

  // --- skip NL-shell indexing ---
  while (getline(fin, line)) {
    if (line.find("NLREP") != string::npos)
      break;
  }
  while (getline(fin, line)) {
    if (line.empty() || line.find_first_not_of(" \t\r") == string::npos)
      break;
  }

  // --- skip N-shell indexing ---
  while (getline(fin, line)) {
    if (line.find("NREP") != string::npos)
      break;
  }
  while (getline(fin, line)) {
    if (line.empty() || line.find_first_not_of(" \t\r") == string::npos)
      break;
  }

  // --- PRTI sections ---
  vd_te.clear();
  vd_alfi.clear();
  vd_alff.clear();
  vd_alft.clear();

  for (unsigned int iprt = 0; iprt < ui_nprnti; iprt++) {
    while (getline(fin, line)) {
      if (line.find("PRTI=") != string::npos)
        break;
    }

    getline(fin, line);

    getline(fin, line);
    if (fin.eof())
      break;
    read_e10(line, 11, 10, vd_te);

    streampos saved = fin.tellg();
    getline(fin, line);
    int extra = read_e10(line, 11, 10, vd_te);
    if (extra == 0) {
      fin.clear();
      fin.seekg(saved);
    }

    unsigned int nte = vd_te.size();
    if (nte == 0) {
      cerr << "Warning: no temperature values found in PRTI block" << endl;
      break;
    }

    for (unsigned int i = 0; i < ui_nlvl; i++) {
      getline(fin, line);
      if (!fin)
        break;
      vector<double> row;
      read_e10(line, 11, 10, row);
      while (row.size() < nte) {
        if (!getline(fin, line))
          break;
        size_t nprev = row.size();
        read_e10(line, 11, 10, row);
        if (row.size() == nprev)
          break;
      }
      if (row.size() > nte)
        row.resize(nte);
      vd_alfi.push_back(row);
    }
  }

  // --- skip PRTF / INREP sections til ALFT ---
  // After PRTI, there are PRTF sections, then INREP, then ALFT
  // Just skip everything until "ALFT( 1)" appears
  while (getline(fin, line)) {
    if (line.find("ALFT(") != string::npos)
      break;
  }
  // Skip the separator line ("---- -------- ...")
  getline(fin, line);

  // Read ALFT data: one line per temperature point
  // Format: (E10.2,1X,20E10.2) — T in cols 0-9, 1X at col 10, ALFT(1..NPRNTI) in cols 11+
  for (size_t i = 0; i < vd_te.size(); i++) {
    getline(fin, line);
    if (!fin)
      break;
    vector<double> vals;
    read_e10(line, 11, ui_nprnti, vals);
    for (size_t j = 0; j < vals.size(); j++)
      vd_alft.push_back(vals[j]);
  }
  if (vd_alft.size() > vd_te.size())
    vd_alft.resize(vd_te.size());

  fin.close();

  // Compute level-summed total
  vd_alf_sum.clear();
  if (!vd_alfi.empty() && !vd_alfi[0].empty()) {
    size_t nte = vd_te.size();
    vd_alf_sum.assign(nte, 0.0);
    for (auto const &row : vd_alfi) {
      for (size_t j = 0; j < row.size() && j < nte; j++) {
        vd_alf_sum[j] += row[j];
      }
    }
  }

  cout << "   read " << vd_alfi.size() << " level records, " << vd_alff.size()
       << " NL-shell records, " << vd_te.size() << " temperature points"
       << (vd_alft.empty() ? "" : ", total DR available") << endl;
}

void ADASadf09::write_ratecoef() {
  if (vd_te.empty()) {
    cerr << "No temperature data available." << endl;
    return;
  }

  double Tmin, Tmax, Tdelta;

  cout << " Give temperature range in K (Tmin, Tmax, delta): ";
  cin >> Tmin >> Tmax >> Tdelta;
  bool logscale_flag = Tmin < 10.0;
  unsigned int npts = (Tmax-Tmin)/Tdelta+1;

  // Build the requested temperature grid
  vector<double> new_te(npts);
  if (logscale_flag) {
    for (unsigned int i = 0; i < npts; i++) {
      new_te[i] = pow(10.0, Tmin + i*Tdelta);
    }
  } else {
    for (unsigned int i = 0; i < npts; i++) {
      new_te[i] = Tmin + i * Tdelta;
    }
  }

  // Interpolate level-summed total and ALFT total to new grid
  vector<double> sum_interp = interpolate(vd_te, vd_alf_sum, new_te);
  vector<double> alft_interp;
  if (!vd_alft.empty()) {
    alft_interp = interpolate(vd_te, vd_alft, new_te);
  }

  // Derive output filename
  string outname = s_filename;
  size_t dot = outname.rfind('.');
  if (dot != string::npos)
    outname.resize(dot);
  outname += "_alpha.dat";

  ofstream fout(outname);
  if (!fout) {
    cerr << "ERROR: cannot open output file " << outname << endl;
    return;
  }

  time_t rawtime;
  time(&rawtime);

  fout << "####################################################################"
          "############"
       << endl;
  fout << "### Total DR rate coefficient (cm^3/s) from ADAS ADF09 file" << endl;
  fout << "###" << endl;
  fout << "###  hydrocal revision     : " << HYDROCAL_REVISION << endl;
  fout << "###               filename : " << outname << endl;
  fout << "###      start date & time : " << asctime(localtime(&rawtime));
  fout << "###             ADF09 file : " << s_filename << endl;
  fout << "###     isoelectronic seq. : " << c3_isoseq << endl;
  fout << "###       nuclear charge Z : " << ui_nuccharge << endl;
  fout << "###               coupling : " << c3_coupling << endl;
  fout << "###        resolved levels : " << ui_nlvl << endl;
  fout << "###     original temp. pts : " << vd_te.size() << endl;
  if (logscale_flag) {
    fout << "###          log10(Tmin/K) : " << Tmin << endl;
    fout << "###          log10(Tmax/K) : " << Tmax << endl;
    fout << "###                Tdelta  : " << Tdelta << endl;
  }
  else {
    fout << "###               Tmin (K) : " << Tmax << endl;
    fout << "###               Tmax (K) : " << Tmax << endl;
    fout << "###             Tdelta (K) : " << Tdelta << endl;
  }
  fout << "###     no. of grid points : " << npts << endl;
  fout << "###               column 1 : Temperature (K)" << endl;
  fout << "###               column 2 : Temperature (eV)" << endl;
  fout << "###               column 3 : DR rate coefficient summed over resolved levels "
          "(cm^3/s)"
       << endl;
  if (!vd_alft.empty()) {
    fout << "###               column 4 : Total DR rate coefficient including unresolved Rydberg levels"
            "table (cm^3/s)"
         << endl;
  }
  fout << "####################################################################"
          "############"
       << endl;
  fout << "# T (K) [1]   kT (eV) [2]   alpha_sum (cm^3/s) [3]";
  if (!vd_alft.empty())
    fout << "   alpha_alft (cm^3/s) [4]";
  fout << endl;
  fout << "###-----------------------------------------------------------------"
          "------------"
       << endl;

  fout << scientific << setprecision(4);
  for (unsigned int i = 0; i < npts; i++) {
    fout << setw(12) << new_te[i] << ", " << new_te[i]*hydroconst::kB_eV_K << ",  " << setw(14) << sum_interp[i];
    if (!vd_alft.empty())
      fout << ",  " << setw(14) << alft_interp[i];
    fout << endl;
  }

  fout.close();
  cout << "   wrote " << npts << " data points to " << outname << endl;
}
