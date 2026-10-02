/**
    @file hydrocalGUI.cxx
    @brief GUI for hydrocal

    @par CREATION
    @author Stefan Schippers
    @date 2021-02-12

    @par VERSION
    @verbatim
    $Id: hydrocalGUI.cxx 2037 2026-07-17 15:21:56Z iamp $
// SPDX-License-Identifier: MIT
    @endverbatim
 */
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cctype>
#include <signal.h>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>

#include <FL/Fl.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Progress.H>
#include <FL/Fl_Text_Buffer.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Window.H>
#include <FL/fl_ask.H>
#include <FL/fl_file_chooser.H>

#include "buildinfo.h"

using namespace std;

Fl_Window *HCwin = nullptr;

const int NO_INPUTS = 30;
Fl_Value_Input *HCinput[NO_INPUTS] = {}; ///< array of input values
double HCvalue[NO_INPUTS];

const int NO_CHECKS = 3;
Fl_Check_Button *HCcheck[NO_CHECKS] = {};

const int NO_STRINGS = 2;
Fl_Input *HCstring_input[NO_STRINGS];

const int NO_BUTTONS = 7;
Fl_Button *HCbutton[NO_BUTTONS] = {};

Fl_Button *HCheader = nullptr, *HCexecute = nullptr, *HCclose = nullptr,
          *HCreset = nullptr, *HChelp = nullptr;
Fl_Text_Display *HCtextdisplay = NULL;
Fl_Text_Buffer *HCtextbuffer = NULL;
Fl_Progress *HCprogress = nullptr;
Fl_Button *HCplot = nullptr, *HCcancel = nullptr;
Fl_Button *HCsave = nullptr, *HCload = nullptr;
std::atomic<pid_t> hydrocal_pid{-1};

enum HCchoices {
  none = -1,
  FranckCondon,
  PIPEMonteCarlo,
  CRYMonteCarlo,
  CSRMonteCarlo,
  ESRMonteCarlo,
  DRalphaMB,
  RRalphaMB
};

int HCchoice = none;

enum StorageRingIDs // copied from Cooler.h
{
  NO_STORAGE_RING,
  CRY_STORAGE_RING, // CRYRING
  CSR_STORAGE_RING,
  ESR_STORAGE_RING,
  CROSSED_BEAMS
};

string hydrocal_cache_dir, hydrocal_infile, hydrocal_outfile,
    hydrocal_emptyfile;
string hydrocal_start_dir;

bool read_data(const char *filename);
bool write_data(int nstrings, int nchecks, int ndata, const char *filename);

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void FranckCondon_callbackC(Fl_Widget * /*w*/, void *) {
  static int i0max = 0;
  if (i0max >= NO_INPUTS)
    i0max = NO_INPUTS - 1;

  // distribution of initial vinrational levels
  if (HCcheck[0]->value()) {
    for (int i = 1; i <= i0max; i++) {
      if ((26 + i) >= NO_INPUTS)
        continue;
      HCinput[26 + i]->hide();
    }
    i0max = 1;
  } else {
    for (int i = 1; i <= HCinput[8]->value(); i++) {
      if ((26 + i) >= NO_INPUTS)
        continue;
      string helpstr = "frac v=" + to_string(i);
      HCinput[26 + i]->copy_label(helpstr.c_str());
      HCinput[26 + i]->show();
    }
    for (int i = HCinput[8]->value() + 1; i <= i0max; i++) {
      if ((26 + i) >= NO_INPUTS)
        continue;
      HCinput[26 + i]->hide();
    }
    i0max = HCinput[8]->value();
  }

  // rotations
  if (HCinput[17]->value() > 0.5) {
    for (int i = 1; i <= 5; i++) {
      HCinput[17 + i]->show();
    }
  } else {
    for (int i = 1; i <= 5; i++) {
      HCinput[17 + i]->hide();
    }
  }

  // error analysis
  if (HCcheck[1]->value()) {
    for (int i = 1; i <= 3; i++) {
      HCinput[23 + i]->show();
    }
  } else {
    for (int i = 1; i <= 3; i++) {
      HCinput[23 + i]->hide();
    }
  }
}

////////////////////////////////////////////////////////////////////////
void FranckCondon_callback(Fl_Widget * /*w*/, void *) {
  HCload->show();
  HCsave->show();
  HCchoice = FranckCondon;
  HCheader->copy_label("hydrocal: Franck-Condon");
  read_data("FranckCondon");
  HCstring_input[0]->show();
  HCstring_input[0]->label("output file (*.fcf)");
  HCcheck[0]->label("Boltzmann also for vibrations");
  HCcheck[0]->show();
  HCcheck[1]->label("Error analysis");
  HCcheck[1]->show();
  HCinput[0]->label("atomic mass of atom 1");
  HCinput[1]->label("atomic mass of atom 2");
  HCinput[2]->label("lower Re (nm)");
  HCinput[3]->label("lower hbar*omega (eV)");
  HCinput[4]->label("lower hbar*omegachi (eV)");
  HCinput[5]->label("upper Re (nm)");
  HCinput[6]->label("upper hbar*omega (eV)");
  HCinput[7]->label("upper hbar*omegachi (eV)");
  HCinput[8]->label("lower max v");
  HCinput[9]->label("upper max v");
  HCinput[10]->label("Lorentzian width (eV)");
  HCinput[11]->label("Gaussian width (eV)");
  HCinput[12]->label("excitation energy (eV)");
  HCinput[13]->label("Emin (eV)");
  HCinput[14]->label("Emax (eV)");
  HCinput[15]->label("Edelta (eV)");
  HCinput[16]->label("temperature (K)");
  HCinput[17]->label("max Jrot");
  HCinput[18]->label("lower Be (1/cm)");
  HCinput[19]->label("lower Lambda");
  HCinput[20]->label("upper Lambda");
  HCinput[21]->label("nuclear spin");
  HCinput[22]->label("symmetry (-1 / +1)");

  HCinput[23]->label("delta Re (nm)");
  HCinput[24]->label("delta hbar*omega (eV)");
  HCinput[25]->label("delta hbar*omegachi (eV)");

  for (int i = 0; i <= 21; i++)
    HCinput[i]->show();

  HCinput[8]->callback(FranckCondon_callbackC);
  HCinput[17]->callback(FranckCondon_callbackC);
  HCcheck[0]->callback(FranckCondon_callbackC);
  HCcheck[1]->callback(FranckCondon_callbackC);
  FranckCondon_callbackC(NULL, NULL);
}

////////////////////////////////////////////////////////////////////////
void FranckCondon_store_data() {
  fstream fout(hydrocal_infile, fstream::out);
  fout << "9" << endl;
  fout << "2" << endl;
  fout << "a" << endl;
  for (int i = 0; i <= 7; i++)
    fout << HCinput[i]->value() << endl;
  int v1max = (int)round(HCinput[8]->value());
  int v2max = (int)round(HCinput[9]->value());
  fout << v1max << endl;
  fout << v2max << endl;
  fout << "0" << endl;
  if (HCcheck[1]->value()) {
    fout << "y" << endl;
    fout << HCinput[24]->value() << endl;
    fout << HCinput[25]->value() << endl;
    fout << HCinput[23]->value() << endl;
  } else {
    fout << "n" << endl;
  }
  fout << "y" << endl;
  fout << HCstring_input[0]->value() << endl;
  fout << HCinput[16]->value() << endl;
  int Jmax = (int)round(HCinput[17]->value());
  fout << Jmax << endl;
  if (Jmax > 0) {
    fout << HCinput[18]->value() << endl;
    fout << HCinput[19]->value() << endl;
    fout << HCinput[20]->value() << endl;
    fout << HCinput[21]->value() << endl;
    if (HCinput[22]->value() < 0.0) {
      fout << "a" << endl;
    } else {
      fout << "s" << endl;
    }
  }
  if (v1max > 0) {
    if (HCcheck[0]->value()) { // temperature for Boltzmann distribution
      fout << "y" << endl;
    } else {
      fout << "n" << endl;
      for (int i = 1; i <= v1max; i++) {
        fout << HCinput[26 + i]->value() << endl;
      }
    }
  }
  for (int i = 10; i <= 12; i++) {
    fout << HCinput[i]->value() << endl;
  }
  fout << HCinput[13]->value() << " " << HCinput[14]->value() << " "
       << HCinput[15]->value() << endl;
  fout << "0" << endl;
  fout << "0" << endl;
  fout.close();

  write_data(1, 2, NO_INPUTS, "FranckCondon");
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void PIPEMonteCarlo_callbackC(Fl_Widget * /*w*/, void *) {
  if (HCcheck[1]->value()) {
    HCstring_input[0]->copy_label("    output file (*.pid)");
  } else {
    HCstring_input[0]->copy_label("output file (*.pmc)");
  }
  if (HCcheck[0]->value()) {
    HCinput[5]->label("core-hole lifetime (fs)");
    HCinput[6]->label("equilibrium distance (eV)");
    HCinput[7]->show();
  } else {
    HCinput[5]->label("min KER (eV)");
    HCinput[6]->label("max KER (eV)");
    HCinput[7]->hide();
  }
}

////////////////////////////////////////////////////////////////////////
void PIPEMonteCarlo_callback(Fl_Widget * /*w*/, void *) {
  HCload->show();
  HCsave->show();
  HCchoice = PIPEMonteCarlo;
  HChelp->hide(); // no help available
  HCheader->copy_label("hydrocal: PIPE Monte-Carlo");
  read_data("PIPEMonteCarlo");
  HCstring_input[0]->show();
  HCstring_input[0]->label("output file (*.pmc)");
  HCcheck[0]->label("KER from core-hole lifetime");
  HCcheck[0]->show();
  HCcheck[1]->label("CST input");
  HCcheck[1]->show();
  HCinput[0]->label("ion energy (keV)");
  HCinput[1]->label("  mass of fragment 1");
  HCinput[2]->label("charge of fragment 1");
  HCinput[3]->label("  mass of fragment 2");
  HCinput[4]->label("charge of fragment 2");
  HCinput[5]->label("min KER (eV)");
  HCinput[6]->label("max KER (eV)");
  HCinput[7]->label("hbar*omega (eV)");
  HCinput[8]->label("anisotropy (-2<= a2 <= 1)");
  HCinput[9]->label("FWHM ion beam (mm)");
  HCinput[10]->label("divergence (mrad)");
  HCinput[11]->label("momentum spread dp/p");
  HCinput[12]->label("x-Trinos (mm)");
  HCinput[13]->label("y-Trinos (mm)");
  HCinput[14]->label("x-Detector (mm)");
  HCinput[15]->label("MC iterations");

  for (int i = 0; i <= 15; i++)
    HCinput[i]->show();
  HCcheck[0]->callback(PIPEMonteCarlo_callbackC);
  HCcheck[1]->callback(PIPEMonteCarlo_callbackC);
  PIPEMonteCarlo_callbackC(NULL, NULL);
}

////////////////////////////////////////////////////////////////////////
void PIPEMonteCarlo_store_data() {
  fstream fout(hydrocal_infile, fstream::out);
  fout << "9" << endl;
  fout << "1" << endl;
  if (HCcheck[1]->value()) {
    fout << "c" << endl;
  } else {
    fout << "f" << endl;
  }
  fout << "n" << endl;
  for (int i = 0; i <= 4; i++)
    fout << HCinput[i]->value() << endl;
  if (HCcheck[0]->value()) {
    fout << "y" << endl;
  } else {
    fout << "n" << endl;
  }
  fout << HCinput[5]->value() << endl;
  fout << HCinput[6]->value() << endl;
  if (HCcheck[0]->value())
    fout << HCinput[7]->value() << endl;
  for (int i = 8; i <= 14; i++) {
    if (i == 8) {
      fout << (int)round(HCinput[i]->value()) << endl;
    } else {
      fout << HCinput[i]->value() << endl;
    }
  }
  if (!HCcheck[1]->value())
    fout << "B" << endl;
  fout << (int)(HCinput[15]->value()) << endl;
  fout << HCstring_input[0]->value() << endl;
  fout << "0" << endl;
  fout << "0" << endl;
  fout.close();

  write_data(1, 2, 16, "PIPEMonteCarlo");
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void CoolerMonteCarlo_callbackC(Fl_Widget * /*w*/, void *) {
  if (HCcheck[2]->value()) { // output of energy distribution
    HCinput[21]->hide();
    HCinput[20]->show();
    HCstring_input[0]->label("output file (*.mcd)");
  } else { // output of convolution
    HCinput[21]->show();
    HCinput[20]->hide();
    HCstring_input[0]->label("output file (*.mcr)");
  }

  string infilename = HCstring_input[1]->value();
  if (infilename == "none") { // RR only
    HCinput[21]->hide();
    HCinput[12]->show();
  } else { // check for autostructure input
    if (!HCcheck[2]->value())
      HCinput[21]->show();
    string ext("bdr");
    size_t found = infilename.rfind(ext);
    if (found != string::npos) {
      // autostructure input file (binned cross sections)
      HCinput[10]->label("Emin (eV)");
      HCinput[11]->label("Emax (eV)");
      HCinput[12]->hide();
    } else {
      // general input file
      HCinput[10]->label("Emin (eV or log)");
      HCinput[11]->label("Emax (eV or log)");
      HCinput[12]->show();
    }
  }
}
////////////////////////////////////////////////////////////////////////
void CoolerMonteCarlo_callback(Fl_Widget * /*w*/, void *srid) {
  HCload->show();
  HCsave->show();
  int storage_ring_ID = *(int *)srid;

  if (storage_ring_ID == CRY_STORAGE_RING) {
    HCchoice = CRYMonteCarlo;
    HCheader->copy_label("hydrocal: CRYRING Monte-Carlo");
    read_data("CRYMonteCarlo");
  } else if (storage_ring_ID == CSR_STORAGE_RING) {
    HCchoice = CSRMonteCarlo;
    HCheader->copy_label("hydrocal: CSR Monte-Carlo");
    read_data("CSRMonteCarlo");
  } else if (storage_ring_ID == ESR_STORAGE_RING) {
    HCchoice = ESRMonteCarlo;
    HCheader->copy_label("hydrocal: ESR Monte-Carlo");
    read_data("ESRMonteCarlo");
  } else {
    return;
  }
  HCstring_input[1]->show();
  HCstring_input[1]->label("theory input file (*.*, use 'none' for RR only)");
  HCstring_input[1]->callback(CoolerMonteCarlo_callbackC);
  HCstring_input[0]->show();
  HCstring_input[0]->label("output file (*.mcr)");
  HCcheck[0]->label("scanning with drift tube");
  if (storage_ring_ID == CRY_STORAGE_RING) {
    HCcheck[0]->show();
  } else {
    HCcheck[0]->hide();
  }
  HCcheck[1]->label("electrons slower than ions");
  HCcheck[1]->show();
  HCcheck[2]->label("output of energy distribution");
  HCcheck[2]->show();
  HCcheck[2]->callback(CoolerMonteCarlo_callbackC);
  HCinput[0]->label("expansion factor");
  HCinput[1]->label("electron current (mA)");
  HCinput[2]->label("cathode voltage (V)");
  HCinput[3]->label("cooling energy (eV)");
  HCinput[4]->label("ion mass (u)");
  HCinput[5]->label("ion charge");
  HCinput[6]->label("ion momentum spread (dp/p)");
  HCinput[7]->label("ion kT perp (meV)");
  HCinput[8]->label("electron kT par (meV)");
  HCinput[9]->label("electron kT perp (meV)");
  HCinput[10]->label("Emin (eV or log)");
  HCinput[11]->label("Emax (eV or log)");
  HCinput[12]->label("Edelta (eV or log)");
  HCinput[13]->label("Emax for output (eV)");
  HCinput[14]->label("RR nmin");
  HCinput[15]->label("RR lmin");
  HCinput[16]->label("RR nmax");
  HCinput[17]->label("RR nele in (nmin,lmin)");
  HCinput[18]->label("level number");
  HCinput[19]->label("MC iterations");
  HCinput[20]->label("relative energy (eV)");
  HCinput[21]->label("DR energy shift (eV)");

  for (int i = 0; i <= 21; i++)
    HCinput[i]->show();
  CoolerMonteCarlo_callbackC(NULL, NULL);
}

////////////////////////////////////////////////////////////////////////
void CoolerMonteCarlo_store_data(int storage_ring_ID) {
  bool RRonly_flag = (string(HCstring_input[1]->value()) == "none");
  fstream fout(hydrocal_infile, fstream::out);
  fout << "6" << endl;
  fout << "2" << endl;
  fout << storage_ring_ID << endl;
  fout << HCinput[2]->value() << endl; // cathode voltage
  if (storage_ring_ID != ESR_STORAGE_RING) {
    fout << HCinput[0]->value()
         << endl; // expansion factor (ESRinit hardcodes it)
  }
  fout << HCinput[1]->value() << endl; // electron current (mA)
  if (storage_ring_ID == CRY_STORAGE_RING) {
    if (HCcheck[0]->value()) {
      fout << "d" << endl; // scanning with drifttube
    } else {
      fout << "c" << endl; // scanning with cathode
    }
  }
  fout << "y" << endl; // account for variation of drifttube potential
  if (storage_ring_ID != CSR_STORAGE_RING) {
    fout << "y"
         << endl; // account for angular variation (CSRinit doesn't read this)
  }
  if (storage_ring_ID == ESR_STORAGE_RING) {
    fout << "0.48" << endl; // offset angle (mrad)
  }
  for (int i = 3; i <= 7; i++) {
    fout << HCinput[i]->value() << endl;
  }
  fout << HCinput[8]->value() << " " << HCinput[9]->value()
       << endl; // ebeam temperatures

  fout << HCstring_input[1]->value() << endl; // DR theory file
  if (!RRonly_flag) {
    fout << HCinput[21]->value() << endl; // overall DR energy shift
    fout << HCinput[18]->value() << endl; // level number
    if (!HCinput[12]->visible()) { // autostructure input file (no Edelta shown)
      fout << "y" << endl;         // extend energy range
    }
  }

  fout << HCinput[10]->value() << " " << HCinput[11]->value() << endl;
  if (HCinput[12]->visible()) { // general input file
    fout << HCinput[12]->value() << endl;
  }
  fout << HCinput[13]->value() << endl; // max energy for output
  fout << HCinput[14]->value() << " " << HCinput[15]->value() << " "
       << HCinput[16]->value() << endl;
  fout << HCinput[17]->value() << endl;
  if (HCcheck[1]->value()) {
    fout << "s" << endl;
  } else {
    fout << "f" << endl;
  }
  if (HCcheck[2]->value()) {
    fout << "d" << endl;
    fout << HCinput[20]->value() << endl;
  } else {
    fout << "r" << endl;
  }
  fout << (int)(HCinput[19]->value()) << endl;
  fout << HCstring_input[0]->value() << endl;
  fout << "e" << endl;
  fout.close();
  if (storage_ring_ID == CRY_STORAGE_RING) {
    write_data(2, 3, 21, "CRYMonteCarlo");
  } else if (storage_ring_ID == CSR_STORAGE_RING) {
    write_data(2, 3, 21, "CSRMonteCarlo");
  } else if (storage_ring_ID == ESR_STORAGE_RING) {
    write_data(2, 3, 21, "ESRMonteCarlo");
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void DRalphaMB_callbackC(Fl_Widget * /*w*/, void *) {
  bool autos_flag = false;
  string extension = HCstring_input[1]->value();
  extension.erase(0, extension.find('.'));
  if (extension == ".bdr") {
    HCstring_input[0]->label("output file (*.cdr)");
    autos_flag = true;
  } else if (extension == ".brr") {
    HCstring_input[0]->label("output file (*.crr)");
    autos_flag = true;
  }

  if (autos_flag) {
    HCcheck[0]->show();
    if (HCcheck[0]->value()) {
      HCinput[5]->show();
    } else {
      HCinput[5]->hide();
    }
  } else {
    HCcheck[0]->hide();
    HCinput[5]->hide();
  }
}
////////////////////////////////////////////////////////////////////////
void DRalphaMB_callback(Fl_Widget * /*w*/, void *) {
  HCload->show();
  HCsave->show();
  HCchoice = DRalphaMB;
  HCheader->copy_label("hydrocal: Merged-beams DR rate coefficient");
  read_data("DRalphaMB");
  HCcheck[0]->label("more than 1 initial level");
  HCcheck[0]->callback(DRalphaMB_callbackC);
  HCstring_input[0]->show();
  HCstring_input[0]->label("output file (*.cdr)");
  HCstring_input[1]->show();
  HCstring_input[1]->label("DR theory file (*.*)");
  HCstring_input[1]->callback(DRalphaMB_callbackC);

  HCinput[0]->label("Emin (eV or log)");
  HCinput[1]->label("Emax (eV or log)");
  HCinput[2]->label("Edelta (eV or log)");
  HCinput[3]->label("electron kT par (meV)");
  HCinput[4]->label("electron kT perp (meV)");
  HCinput[5]->label("initial level number");

  for (int i = 0; i <= 5; i++)
    HCinput[i]->show();
  DRalphaMB_callbackC(NULL, NULL);
}

////////////////////////////////////////////////////////////////////////
void DRalphaMB_store_data() {
  string extension = HCstring_input[1]->value();
  size_t dot = extension.find_last_of('.');
  if (dot != string::npos)
    extension = extension.substr(dot);
  bool autos_flag = (extension == ".bdr" || extension == ".brr");

  fstream fout(hydrocal_infile, fstream::out);
  fout << "6" << endl;
  fout << "1" << endl;
  fout << "5" << endl;
  fout << HCstring_input[1]->value() << endl;
  if (autos_flag) {
    if (HCcheck[0]->value())
      fout << HCinput[5]->value() << endl;
  } else {
    fout << "1" << endl;
    fout << "1" << endl;
  }
  fout << "0" << endl;
  fout << HCinput[0]->value() << " " << HCinput[1]->value() << " "
       << HCinput[2]->value() << endl;
  fout << HCinput[3]->value() << " " << HCinput[4]->value() << endl;
  fout << "y" << endl;
  fout << HCstring_input[0]->value() << endl;
  fout << "n" << endl;
  fout << "0" << endl;
  fout << "0" << endl;
  fout.close();

  write_data(2, 1, 6, "DRalphaMB");
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void RRalphaMB_callbackC(Fl_Widget * /*w*/, void *) {
  if (HCcheck[0]->value()) {
    HCstring_input[1]->show();
  } else {
    HCstring_input[1]->hide();
  }
}

////////////////////////////////////////////////////////////////////////
void RRalphaMB_callback(Fl_Widget * /*w*/, void *) {
  HCload->show();
  HCsave->show();

  HCchoice = RRalphaMB;
  HCheader->copy_label("hydrocal: Merged-beams RR rate coefficient");
  read_data("RRalphaMB");
  HCstring_input[0]->show();
  HCstring_input[0]->label("output file (*.crr)");
  HCstring_input[1]->label("fraction table (*.fnl)");

  HCcheck[0]->label("surviving fractions from file");
  HCcheck[0]->show();
  HCcheck[0]->callback(RRalphaMB_callbackC);

  HCinput[0]->label("RR Zeff");
  HCinput[1]->label("RR nmin");
  HCinput[2]->label("RR lmin");
  HCinput[3]->label("RR nmax");
  HCinput[4]->label("RR nele in (nmin,lmin)");
  HCinput[5]->label("Emin (eV or log)");
  HCinput[6]->label("Emax (eV or log)");
  HCinput[7]->label("Edelta (eV or log)");
  HCinput[8]->label("electron kT par (meV)");
  HCinput[9]->label("electron kT perp (meV)");

  for (int i = 0; i <= 9; i++)
    HCinput[i]->show();
  RRalphaMB_callbackC(NULL, NULL);
}

////////////////////////////////////////////////////////////////////////
void RRalphaMB_store_data() {
  fstream fout(hydrocal_infile, fstream::out);
  fout << "6" << endl;
  fout << "1" << endl;
  fout << "1" << endl;
  fout << HCinput[0]->value() << " " << HCinput[1]->value() << " "
       << HCinput[2]->value() << " " << HCinput[3]->value() << endl;
  fout << HCinput[4]->value() << endl;
  if (HCcheck[0]->value()) {
    fout << "y" << endl;
    fout << HCstring_input[1]->value() << endl;
  } else {
    fout << "n" << endl;
  }
  fout << "0" << endl;
  fout << HCinput[5]->value() << " " << HCinput[6]->value() << " "
       << HCinput[7]->value() << endl;
  fout << HCinput[8]->value() << " " << HCinput[9]->value() << endl;
  fout << "y" << endl;
  fout << HCstring_input[0]->value() << endl;
  fout << "n" << endl;
  fout << "0" << endl;
  fout << "0" << endl;
  fout.close();

  write_data(2, 1, 10, "RRalphaMB");
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
bool read_data(const char *filename) {
  for (int i = 0; i < NO_STRINGS; i++) {
    HCstring_input[i]->hide();
    HCstring_input[i]->value("");
  }
  for (int i = 0; i < NO_CHECKS; i++) {
    HCcheck[i]->hide();
    HCcheck[i]->value(0);
  }
  for (int i = 0; i < NO_INPUTS; i++) {
    HCinput[i]->hide();
    HCinput[i]->value(0);
  }
  for (int i = 0; i < NO_BUTTONS; i++)
    if (i != HCchoice)
      HCbutton[i]->deactivate();

  string fn = hydrocal_cache_dir + filename;
  fstream fin(fn, fstream::in);
  if (!fin.is_open())
    return 0;

  int nstrings = 0, nchecks = 0, ndata = 0;
  fin >> nstrings;
  if (!fin.good())
    return false;
  fin >> nchecks;
  if (!fin.good())
    return false;
  fin >> ndata;
  if (!fin.good())
    return false;
  double value;
  string helpstr;
  for (int i = 0; i < nstrings; i++) {
    fin >> helpstr;
    if (!fin.good())
      return false;
    HCstring_input[i]->value(helpstr.data());
  }
  for (int i = 0; i < nchecks; i++) {
    if (i >= NO_CHECKS)
      break;
    fin >> value;
    if (!fin.good())
      return false;
    HCcheck[i]->value(value);
  }
  for (int i = 0; i < ndata; i++) {
    if (i >= NO_INPUTS)
      break;
    fin >> value;
    if (!fin.good())
      return false;
    HCinput[i]->value(value);
  }
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////
bool write_data(int nstrings, int nchecks, int ndata, const char *filename) {
  string fn = hydrocal_cache_dir + filename;
  fstream fout(fn, fstream::out);
  if (!fout.is_open())
    return false;
  fout << nstrings << endl;
  fout << nchecks << endl;
  int wdata = ndata;
  if (wdata > NO_INPUTS)
    wdata = NO_INPUTS;
  fout << wdata << endl;
  for (int i = 0; i < nstrings; i++)
    fout << HCstring_input[i]->value() << endl;
  for (int i = 0; i < nchecks; i++) {
    int checkvalue = HCcheck[i]->value() ? 1 : 0;
    fout << checkvalue << endl;
  }
  for (int i = 0; i < wdata; i++)
    fout << HCinput[i]->value() << endl;
  return true;
}

////////////////////////////////////////////////////////////////////////
const char *get_mode_name() {
  switch (HCchoice) {
  case FranckCondon:
    return "FranckCondon";
  case PIPEMonteCarlo:
    return "PIPEMonteCarlo";
  case CRYMonteCarlo:
    return "CRYMonteCarlo";
  case CSRMonteCarlo:
    return "CSRMonteCarlo";
  case ESRMonteCarlo:
    return "ESRMonteCarlo";
  case DRalphaMB:
    return "DRalphaMB";
  case RRalphaMB:
    return "RRalphaMB";
  default:
    return "none";
  }
}

///////////////////////////////////////////////////////////////////////////////////////
bool write_values_to_file(const char *fullpath) {
  fstream fout(fullpath, fstream::out);
  if (!fout.is_open())
    return false;
  fout << get_mode_name() << endl;
  fout << NO_STRINGS << endl;
  fout << NO_CHECKS << endl;
  fout << NO_INPUTS << endl;
  for (int i = 0; i < NO_STRINGS; i++)
    fout << HCstring_input[i]->value() << endl;
  for (int i = 0; i < NO_CHECKS; i++) {
    int checkvalue = HCcheck[i]->value() ? 1 : 0;
    fout << checkvalue << endl;
  }
  for (int i = 0; i < NO_INPUTS; i++)
    fout << HCinput[i]->value() << endl;
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////
bool read_values_from_file(const char *fullpath) {
  fstream fin(fullpath, fstream::in);
  if (!fin.is_open())
    return false;
  string saved_mode;
  fin >> saved_mode;
  if (!fin.good())
    return false;
  if (saved_mode != get_mode_name()) {
    fl_alert("Preset file is for mode '%s', not '%s'.", saved_mode.c_str(),
             get_mode_name());
    return false;
  }
  int nstrings = 0, nchecks = 0, ndata = 0;
  fin >> nstrings >> nchecks >> ndata;
  if (!fin.good())
    return false;
  string helpstr;
  double value;
  for (int i = 0; i < nstrings && i < NO_STRINGS; i++) {
    fin >> helpstr;
    HCstring_input[i]->value(helpstr.data());
  }
  for (int i = 0; i < nchecks && i < NO_CHECKS; i++) {
    fin >> value;
    HCcheck[i]->value(value);
  }
  for (int i = 0; i < ndata && i < NO_INPUTS; i++) {
    fin >> value;
    HCinput[i]->value(value);
  }
  return true;
}

////////////////////////////////////////////////////////////////////////
void refresh_mode_ui() {
  switch (HCchoice) {
  case FranckCondon:
    FranckCondon_callbackC(NULL, NULL);
    break;
  case PIPEMonteCarlo:
    PIPEMonteCarlo_callbackC(NULL, NULL);
    break;
  case CRYMonteCarlo:
  case CSRMonteCarlo:
  case ESRMonteCarlo:
    CoolerMonteCarlo_callbackC(NULL, NULL);
    break;
  case DRalphaMB:
    DRalphaMB_callbackC(NULL, NULL);
    break;
  case RRalphaMB:
    RRalphaMB_callbackC(NULL, NULL);
    break;
  default:
    break;
  }
}

////////////////////////////////////////////////////////////////////////
void save_callback(Fl_Widget *, void *) {
  if (HCchoice == none) {
    fl_alert("Please select a mode first.");
    return;
  }
  const char *fn =
      fl_file_chooser("Save settings", "*.hcg", hydrocal_start_dir.c_str());
  if (!fn)
    return;
  if (!write_values_to_file(fn)) {
    fl_alert("Could not save settings to '%s'.", fn);
    return;
  }
  string saved(fn);
  size_t pos = saved.rfind('/');
  if (pos != string::npos)
    hydrocal_start_dir = saved.substr(0, pos + 1);
}

////////////////////////////////////////////////////////////////////////
void load_callback(Fl_Widget *, void *) {
  if (HCchoice == none) {
    fl_alert("Please select a mode first.");
    return;
  }
  const char *fn =
      fl_file_chooser("Load settings", "*.hcg", hydrocal_start_dir.c_str());
  if (!fn)
    return;
  if (!read_values_from_file(fn)) {
    fl_alert("Could not load settings from '%s'.", fn);
    return;
  }
  string loaded(fn);
  size_t pos = loaded.rfind('/');
  if (pos != string::npos)
    hydrocal_start_dir = loaded.substr(0, pos + 1);
  refresh_mode_ui();
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief call back function from reset button
 */
void reset_callback(Fl_Widget * /*w*/, void *) {
  HCload->hide();
  HCsave->hide();
  HCchoice = none;
  HChelp->show();
  HCheader->copy_label("hydrocal GUI");
  for (int i = 0; i < NO_STRINGS; i++) {
    HCstring_input[i]->hide();
    HCstring_input[i]->callback();
  }
  for (int i = 0; i < NO_CHECKS; i++) {
    HCcheck[i]->hide();
    HCcheck[i]->callback();
  }
  for (int i = 0; i < NO_INPUTS; i++) {
    HCinput[i]->hide();
    HCinput[i]->callback();
  }
  for (int i = 0; i < NO_BUTTONS; i++)
    HCbutton[i]->activate();
  HCplot->hide();
  HCtextbuffer->loadfile(hydrocal_emptyfile.c_str());
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief Sets the text display content (thread-safe via Fl::awake)
 */
static void hc_set_text(void *data) {
  string *s = static_cast<string *>(data);
  HCtextbuffer->text(s->c_str());
  delete s;
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief Appends text to the display (thread-safe via Fl::awake)
 */
static void hc_append_text(void *data) {
  string *s = static_cast<string *>(data);
  HCtextbuffer->append(s->c_str());
  delete s;
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief Animates the progress bar (indeterminate pulse)
 */
static void progress_anim_cb(void *) {
  if (!HCprogress->visible())
    return;
  static double val = 0;
  static int dir = 1;
  val += dir * 6;
  if (val >= 94)
    dir = -1;
  if (val <= 6)
    dir = 1;
  HCprogress->value(val);
  HCprogress->label("running ...");
  Fl::repeat_timeout(0.08, progress_anim_cb);
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief Deactivates GUI buttons and starts progress bar during calculation
 */
static void hc_deactivate_buttons(void *) {
  HCexecute->deactivate();
  HCclose->deactivate();
  HCreset->deactivate();
  HCplot->hide();
  HCcancel->show();
  HCcancel->activate();
  HCcancel->label("cancel");
  HCprogress->show();
  HCprogress->value(0);
  Fl::add_timeout(0.08, progress_anim_cb);
  HCwin->redraw();
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief Reactivates GUI buttons and stops progress bar after calculation
 */
static void hc_reactivate_buttons(void *) {
  HCexecute->activate();
  HCclose->activate();
  HCreset->activate();
  HCplot->show();
  HCcancel->hide();
  HCprogress->hide();
  HCwin->redraw();
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief Cancels the running hydrocal calculation
 */
void cancel_callback(Fl_Widget *, void *) {
  if (!HCcancel)
    return;
  pid_t pid = hydrocal_pid;
  if (pid > 0) {
    kill(pid, SIGTERM);
    HCcancel->deactivate();
    HCcancel->label("cancelling...");
    HCwin->redraw();
  }
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief Executes hydrocal in a detached thread, piping output to the GUI
 */
void run_hydrocal() {
  int pipefd[2];
  if (pipe(pipefd) == -1) {
    Fl::awake(hc_set_text, new string("ERROR: Could not create pipe!\n"));
    Fl::awake(hc_reactivate_buttons);
    return;
  }

  pid_t pid = fork();
  if (pid == -1) {
    close(pipefd[0]);
    close(pipefd[1]);
    Fl::awake(hc_set_text, new string("ERROR: Could not fork!\n"));
    Fl::awake(hc_reactivate_buttons);
    return;
  }

  if (pid == 0) {
    // child: redirect stdin/out/err and exec hydrocal
    close(pipefd[0]);
    int fd = open(hydrocal_infile.c_str(), O_RDONLY);
    if (fd != -1) {
      dup2(fd, STDIN_FILENO);
      close(fd);
    }
    dup2(pipefd[1], STDOUT_FILENO);
    dup2(pipefd[1], STDERR_FILENO);
    close(pipefd[1]);
    chdir(hydrocal_start_dir.c_str());
    execlp("hydrocal", "hydrocal", NULL);
    _exit(1);
  }

  // parent: set PID before revealing the cancel button
  close(pipefd[1]);
  hydrocal_pid = pid;
  Fl::awake(hc_set_text, new string("starting up hydrocal ...\n\n"));
  Fl::awake(hc_deactivate_buttons);

  FILE *pipe = fdopen(pipefd[0], "r");
  if (!pipe) {
    close(pipefd[0]);
    hydrocal_pid = 0;
    Fl::awake(hc_set_text, new string("ERROR: fdopen failed!\n"));
    Fl::awake(hc_reactivate_buttons);
    return;
  }
  char buffer[4096];
  while (fgets(buffer, sizeof(buffer), pipe)) {
    Fl::awake(hc_append_text, new string(buffer));
  }
  fclose(pipe);

  int status;
  waitpid(pid, &status, 0);
  hydrocal_pid = 0;

  if (WIFSIGNALED(status) && WTERMSIG(status) == SIGTERM) {
    Fl::awake(hc_append_text, new string("\n\nhydrocal has been cancelled!\n"));
  } else {
    Fl::awake(hc_append_text, new string("\n\nhydrocal has terminated!\n"));
  }
  Fl::awake(hc_reactivate_buttons);
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief Returns the expected output filename for the current mode
 */
static string get_output_filename() {
  string base = HCstring_input[0]->value();
  if (base.empty())
    return "";
  switch (HCchoice) {
  case CRYMonteCarlo:
  case CSRMonteCarlo:
  case ESRMonteCarlo: {
    string ext = HCcheck[2]->value() ? ".mcd" : ".mcr";
    if (base.size() > ext.size() &&
        base.substr(base.size() - ext.size()) == ext)
      return base;
    return base + ext;
  }
  default:
    return base;
  }
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief Launches gnuplot to visualize the calculation output
 */
void plot_callback(Fl_Widget *, void *) {
  string file = get_output_filename();
  if (file.empty() || HCchoice == none) {
    fl_alert("No output file to plot. Run a calculation first.");
    return;
  }

  ifstream ftest(file);
  if (!ftest.good()) {
    fl_alert("Output file '%s' not found.\nCheck that the calculation finished "
             "successfully.",
             file.c_str());
    return;
  }
  ftest.close();

  FILE *gp = popen("gnuplot -persist", "w");
  if (!gp) {
    fl_alert("Could not start gnuplot. Is it installed?");
    return;
  }

  fprintf(gp, "set datafile separator ','\n");
  switch (HCchoice) {
  case FranckCondon:
    fprintf(gp, "set xlabel 'Energy (eV)'\n");
    fprintf(gp, "set ylabel 'FC factor'\n");
    fprintf(gp, "plot '%s' u 1:2 w l title 'Franck-Condon'\n", file.c_str());
    break;
  case CRYMonteCarlo:
  case CSRMonteCarlo:
  case ESRMonteCarlo:
    if (HCcheck[2]->value()) {
      fprintf(gp, "set xlabel 'CM energy (eV)'\n");
      fprintf(gp, "set ylabel 'distribution (1/eV)'\n");
      fprintf(gp, "plot '%s' u 1:2 w l title 'energy distribution'\n",
              file.c_str());
    } else {
      bool RRonly = (string(HCstring_input[1]->value()) == "none");
      fprintf(gp, "set xlabel 'CM energy (eV)'\n");
      fprintf(gp, "set ylabel 'alpha (cm^3/s)'\n");
      fprintf(gp, "set logscale y\n");
      if (RRonly) {
        fprintf(gp,
                "plot '%s' u 2:4 w l title 'RR', '' u 2:5 w l title 'total'\n",
                file.c_str());
      } else {
        fprintf(gp,
                "plot '%s' u 2:3 w l title 'DR', '' u 2:4 w l title 'RR', '' u "
                "2:5 w l title 'total'\n",
                file.c_str());
      }
    }
    break;
  case DRalphaMB:
  case RRalphaMB:
    fprintf(gp, "set xlabel 'Energy (eV)'\n");
    fprintf(gp, "set ylabel 'alpha (cm^3/s)'\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "plot '%s' u 1:2 w l title 'rate coefficient'\n", file.c_str());
    break;
  case PIPEMonteCarlo:
    fprintf(gp, "set xlabel 'position (mm)'\n");
    fprintf(gp, "set ylabel 'counts'\n");
    fprintf(gp, "plot '%s' u 1:2 w l title 'PIPE MC'\n", file.c_str());
    break;
  default:
    fprintf(gp, "plot '%s' u 1:2 w l title 'data'\n", file.c_str());
  }
  pclose(gp);
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief Callback with is executed whenever the text buffer is modified
 */
void text_buffer_modified_callback(int /*pos*/, int /*nInserted*/,
                                   int /*nDeleted*/, int /*nRestyled*/,
                                   const char * /*deletedText*/,
                                   void * /*cbArg*/) {
  // cout << "text buffer modified at position " << pos << ", " << nInserted <<
  // " characters inserted, " << nDeleted << " characters deleted." << endl;
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief Calls the hydrocal code
 */
void execute_hydrocal_callback(Fl_Widget *, void *) {
  switch (HCchoice) {
  case FranckCondon:
    FranckCondon_store_data();
    break;
  case PIPEMonteCarlo:
    PIPEMonteCarlo_store_data();
    break;
  case CRYMonteCarlo:
    CoolerMonteCarlo_store_data(CRY_STORAGE_RING);
    break;
  case CSRMonteCarlo:
    CoolerMonteCarlo_store_data(CSR_STORAGE_RING);
    break;
  case ESRMonteCarlo:
    CoolerMonteCarlo_store_data(ESR_STORAGE_RING);
    break;
  case RRalphaMB:
    RRalphaMB_store_data();
    break;
  case DRalphaMB:
    DRalphaMB_store_data();
    break;
  default:
    return;
  }

  thread trh(run_hydrocal);
  trh.detach();
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief brings up the help window
 */
void help_callback(Fl_Widget *, void *) {
  switch (HCchoice) {
  case none:
    fl_message("Click on one of the blue buttons to proceed.\n \nReturn to top "
               "menu by clicking on 'reset'.");
    break;
  case ESRMonteCarlo:
  case CSRMonteCarlo:
  case CRYMonteCarlo:
    fl_message(
        "theory input file:\n*.bdr: DR input file (binned cross sections) from "
        "AUTOSTRUCTRUE\n*.*  : "
        "individual DR resonances, min 2, max 4 columns (energy strength width "
        "...)\n\nLinear energy range: "
        "e.g. Emin=0 Emax=200 Edelta=0.5    ->    0 eV to 200 eV in 0.5 eV "
        "steps\n\nLog energy range: Emin=-5 "
        "Emax=2 Edelta=0.1    ->    1E-5 eV to 1E2 eV 10 points per decade\n");
    break;
  case FranckCondon:
    fl_message("Boltzmann also for vibrations: Rotational temperature is also "
               "used for populations of vibrational "
               "levels.\n                         If not checked: Individual "
               "input of vibrational level "
               "populations.\n\n \nmax Nrot = 0: no rotational broadening");
    break;
  case DRalphaMB:
    fl_message("theory input file:\n*.bdr: DR input file (binned cross "
               "sections) from AUTOSTRUCTRUE\n*.brr: RR input file "
               "(binned cross sections) from AUTOSTRUCTURE\n*.*  : individual "
               "DR resonances, min 2, max 4 columns (energy "
               "strength width ...) \n\nLinear energy range: e.g. Emin=0 "
               "Emax=200 Edelta=0.5    ->    0 eV to 200 eV in 0.5 "
               "eV steps\n\nLog energy range: Emin=-5 Emax=2 Edelta=0.1    ->  "
               "  1E-5 eV to 1E2 eV 10 points per decade\n");
    break;
  case RRalphaMB:
    fl_message("Linear energy range: e.g. Emin=0 Emax=200 Edelta=0.5    ->    "
               "0 eV to 200 eV in 0.5 eV steps\n\nLog "
               "energy range: Emin=-5 Emax=2 Edelta=0.1    ->    1E-5 eV to "
               "1E2 eV 10 points per decade\n");
    break;
  }
}

////////////////////////////////////////////////////////////////////////
/**
 * @brief Quits the hydrocal GUI program
 */
void gui_close_callback(Fl_Widget *, void *) {
  pid_t pid = hydrocal_pid.load();
  if (pid > 0) {
    kill(pid, SIGTERM);
    int status;
    waitpid(pid, &status, 0);
  }
  exit(1);
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
/**
 * @brief Main of hydrocal GUI
 *
 * @param argc Number of command line parameters
 * @param argv Array of command line parameter strings
 */
int main(int /*argc*/, char * /*argv*/[]) {
  double w = 1700;
  double h = 830;

  // make sure that the storage directory $(HOME)/.hydrocalgui/ exists
  // this directory is used for storing the inputs to the various menu items
  char *home = getenv("HOME");
  if (!home) {
    cout << "ERROR: Environment variable $HOME not defined" << endl;
    exit(1);
  }
  hydrocal_cache_dir = string(home) + "/.hydrocalgui/";
  mkdir(hydrocal_cache_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd)))
    hydrocal_start_dir = string(cwd) + "/";
  else
    hydrocal_start_dir = "./";
  hydrocal_infile = hydrocal_cache_dir + "hcin";
  hydrocal_outfile = hydrocal_cache_dir + "hcout";
  hydrocal_emptyfile = hydrocal_cache_dir + "empty";
  fstream fempty(hydrocal_emptyfile, fstream::out);
  fempty << " ";
  fempty.close();

  Fl::lock(); // enable lock mechanism for multi-threading
  //**************************************************************************
  if (!HCwin) {
    HCwin = new Fl_Double_Window(int(w), int(h));
    string label = "hydrocal GUI - revision " + string(HYDROCAL_REVISION); // from buildinfo.h included above
    if (!isdigit(label.back())) label.resize(label.size()-1); // cut off the last charcater if it is not numeric
    HCwin->copy_label(label.data());
    HCwin->iconlabel("HC");
    Fl_Group::current()->resizable(HCwin);

    int x_ref = 250, y_ref = 100;

    int width = 100, height = 22;

    int xx_ref = x_ref, yy_ref = 0;

    HCstring_input[1] =
        new Fl_Input(xx_ref + 4 * width, y_ref, 3 * width, height);
    HCstring_input[1]->hide();

    y_ref += 2 * height;

    for (int i = 0; i < NO_CHECKS; i++) {
      y_ref += (i / 2) * height;
      HCcheck[i] = new Fl_Check_Button(xx_ref - 200 + (i % 2 + 1) * 300, y_ref,
                                       2.5 * width, height, "check");
      HCcheck[i]->value(0);
      HCcheck[i]->hide();
    }

    y_ref += 2 * height;

    for (int i = 0; i < NO_INPUTS; i++) {
      if ((i % 16) == 0) {
        yy_ref = y_ref;
        xx_ref += 300;
      }
      HCinput[i] = new Fl_Value_Input(xx_ref, yy_ref, width, height);
      HCinput[i]->value(0);
      HCinput[i]->hide();
      yy_ref += 1.5 * height;
    }

    y_ref += 25 * height;
    HCstring_input[0] =
        new Fl_Input(x_ref + 4 * width, y_ref, 3 * width, height);
    HCstring_input[0]->hide();

    HCchoice = -1;

    //*****************************************************
    HCheader = new Fl_Button(500, 20, 700, 50, "hydrocal GUI");
    HCheader->color(FL_GRAY);
    HCheader->box(FL_FLAT_BOX);
    HCheader->labelsize(24);
    HCheader->labelfont(FL_HELVETICA_BOLD);
    HCheader->labelcolor(FL_DARK_RED);

    x_ref = 30;
    y_ref = 30;
    int y_delta = 50;

    HCreset = new Fl_Button(x_ref, y_ref, 120, 30, "reset");
    HCreset->callback(reset_callback);
    HCreset->color(FL_GRAY);
    HCreset->box(FL_UP_BOX);
    HCreset->labelsize(20);
    HCreset->labelfont(FL_HELVETICA);
    HCreset->labelcolor(FL_BLACK);

    HChelp = new Fl_Button(x_ref + 140, y_ref, 120, 30, "help");
    HChelp->callback(help_callback);
    HChelp->color(FL_YELLOW);
    HChelp->box(FL_UP_BOX);
    HChelp->labelsize(20);
    HChelp->labelfont(FL_HELVETICA);
    HChelp->labelcolor(FL_BLACK);

    HCload = new Fl_Button(x_ref + 300, y_ref, 100, 30, "load");
    HCload->callback(load_callback);
    HCload->color(FL_GRAY);
    HCload->box(FL_UP_BOX);
    HCload->labelsize(18);
    HCload->labelfont(FL_HELVETICA);
    HCload->labelcolor(FL_BLACK);
    HCload->hide();

    HCsave = new Fl_Button(x_ref + 410, y_ref, 100, 30, "save");
    HCsave->callback(save_callback);
    HCsave->color(FL_GRAY);
    HCsave->box(FL_UP_BOX);
    HCsave->labelsize(18);
    HCsave->labelfont(FL_HELVETICA);
    HCsave->labelcolor(FL_BLACK);
    HCsave->hide();

    y_ref += y_delta;
    width = 260;
    height = 30;

    HCbutton[0] = new Fl_Button(x_ref, y_ref, width, height, "Franck-Condon");
    HCbutton[0]->callback(FranckCondon_callback);

    y_ref += y_delta;
    HCbutton[1] =
        new Fl_Button(x_ref, y_ref, width, height, "PIPE Monte-Carlo");
    HCbutton[1]->callback(PIPEMonteCarlo_callback);

    y_ref += y_delta;
    const int CRYid = CRY_STORAGE_RING;
    HCbutton[2] =
        new Fl_Button(x_ref, y_ref, width, height, "CRYRING Monte-Carlo");
    HCbutton[2]->callback(CoolerMonteCarlo_callback, (void *)(&CRYid));

    y_ref += y_delta;
    const int CSRid = CSR_STORAGE_RING;
    HCbutton[3] = new Fl_Button(x_ref, y_ref, width, height, "CSR Monte-Carlo");
    HCbutton[3]->callback(CoolerMonteCarlo_callback, (void *)(&CSRid));

    y_ref += y_delta;
    const int ESRid = ESR_STORAGE_RING;
    HCbutton[4] = new Fl_Button(x_ref, y_ref, width, height, "ESR Monte-Carlo");
    HCbutton[4]->callback(CoolerMonteCarlo_callback, (void *)(&ESRid));

    y_ref += y_delta;
    HCbutton[5] =
        new Fl_Button(x_ref, y_ref, width, height, "DR rate coefficient (MB)");
    HCbutton[5]->callback(DRalphaMB_callback);

    y_ref += y_delta;
    HCbutton[6] =
        new Fl_Button(x_ref, y_ref, width, height, "RR rate coefficient (MB)");
    HCbutton[6]->callback(RRalphaMB_callback);

    for (int i = 0; i < NO_BUTTONS; i++) {
      HCbutton[i]->color(FL_DARK_BLUE);
      HCbutton[i]->box(FL_UP_BOX);
      HCbutton[i]->labelsize(16);
      HCbutton[i]->labelfont(FL_HELVETICA_BOLD);
      HCbutton[i]->labelcolor(FL_WHITE);
    }

    x_ref = 1300;

    HCexecute = new Fl_Button(x_ref, 30, 100, 30, "execute");
    HCexecute->callback(execute_hydrocal_callback);
    HCexecute->color(FL_GREEN);
    HCexecute->box(FL_PLASTIC_UP_BOX);
    HCexecute->labelsize(18);
    HCexecute->labelfont(FL_HELVETICA_BOLD);
    HCexecute->labelcolor(FL_BLACK);

    HCplot = new Fl_Button(x_ref + 110, 30, 100, 30, "plot");
    HCplot->callback(plot_callback);
    HCplot->color(FL_CYAN);
    HCplot->box(FL_PLASTIC_UP_BOX);
    HCplot->labelsize(18);
    HCplot->labelfont(FL_HELVETICA_BOLD);
    HCplot->labelcolor(FL_BLACK);
    HCplot->hide();

    HCcancel = new Fl_Button(x_ref + 220, 30, 100, 30, "cancel");
    HCcancel->callback(cancel_callback);
    HCcancel->color(FL_RED);
    HCcancel->box(FL_PLASTIC_UP_BOX);
    HCcancel->labelsize(18);
    HCcancel->labelfont(FL_HELVETICA_BOLD);
    HCcancel->labelcolor(FL_BLACK);
    HCcancel->hide();

    HCclose = new Fl_Button(1560, 770, 120, 30, "close");
    HCclose->callback(gui_close_callback);
    HCclose->color(FL_RED);
    HCclose->box(FL_ROUND_UP_BOX);
    HCclose->labelsize(24);
    HCclose->labelfont(FL_HELVETICA_BOLD);
    HCclose->labelcolor(FL_BLACK);

    HCtextdisplay = new Fl_Text_Display(1000, 80, 690, 660);
    HCtextdisplay->color(FL_GRAY);
    HCtextdisplay->textfont(FL_COURIER);
    HCtextdisplay->textsize(11);
    HCtextbuffer = new Fl_Text_Buffer();
    // HCtextbuffer->add_modify_callback(text_buffer_modified_callback,
    // nullptr);
    HCtextdisplay->buffer(HCtextbuffer);

    HCprogress = new Fl_Progress(1000, 745, 690, 20);
    HCprogress->minimum(0);
    HCprogress->maximum(100);
    HCprogress->value(0);
    HCprogress->selection_color(FL_DARK_BLUE);
    HCprogress->labelcolor(FL_BLACK);
    HCprogress->labelsize(10);
    HCprogress->hide();

    HCwin->end();
  }
  if (!HCwin->shown())
    HCwin->show();
  HCwin->redraw();
  HCchoice = none;

  // Run forever and wait for user interaction. The program terminates if the
  // user hits the "quit" button.
  return Fl::run();
} // end main
