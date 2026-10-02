/**
 * @file readstrz.ccx
 *
 * @brief Conversion of files in Strahlenzentrum (STRZ) data format
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: readstrz.cxx 2039 2026-07-20 07:57:32Z iamp $
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */

// #include <stdfloat> //float type of defined size (requires at least C++23)
#include <chrono> // clocks and time
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#ifdef _WIN32
#include <stdlib.h>
#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN 1234
#define BIG_ENDIAN 4321
#endif
#define le16toh(x) (x)
#define be16toh(x) _byteswap_ushort(x)
#define le32toh(x) (x)
#define be32toh(x) _byteswap_ulong(x)
#define le64toh(x) (x)
#define be64toh(x) _byteswap_uint64(x)
#else
#include <endian.h>
#endif
#include "matrix.h"
#include "svnrevision.h"

using namespace std;

namespace STRZ {

////////////////////////////////////////////////////////////////////////////////
/**
 * @struct HEADERstddat
 * @brief File header of abs and scan data files in STRZ format
 *
 * for details see http://www.strz.uni-giessen.de/ExpHelp/esw/esw.html
 */
#define lIDHDR 8
#define lHDLEN 1
#define lEXPMNT 6
#define lIDPRG 8
#define lSTDAT 9
#define lSTTIM 8
#define lSPDAT 9
#define lSPTIM 8
#define lSPENAM 8
#define lSPTYPE 4
#define lROWS 6
#define lCOLS 6
#define lBYTES 1
#define lHDFREE 4
#define lRESRV 38
#define lLTXT 4
#define lTEXT 80

typedef struct {
  char idhdr[lIDHDR];   /* Identification of header: "STRZ-VXW" */
  char hdlen[lHDLEN];   /* Length of header: "1" */
  char expmnt[lEXPMNT]; /* Experiment */
  char idprg[lIDPRG];   /* ID of generating Program: "ESW " */
  char stdat[lSTDAT];   /* Date of start */
  char sttim[lSTTIM];   /* Time of start */
  char spdat[lSPDAT];   /* Date of stop */
  char sptim[lSPTIM];   /* Time of stop */
  char spenam[lSPENAM]; /* Name of spectrum */
  char sptype[lSPTYPE]; /* Type of spectrum: "MCA2" */
  char rows[lROWS];     /* Number of rows: "     5" */
  char cols[lCOLS];     /* Channels/row: " <var>" */
  char bytes[lBYTES];   /* Bytes/channel: "4" */
  char hdfree[lHDFREE]; /* First free byte in header (0,...) */
  char resrv[lRESRV];   /* Reserved */
  char ltxt[lLTXT];     /* Length of text: "80" */
  char text[lTEXT];     /* Text */
} HEADERstddat;

////////////////////////////////////////////////////////////////////////////////
/**
 * @enum StrzDataFormat
 * @brief Data formats of abs and scan data files
 */
enum StrzDataFormat {
  ESI_VAX = 1,  ///< old abs data format on VAX computer
  ESW_VAX = 2,  ///< new abs data format on VAX computer
  ESS_VAX = 3,  ///< scan data format on VAX computer
  ESW_VXW = 4,  ///< new abs data format on VX-works computer
  ESS_VXW = 5,  ///< scan data format on motorola VX-works computer
  MASS_VXW = 6, ///< mass data format on motorola VX-works computer
  ESW_2021 = 7, ///< abs data format 2021 (multi-mode e-gan, VX works)
  ESS_2021 = 8, ///< scan data format 2021 (multi-mode e-gan, VX works)
  MASS_VXI = 9  ///< mass data format on intel VX-works computer
};

uint16_t read_uint16(ifstream &in, int byteorder) {
  uint16_t tmp = 0;
  in.read(reinterpret_cast<char *>(&tmp), 2);
  return byteorder == LITTLE_ENDIAN ? le16toh(tmp) : be16toh(tmp);
}

uint32_t read_uint32(ifstream &in, int byteorder) {
  uint32_t tmp = 0;
  in.read(reinterpret_cast<char *>(&tmp), 4);
  return byteorder == LITTLE_ENDIAN ? le32toh(tmp) : be32toh(tmp);
}

uint64_t read_uint64(ifstream &in, int byteorder) {
  uint64_t tmp = 0;
  in.read(reinterpret_cast<char *>(&tmp), 8);
  return byteorder == LITTLE_ENDIAN ? le64toh(tmp) : be64toh(tmp);
}

float read_float32(ifstream &in, int byteorder) {
  uint32_t tmp = read_uint32(in, byteorder);
  float val = 0.0;
  unsigned int size = 4;
  if (size > sizeof(val))
    size = sizeof(val);
  memcpy(&val, &tmp, size);
  return val; // return type float32_t would be safer, but this would requires
              // C++23
}

double read_float64(ifstream &in, int byteorder) {
  uint64_t tmp = read_uint64(in, byteorder);
  double val = 0.0;
  unsigned int size = 8;
  if (size > sizeof(val))
    size = sizeof(val);
  memcpy(&val, &tmp, size);
  return val; // return type float64_t would be safer, but this would requires
              // C++23
}

float read_vax_float32(ifstream &in) {
  uint32_t fourbytes = 0;
  in.read(reinterpret_cast<char *>(&fourbytes), 4);
  uint16_t twobytes0 = static_cast<uint16_t>(fourbytes & 0x0000FFFF);
  uint16_t sign = twobytes0 & 0x80;
  uint16_t twobytes1 = static_cast<uint16_t>((fourbytes & 0xFFFF0000) >> 16);
  twobytes0 = (twobytes0 & 0x7FFF);
  if (twobytes0 < 0x100)
    return 0.0; // underflow
  twobytes0 -= 0x100;
  twobytes0 |= sign;
  uint8_t byte0 = static_cast<uint8_t>(twobytes0 & 0x00FF);
  uint8_t byte1 = static_cast<uint8_t>((twobytes0 & 0xFF00) >> 8);
  uint8_t byte2 = static_cast<uint8_t>(twobytes1 & 0x00FF);
  uint8_t byte3 = static_cast<uint8_t>((twobytes1 & 0xFF00) >> 8);
  fourbytes =
      byte2 +
      0x100 * (byte3 + 0x100 * (byte0 + 0x100 * byte1)); // little  endian
  float f;
  memcpy(&f, &fourbytes, sizeof(f));
  return f;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Returns the electron current range
 *
 * @param rangeID ID of the range to be used
 * @return The full-scale value of the selected range in A
 */
double EcurRange(double rangeID) { // returns electron current range in A
  double range = -1.0;
  int selector = round(rangeID);
  /** @code */
  switch (selector) {
  case -1:
    range = 1.5;
    break; // 1500 mA
  case 0:
    range = 0.5;
    break; // 500 mA
  case 1:
    range = 0.15;
    break; // 150 mA
  case 2:
    range = 0.05;
    break; // 50 mA
  case 3:
    range = 0.015;
    break; // 15 mA
  case 4:
    range = 5e-3;
    break; // 5 mA
  case 5:
    range = 1.5e-3;
    break; // 1.5 mA
  case 6:
    range = 500e-6;
    break; // 500 muA
  case 7:
    range = 150e-6;
    break; // 150 muA
  default:
    range = -1.0;
    break; // indicates wrong selector
  }
  return range;
  /** @endcode */
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Returns the ion current range
 *
 * @param rangeID ID of the range to be used
 * @return The full-scale value of the selected range in A
 */
double IcurRange(double rangeID) { // returns ion current range in Amperes
  double range = -1.0;
  int selector = round(rangeID);
  /** @code */
  switch (selector) {
  case 0:
    range = 100e-6;
    break; // 100 muA
  case 1:
    range = 30e-6;
    break; // 30 muA
  case 2:
    range = 10e-6;
    break; // 10 muA
  case 3:
    range = 3e-6;
    break; // 3 muA
  case 4:
    range = 1e-6;
    break; // 1 muA
  case 5:
    range = 300e-9;
    break; // 300 nA
  case 6:
    range = 100e-9;
    break; // 100 nA
  case 7:
    range = 30e-9;
    break; // 30 nA
  case 8:
    range = 10e-9;
    break; // 10 nA
  case 9:
    range = 3e-9;
    break; // 3 nA
  case 10:
    range = 1e-9;
    break; // 1 nA
  case 11:
    range = 300e-12;
    break; // 300 pA
  case 12:
    range = 100e-12;
    break; // 100 pA
  case 13:
    range = 30e-12;
    break; // 30 pA
  case 14:
    range = 10e-12;
    break; // 10 pA
  case 15:
    range = 3e-12;
    break; // 3 pA
  default:
    range = -1.0;
    break; // indicates wrong selector
  }
  return range;
  /** @endcode */
}

/////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Reads binary data files from Giessen electron-ion crossed-beams setup
 * and converts to ascii
 *
 * @return 0 file processed ok
 * @return 1 file does not exist
 * @return 2 unknown data format detected
 * @return 3 end of reading files requested
 */
int ConvertStrzDataFile(void) {
  string filename = "";
  cout << endl << " Give the filename of the STRZ data file ('q' quits) ...: ";
  cin >> filename;
  if (filename == "q")
    return 3;
  ifstream strzfile(filename, ios::binary);
  if (!strzfile) {
    cout << endl << " ERROR:: Could not open " << filename << endl;
    return 1;
  }
  string ascfilename = filename + ".asc";

  // time stamp
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  struct Parameter {
    string name;
    string value;
  };
  vector<Parameter> params;
  params.push_back(
      {"hydrocal SVN revision", SVNrevision}); // defined in svnrevision.h
  string outstr = asctime(timeinfo);
  outstr.resize(outstr.size() - 1); // remove newline character
  params.push_back({"date & time of conversion", outstr});
  params.push_back({"STRZ filename", filename});
  params.push_back({"ascii filename", ascfilename});

  // first the standard header is read
  HEADERstddat headerstd;
  strzfile.read(reinterpret_cast<char *>(&headerstd), sizeof(headerstd));

  outstr = headerstd.text;
  outstr.resize(lTEXT);
  params.push_back({"comment", outstr});

  outstr = headerstd.idhdr;
  outstr.resize(lIDHDR);
  params.push_back({"identification of header", outstr});
  bool vax_data_flag =
      (outstr.find("VAX") < string::npos); // indicates VAX data format
  bool vxi_data_flag =
      (outstr.find("VXI") < string::npos); // indicates VXI data format
  bool vxw_data_flag =
      (outstr.find("VXW") <
       string::npos); // indicates VXW data format  params.push_back({"hydrocal
                      // date & time", asctime (timeinfo)});

  int byteorder = (vax_data_flag || vxi_data_flag) ? LITTLE_ENDIAN : BIG_ENDIAN;

  outstr = headerstd.hdlen;
  outstr.resize(lHDLEN);
  params.push_back({"length of header", outstr});

  outstr = headerstd.expmnt;
  outstr.resize(lEXPMNT);
  params.push_back({"name of experiment", outstr});

  outstr = headerstd.idprg;
  outstr.resize(lIDPRG);
  params.push_back({"ID of generating program", outstr});

  unsigned int data_format = 0;
  bool abs_data_flag = false;
  bool scan_data_flag = false;
  bool mass_data_flag = false;
  if (outstr.find("ESI") < string::npos) {
    abs_data_flag = true;
    data_format = ESI_VAX;
    params.push_back({"data format", "ESI-VAX"});
  } else if (outstr.find("ESW") < string::npos) {
    abs_data_flag = true;
    if (vax_data_flag) {
      data_format = ESW_VAX;
      params.push_back({"data format", "ESW-VAX"});
    } else {
      if (outstr.find("2021") < string::npos) {
        data_format = ESW_2021;
        params.push_back({"data format", "ESW-2021"});
      } else {
        data_format = ESW_VXW;
        params.push_back({"data format", "ESW-VXW"});
      }
    }
  } else if (outstr.find("ESS") < string::npos) {
    scan_data_flag = true;
    if (vax_data_flag) {
      data_format = ESS_VAX;
      params.push_back({"data format", "ESS-VAX"});
    } else {
      if (outstr.find("2021") < string::npos) {
        data_format = ESS_2021;
        params.push_back({"data format", "ESS-2021"});
      } else {
        data_format = ESS_VXW;
        params.push_back({"data format", "ESS-VXW"});
      }
    }
  } else if ((outstr.find("MASS") < string::npos) ||
             (outstr.find("MAG2") < string::npos)) {
    mass_data_flag = true;
    if (vax_data_flag) {
      cout << endl << " ERROR:: Could not identify data format!" << endl;
      return 2; // unknown data format
    } else if (vxw_data_flag) {
      data_format = MASS_VXW;
      params.push_back({"data format", "MASS-VXW"});
    } else if (vxi_data_flag) {
      data_format = MASS_VXI;
      params.push_back({"data format", "MASS-VXI"});
    }
  } else {
    cout << endl << " ERROR:: Could not identify data format!" << endl;
    for (const Parameter &p : params) {
      cout << " " << setw(50) << p.name << ": " << p.value << endl;
    }
    return 2; // unknown data format
  }

  if (abs_data_flag) {
    params.push_back({"type of measurement", "ABS"});
  } else if (scan_data_flag) {
    params.push_back({"type of measurement", "SCAN"});
  } else if (mass_data_flag) {
    params.push_back({"type of measurement", "MASS"});
  }

  outstr = headerstd.stdat;
  outstr.resize(lSTDAT);
  params.push_back({"start date", outstr});

  outstr = headerstd.sttim;
  outstr.resize(lSTTIM);
  params.push_back({"start time", outstr});

  outstr = headerstd.spdat;
  outstr.resize(lSPDAT);
  params.push_back({"stop date", outstr});

  outstr = headerstd.sptim;
  outstr.resize(lSPTIM);
  params.push_back({"stop time", outstr});

  outstr = headerstd.spenam;
  outstr.resize(lSPENAM);
  params.push_back({"name of spectrum", outstr});

  outstr = headerstd.sptype;
  outstr.resize(lSPTYPE);
  params.push_back({"spec type", outstr});

  outstr = headerstd.rows;
  outstr.resize(lROWS);
  unsigned int nrows = stoi(outstr);
  params.push_back({"number of rows", to_string(nrows)});

  outstr = headerstd.cols;
  outstr.resize(lCOLS);
  unsigned int ncols = stoi(outstr);
  params.push_back({"number of channels per row", to_string(ncols)});

  outstr = headerstd.bytes;
  outstr.resize(lBYTES);
  unsigned int nbytes = stoi(outstr);
  params.push_back({"number of bytes per channel", to_string(nbytes)});

  // ---------------------- end of standard header ------------------

  // ---------------------- start of special header -----------------

  // global variables
  unsigned int gSHDguntyp = 0;
  float gSHDgunpar[10];
  for (unsigned int i = 0; i < 10; i++)
    gSHDgunpar[i] = 0.0;
  double icur_range = 0, ecur_range = 0, time_base = 0;

  ///////////////////////////////////////////////////// ESI_VAX
  /////////////////////////////////
  if (data_format == ESI_VAX) {
    read_uint16(strzfile, byteorder); // 2 bytes
    uint32_t SHDspelen = read_uint32(strzfile, byteorder); // 4 bytes
    uint16_t SHDblocks = read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDblockm = read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDblockd = read_uint16(strzfile, byteorder); // 2 bytes
    uint32_t SHDrltcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDlftcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint16_t SHDdatcnt = read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDoutcnt = read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDct1cnt = read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDct2cnt = read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDct3cnt = read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDseqcnt = read_uint16(strzfile, byteorder); // 2 bytes
    uint32_t SHDrejcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDerrcnt = read_uint32(strzfile, byteorder); // 4 bytes
    read_uint32(strzfile, byteorder); // 4 bytes
    read_uint16(strzfile, byteorder); // 2 bytes
    float EXPelo = read_vax_float32(strzfile);             //  0
    float EXPecurr = read_vax_float32(strzfile);           //  1
    float EXPionQ = read_vax_float32(strzfile);            //  2
    float EXPionA = read_vax_float32(strzfile);            //  3
    float EXPionE = read_vax_float32(strzfile);            //  4
    float EXPdivSIG = read_vax_float32(strzfile);          //  5
    float EXPdivIcur = read_vax_float32(strzfile);         //  6
    float EXPconv = read_vax_float32(strzfile);            //  7
    float EXPconv1 = read_vax_float32(strzfile);           //  8
    float EXPconv10 = read_vax_float32(strzfile);          //  9
    float EXPconv100 = read_vax_float32(strzfile);         // 10
    float EXPeff = read_vax_float32(strzfile);             // 11
    float EXPdecmm = read_vax_float32(strzfile);           // 12
    float EXPdecdiv = read_vax_float32(strzfile);          // 13

    params.push_back({"electron energy (eV)", to_string(EXPelo)});
    params.push_back({"electron current (A)", to_string(EXPecurr)});
    params.push_back({"ion charge", to_string(EXPionQ)});
    params.push_back({"ion mass (u)", to_string(EXPionA)});
    params.push_back({"ion energy (keV)", to_string(EXPionE)});
    params.push_back({"detection efficiency (%)", to_string(EXPeff)});
    params.push_back({"divider of signal", to_string(EXPdivSIG)});
    params.push_back({"divider of ion current", to_string(EXPdivIcur)});
    params.push_back(
        {"ion current conversion (0.1 nA/KHz)", to_string(EXPconv)});
    params.push_back(
        {"ion current conversion (1.0 nA/KHz)", to_string(EXPconv1)});
    params.push_back(
        {"ion current conversion (10.0 nA/KHz)", to_string(EXPconv10)});
    params.push_back(
        {"ion current conversion (100.0 nA/KHz)", to_string(EXPconv100)});
    params.push_back({"decoder counts divided by", to_string(EXPdecdiv)});
    params.push_back({"position step size (mm)", to_string(EXPdecmm)});

    params.push_back({"length of single spectrum", to_string(SHDspelen)});
    params.push_back({"blocks s", to_string(SHDblocks)});
    params.push_back({"blocks m", to_string(SHDblockm)});
    params.push_back({"blocks d", to_string(SHDblockd)});
    params.push_back({"real time", to_string(SHDrltcnt)});
    params.push_back({"lifetime", to_string(SHDlftcnt)});
    params.push_back({"processed positions", to_string(SHDdatcnt)});
    params.push_back({"positions out of range", to_string(SHDoutcnt)});
    params.push_back({"processed counter 1 data", to_string(SHDct1cnt)});
    params.push_back({"processed counter 2 data", to_string(SHDct2cnt)});
    params.push_back({"processed counter 3 data", to_string(SHDct3cnt)});
    params.push_back({"rejected data", to_string(SHDrejcnt)});
    params.push_back({"sequence errors", to_string(SHDseqcnt)});
    params.push_back({"error counter", to_string(SHDerrcnt)});

  } /////////////////////////////////////////////////// ESW_VAX
    /////////////////////////////////
  else if (data_format == ESW_VAX) {
    read_uint16(strzfile, byteorder); // 2 bytes
    read_uint32(strzfile, byteorder); // 4 bytes
    uint16_t SHDblocks = read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDblockm = read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDblockd = read_uint16(strzfile, byteorder); // 2 bytes
    uint32_t SHDrltcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDlftcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDdatcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDoutcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct1cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct2cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct3cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct4cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDseqcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDfulcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDrejcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDerrcnt = read_uint32(strzfile, byteorder); // 4 bytes
    read_uint32(strzfile, byteorder); // 4 bytes
    read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDslen = read_uint16(strzfile, byteorder);   // 2 bytes
    float EXPelo = read_vax_float32(strzfile);             //  0
    read_vax_float32(strzfile);             //  1
    float EXPionQ = read_vax_float32(strzfile);            //  2
    float EXPionA = read_vax_float32(strzfile);            //  3
    float EXPionE = read_vax_float32(strzfile);            //  4
    float EXPeff = read_vax_float32(strzfile);             //  5
    float EXPtimbas = read_vax_float32(strzfile);          //  6
    float EXPeconvR = read_vax_float32(strzfile);          //  7
    float EXPeconvFS = read_vax_float32(strzfile);         //  8
    float EXPiconvR = read_vax_float32(strzfile);          //  9
    float EXPiconvFS = read_vax_float32(strzfile);         // 10
    float EXPdecdiv = read_vax_float32(strzfile);          // 11
    float EXPdecmm = read_vax_float32(strzfile);           // 12
    char SHDecfprog[13];
    strzfile.read(SHDecfprog, 12);
    SHDecfprog[12] = '\0';
    gSHDguntyp = read_uint16(strzfile, byteorder);
    float SHDulinse = read_vax_float32(strzfile);
    char SHDnotusd[31];
    strzfile.read(SHDnotusd, 30);
    SHDnotusd[30] = '\0';
    float SHDs5scal = read_vax_float32(strzfile);

    // params.push_back({"version of ECF program", SHDecfprog});
    params.push_back({"electron energy (eV)", to_string(EXPelo)});
    params.push_back({"ion charge", to_string(EXPionQ)});
    params.push_back({"ion mass (u)", to_string(EXPionA)});
    params.push_back({"ion energy (keV)", to_string(EXPionE)});
    params.push_back({"detection efficiency (%)", to_string(EXPeff)});
    params.push_back({"ele current range ID", to_string(EXPeconvR)});
    params.push_back(
        {"ele current conversion factor (Hz /f.s.)", to_string(EXPeconvFS)});
    params.push_back({"ion current range ID", to_string(EXPiconvR)});
    params.push_back(
        {"ion current conversion factor (Hz /f.s.)", to_string(EXPiconvFS)});
    params.push_back({"decoder counts divided by", to_string(EXPdecdiv)});
    params.push_back({"position step size (mm)", to_string(EXPdecmm)});
    params.push_back({"lens voltage (V)", to_string(SHDulinse)});
    params.push_back({"scaling factor spectrum 5", to_string(SHDs5scal)});

    params.push_back({"length of single spectrum", to_string(SHDslen)});
    params.push_back({"blocks s", to_string(SHDblocks)});
    params.push_back({"blocks m", to_string(SHDblockm)});
    params.push_back({"blocks d", to_string(SHDblockd)});
    params.push_back({"real time", to_string(SHDrltcnt)});
    params.push_back({"lifetime", to_string(SHDlftcnt)});
    params.push_back({"processed positions", to_string(SHDdatcnt)});
    params.push_back({"positions out of range", to_string(SHDoutcnt)});
    params.push_back({"processed counter 1 data", to_string(SHDct1cnt)});
    params.push_back({"processed counter 2 data", to_string(SHDct2cnt)});
    params.push_back({"processed counter 3 data", to_string(SHDct3cnt)});
    params.push_back({"processed counter 4 data", to_string(SHDct4cnt)});
    params.push_back({"fifo full counter", to_string(SHDfulcnt)});
    params.push_back({"rejected data", to_string(SHDrejcnt)});
    params.push_back({"sequence errors", to_string(SHDseqcnt)});
    params.push_back({"error counter", to_string(SHDerrcnt)});

    ecur_range = EcurRange(EXPeconvR);
    icur_range = IcurRange(EXPiconvR);
    time_base = 1e6 * pow(2, -round(EXPtimbas)); // in Hz

  } /////////////////////////////////////////////////// ESS_VAX
    /////////////////////////////////
  else if (data_format == ESS_VAX) {
    read_uint16(strzfile, byteorder); // 2 bytes
    read_uint32(strzfile, byteorder); // 4 bytes
    uint16_t SHDblocks = read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDblockm = read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDblockd = read_uint16(strzfile, byteorder); // 2 bytes
    uint32_t SHDrltcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDlftcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDdatcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDoutcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct1cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct2cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct3cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct4cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDseqcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDfulcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDrejcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDerrcnt = read_uint32(strzfile, byteorder); // 4 bytes
    read_uint32(strzfile, byteorder); // 4 bytes
    read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDslen = read_uint16(strzfile, byteorder);   // 2 bytes
    float EXPelo = read_vax_float32(strzfile);             //  0
    float EXPehi = read_vax_float32(strzfile);             //  1
    float EXPionQ = read_vax_float32(strzfile);            //  2
    float EXPionA = read_vax_float32(strzfile);            //  3
    float EXPionE = read_vax_float32(strzfile);            //  4
    float EXPeff = read_vax_float32(strzfile);             //  5
    float EXPtimbas = read_vax_float32(strzfile);          //  6
    float EXPeconvR = read_vax_float32(strzfile);          //  7
    float EXPeconvFS = read_vax_float32(strzfile);         //  8
    float EXPiconvR = read_vax_float32(strzfile);          //  9
    float EXPiconvFS = read_vax_float32(strzfile);         // 10
    float EXPtstart = read_vax_float32(strzfile);          // 11
    float EXPtwait = read_vax_float32(strzfile);           // 12
    float EXPtmeas = read_vax_float32(strzfile);           // 13
    char SHDecfprog[13];
    strzfile.read(SHDecfprog, 12);
    SHDecfprog[12] = '\0';
    gSHDguntyp = read_uint16(strzfile, byteorder);
    float SHDulinse = read_vax_float32(strzfile);
    char SHDnotusd[31];
    strzfile.read(SHDnotusd, 30);
    SHDnotusd[30] = '\0';
    float SHDs5scal = read_vax_float32(strzfile);

    params.push_back({"version of ECF program", SHDecfprog});
    params.push_back({"waiting time at startup (s)", to_string(EXPtstart)});
    params.push_back({"min energy (eV)", to_string(EXPelo)});
    params.push_back({"max energy (eV)", to_string(EXPehi)});
    params.push_back(
        {"waiting time before measurement gate (ms)", to_string(EXPtwait)});
    params.push_back(
        {"duration of measurement gate (ms)", to_string(EXPtmeas)});
    params.push_back({"ion charge", to_string(EXPionQ)});
    params.push_back({"ion charge", to_string(EXPionQ)});
    params.push_back({"ion mass (u)", to_string(EXPionA)});
    params.push_back({"ion energy (keV)", to_string(EXPionE)});
    params.push_back({"detection efficiency (%)", to_string(EXPeff)});
    params.push_back(
        {"time base exponent n (2**(-n) MHz)", to_string(EXPtimbas)});
    params.push_back({"ele current range ID", to_string(EXPeconvR)});
    params.push_back(
        {"ele current conversion factor (Hz /f.s.)", to_string(EXPeconvFS)});
    params.push_back({"ion current range ID", to_string(EXPiconvR)});
    params.push_back(
        {"ion current conversion factor (Hz /f.s.)", to_string(EXPiconvFS)});
    params.push_back({"lens voltage (V)", to_string(SHDulinse)});
    params.push_back({"scaling factor spectrum 5", to_string(SHDs5scal)});

    params.push_back({"length of single spectrum", to_string(SHDslen)});
    params.push_back({"blocks s", to_string(SHDblocks)});
    params.push_back({"blocks m", to_string(SHDblockm)});
    params.push_back({"blocks d", to_string(SHDblockd)});
    params.push_back({"real time", to_string(SHDrltcnt)});
    params.push_back({"lifetime", to_string(SHDlftcnt)});
    params.push_back({"processed positions", to_string(SHDdatcnt)});
    params.push_back({"positions out of range", to_string(SHDoutcnt)});
    params.push_back({"processed counter 1 data", to_string(SHDct1cnt)});
    params.push_back({"processed counter 2 data", to_string(SHDct2cnt)});
    params.push_back({"processed counter 3 data", to_string(SHDct3cnt)});
    params.push_back({"processed counter 4 data", to_string(SHDct4cnt)});
    params.push_back({"fifo full counter", to_string(SHDfulcnt)});
    params.push_back({"rejected data", to_string(SHDrejcnt)});
    params.push_back({"sequence errors", to_string(SHDseqcnt)});
    params.push_back({"error counter", to_string(SHDerrcnt)});

    ecur_range = EcurRange(EXPeconvR);
    icur_range = IcurRange(EXPiconvR);
    time_base = 1e6 * pow(2, -round(EXPtimbas)); // in Hz

  } /////////////////////////////////////////////////// ESW_VXW
    /////////////////////////////////
  else if (data_format == ESW_VXW) {
    read_uint16(strzfile, byteorder); // 2 bytes
    uint32_t SHDrltcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDlftcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDdatcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDoutcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct1cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct2cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct3cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct4cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDseqcnt = read_uint32(strzfile, byteorder); // 4 bytes
    read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDfulcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDrejcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDerrcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint16_t SHDslen = read_uint16(strzfile, byteorder);   // 2 bytes
    float EXPelo = read_float32(strzfile, byteorder);      //  0
    read_float32(strzfile, byteorder);      //  1
    float EXPionQ = read_float32(strzfile, byteorder);     //  2
    float EXPionA = read_float32(strzfile, byteorder);     //  3
    float EXPionE = read_float32(strzfile, byteorder);     //  4
    float EXPeff = read_float32(strzfile, byteorder);      //  5
    float EXPtimbas = read_float32(strzfile, byteorder);   //  6
    float EXPeconvR = read_float32(strzfile, byteorder);   //  7
    float EXPeconvFS = read_float32(strzfile, byteorder);  //  8
    float EXPiconvR = read_float32(strzfile, byteorder);   //  9
    float EXPiconvFS = read_float32(strzfile, byteorder);  // 10
    float EXPdecdiv = read_float32(strzfile, byteorder);   // 11
    float EXPdecmm = read_float32(strzfile, byteorder);    // 12
    char SHDecfprog[13];
    strzfile.read(SHDecfprog, 12);
    SHDecfprog[12] = '\0';
    gSHDguntyp = read_uint16(strzfile, byteorder);
    for (int i = 0; i < 10; i++) {
      gSHDgunpar[i] = read_float32(strzfile, byteorder);
    }
    float SHDdeadtm = read_float32(strzfile, byteorder);
    float SHDdtmerr = read_float32(strzfile, byteorder);
    float SHDerrkf = read_float32(strzfile, byteorder);
    float SHDerrec = read_float32(strzfile, byteorder);
    float SHDerric = read_float32(strzfile, byteorder);
    float SHDerrcw = read_float32(strzfile, byteorder);
    float SHDerrde = read_float32(strzfile, byteorder);
    float SHDs5scal = read_float32(strzfile, byteorder);
    uint32_t SHDruntim = read_uint32(strzfile, byteorder);

    params.push_back({"electron energy (eV)", to_string(EXPelo)});
    params.push_back({"ion charge", to_string(EXPionQ)});
    params.push_back({"ion charge", to_string(EXPionQ)});
    params.push_back({"ion mass (u)", to_string(EXPionA)});
    params.push_back({"ion energy (keV)", to_string(EXPionE)});
    params.push_back({"detection efficiency (%)", to_string(EXPeff)});
    params.push_back(
        {"time base exponent n (2**(-n) MHz)", to_string(EXPtimbas)});
    params.push_back({"ele current range ID", to_string(EXPeconvR)});
    params.push_back(
        {"ele current conversion factor (Hz /f.s.)", to_string(EXPeconvFS)});
    params.push_back({"ion current range ID", to_string(EXPiconvR)});
    params.push_back(
        {"ion current conversion factor (Hz /f.s.)", to_string(EXPiconvFS)});
    params.push_back(
        {"counts from angular encoder divided by", to_string(EXPdecdiv)});
    params.push_back({"pulses per 10 mm", to_string(EXPdecmm)});
    params.push_back({"scaling factor spectrum 5", to_string(SHDs5scal)});

    params.push_back({"dead time (us)", to_string(SHDdeadtm)});
    params.push_back({"error dead time (us)", to_string(SHDdtmerr)});
    params.push_back({"dead time (us)", to_string(SHDdeadtm)});
    params.push_back({"error dead time (us)", to_string(SHDdtmerr)});
    params.push_back({"error kinectic factor (%)", to_string(SHDerrkf)});
    params.push_back({"error electron current (%)", to_string(SHDerrec)});
    params.push_back({"error ion current (%)", to_string(SHDerric)});
    params.push_back({"error channel width (%)", to_string(SHDerrcw)});
    params.push_back({"error detection efficiency (%)", to_string(SHDerrde)});

    params.push_back({"length of single spectrum", to_string(SHDslen)});
    params.push_back({"run time", to_string(SHDruntim)});
    params.push_back({"real time", to_string(SHDrltcnt)});
    params.push_back({"lifetime", to_string(SHDlftcnt)});
    params.push_back({"processed positions", to_string(SHDdatcnt)});
    params.push_back({"positions out of range", to_string(SHDoutcnt)});
    params.push_back({"processed counter 1 data", to_string(SHDct1cnt)});
    params.push_back({"processed counter 2 data", to_string(SHDct2cnt)});
    params.push_back({"processed counter 3 data", to_string(SHDct3cnt)});
    params.push_back({"processed counter 4 data", to_string(SHDct4cnt)});
    params.push_back({"fifo full counter", to_string(SHDfulcnt)});
    params.push_back({"rejected data", to_string(SHDrejcnt)});
    params.push_back({"sequence errors", to_string(SHDseqcnt)});
    params.push_back({"error counter", to_string(SHDerrcnt)});

    ecur_range = EcurRange(EXPeconvR);
    icur_range = IcurRange(EXPiconvR);
    time_base = 1e6 * pow(2, -round(EXPtimbas)); // in Hz

  } /////////////////////////////////////////////////// ESS_VXW
    /////////////////////////////////
  else if (data_format == ESS_VXW) {
    read_uint16(strzfile, byteorder); // 2 bytes
    uint32_t SHDrltcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDlftcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDdatcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDoutcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct1cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct2cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct3cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct4cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDseqcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDbovcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDfulcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDrejcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDerrcnt = read_uint32(strzfile, byteorder); // 4 bytes
    read_uint32(strzfile, byteorder); // 4 bytes
    read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDslen = read_uint16(strzfile, byteorder);   // 2 bytes
    float EXPelo = read_float32(strzfile, byteorder);      //  0
    float EXPehi = read_float32(strzfile, byteorder);      //  1
    float EXPionQ = read_float32(strzfile, byteorder);     //  2
    float EXPionA = read_float32(strzfile, byteorder);     //  3
    float EXPionE = read_float32(strzfile, byteorder);     //  4
    float EXPeff = read_float32(strzfile, byteorder);      //  5
    float EXPtimbas = read_float32(strzfile, byteorder);   //  6
    float EXPeconvR = read_float32(strzfile, byteorder);   //  7
    float EXPeconvFS = read_float32(strzfile, byteorder);  //  8
    float EXPiconvR = read_float32(strzfile, byteorder);   //  9
    float EXPiconvFS = read_float32(strzfile, byteorder);  // 10
    float EXPtstart = read_float32(strzfile, byteorder);   // 11
    float EXPtwait = read_float32(strzfile, byteorder);    // 12
    float EXPtmeas = read_float32(strzfile, byteorder);    // 13
    char SHDecfprog[13]{};
    strzfile.read(SHDecfprog, 12);
    SHDecfprog[12] = '\0';
    gSHDguntyp = read_uint16(strzfile, byteorder);
    char SHDguntxt[33];
    strzfile.read(SHDguntxt, 32);
    SHDguntxt[32] = '\0';
    for (int i = 0; i < 10; i++) {
      gSHDgunpar[i] = read_float32(strzfile, byteorder);
    }
    float SHDdeadtm = read_float32(strzfile, byteorder);
    float SHDdtmerr = read_float32(strzfile, byteorder);
    uint32_t SHDruntim = read_uint32(strzfile, byteorder);
    params.push_back({"version of ECF program", SHDecfprog});
    params.push_back({"waiting time at startup (s)", to_string(EXPtstart)});
    if (strncmp(SHDecfprog, "ECFscd.x", 7) == 0) {
      // read extended parameters from expar2 array
      uint32_t EXPnchan1 = read_uint32(strzfile, byteorder);
      uint32_t EXPnchan2 = read_uint32(strzfile, byteorder);
      read_float32(strzfile, byteorder);
      float EXPehi2 = read_float32(strzfile, byteorder);
      read_float32(strzfile, byteorder);
      float EXPelo2 = read_float32(strzfile, byteorder);
      read_float32(strzfile, byteorder);
      float EXPtwait2 = read_float32(strzfile, byteorder);
      read_float32(strzfile, byteorder);
      float EXPtmeas2 = read_float32(strzfile, byteorder);

      params.push_back({"electron gun", SHDguntxt});
      params.push_back({"1st number of channels", to_string(EXPnchan1)});
      params.push_back({"1st min energy (eV)", to_string(EXPelo)});
      params.push_back({"1st max energy (eV)", to_string(EXPehi)});
      params.push_back({"1st waiting time before measurement gate (ms)",
                        to_string(EXPtwait)});
      params.push_back(
          {"1st duration of measurement gate (ms)", to_string(EXPtmeas)});
      params.push_back({"2nd number of channels", to_string(EXPnchan2)});
      params.push_back({"2nd min energy 2 (eV)", to_string(EXPelo2)});
      params.push_back({"2nd max energy 2 (eV)", to_string(EXPehi2)});
      params.push_back({"2nd waiting time before measurement gate (ms)",
                        to_string(EXPtwait2)});
      params.push_back(
          {"2nd duration of measurement gate (ms)", to_string(EXPtmeas2)});
    } else {
      params.push_back({"min energy (eV)", to_string(EXPelo)});
      params.push_back({"max energy (eV)", to_string(EXPehi)});
      params.push_back(
          {"waiting time before measurement gate (ms)", to_string(EXPtwait)});
      params.push_back(
          {"duration of measurement gate (ms)", to_string(EXPtmeas)});
    }
    params.push_back({"ion charge", to_string(EXPionQ)});
    params.push_back({"ion charge", to_string(EXPionQ)});
    params.push_back({"ion mass (u)", to_string(EXPionA)});
    params.push_back({"ion energy (keV)", to_string(EXPionE)});
    params.push_back({"detection efficiency (%)", to_string(EXPeff)});
    params.push_back({"time base (Hz)", to_string(EXPtimbas)});
    params.push_back(
        {"time base exponent n (2**(-n) MHz)", to_string(EXPtimbas)});
    params.push_back({"ele current range ID", to_string(EXPeconvR)});
    params.push_back(
        {"ele current conversion factor (Hz /f.s.)", to_string(EXPeconvFS)});
    params.push_back({"ion current range ID", to_string(EXPiconvR)});
    params.push_back(
        {"ion current conversion factor (Hz /f.s.)", to_string(EXPiconvFS)});
    params.push_back({"dead time (us)", to_string(SHDdeadtm)});
    params.push_back({"error dead time (us)", to_string(SHDdtmerr)});

    params.push_back({"length of single spectrum", to_string(SHDslen)});
    params.push_back({"run time", to_string(SHDruntim)});
    params.push_back({"real time", to_string(SHDrltcnt)});
    params.push_back({"lifetime", to_string(SHDlftcnt)});
    params.push_back({"processed positions", to_string(SHDdatcnt)});
    params.push_back({"positions out of range", to_string(SHDoutcnt)});
    params.push_back({"processed counter 1 data", to_string(SHDct1cnt)});
    params.push_back({"processed counter 2 data", to_string(SHDct2cnt)});
    params.push_back({"processed counter 3 data", to_string(SHDct3cnt)});
    params.push_back({"processed counter 4 data", to_string(SHDct4cnt)});
    params.push_back({"fifo full counter", to_string(SHDfulcnt)});
    params.push_back({"buffer overruns", to_string(SHDbovcnt)});
    params.push_back({"rejected data", to_string(SHDrejcnt)});
    params.push_back({"sequence errors", to_string(SHDseqcnt)});
    params.push_back({"error counter", to_string(SHDerrcnt)});

    ecur_range = EcurRange(EXPeconvR);
    icur_range = IcurRange(EXPiconvR);
    time_base = 1e6 * pow(2, -round(EXPtimbas)); // in Hz
  }

  /////////////////////////////////////////////////// ESW_2021
  ////////////////////////////////
  else if (data_format == ESW_2021) {
    read_uint16(strzfile, byteorder); // 2 bytes
    read_uint16(strzfile, byteorder); // 2 bytes for aligment
    uint32_t SHDrltcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDlftcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDdatcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDoutcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct1cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct2cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct3cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct4cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDseqcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDfulcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDrejcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDerrcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint16_t SHDslen = read_uint16(strzfile, byteorder);   // 2 bytes
    read_uint16(strzfile, byteorder);  // 2 bytes for aligment
    read_uint32(strzfile, byteorder);  // 4 bytes for aligment
    double EXPenergy = read_float64(strzfile, byteorder); //  0+1
    float EXPionQ = read_float32(strzfile, byteorder);   //  2
    float EXPionA = read_float32(strzfile, byteorder);   //  3
    float EXPionE = read_float32(strzfile, byteorder);   //  4
    float EXPeff = read_float32(strzfile, byteorder);    //  5
    float EXPtimbas = read_float32(strzfile, byteorder); //  6
    float EXPeconvR = read_float32(strzfile, byteorder); //  7
    float EXPeconvFS = read_float32(strzfile, byteorder); //  8
    float EXPiconvR = read_float32(strzfile, byteorder);  //  9
    float EXPiconvFS = read_float32(strzfile, byteorder); // 10
    float EXPdecdiv = read_float32(strzfile, byteorder);  // 11
    float EXPdecmm = read_float32(strzfile, byteorder);   // 12
    float EXPwwpot = read_float32(strzfile, byteorder);   // 13
    read_float32(strzfile, byteorder);   // 14
    read_uint32(strzfile, byteorder); // 4 bytes for aligment
    gSHDguntyp = read_uint16(strzfile, byteorder);
    read_uint16(strzfile, byteorder); // 2 bytes for aligment
    for (int i = 0; i < 10; i++) {
      gSHDgunpar[i] = read_float32(strzfile, byteorder);
    }
    float SHDdeadtm = read_float32(strzfile, byteorder);
    float SHDdtmerr = read_float32(strzfile, byteorder);
    float SHDerrkf = read_float32(strzfile, byteorder);
    float SHDerrec = read_float32(strzfile, byteorder);
    float SHDerric = read_float32(strzfile, byteorder);
    float SHDerrcw = read_float32(strzfile, byteorder);
    float SHDerrde = read_float32(strzfile, byteorder);
    float SHDs5scal = read_float32(strzfile, byteorder);
    uint32_t SHDruntim = read_uint32(strzfile, byteorder);

    params.push_back({"electron energy (eV)", to_string(EXPenergy)});
    params.push_back({"ion charge", to_string(EXPionQ)});
    params.push_back({"ion charge", to_string(EXPionQ)});
    params.push_back({"ion mass (u)", to_string(EXPionA)});
    params.push_back({"ion energy (keV)", to_string(EXPionE)});
    params.push_back({"detection efficiency (%)", to_string(EXPeff)});
    params.push_back({"time base (Hz)", to_string(EXPtimbas)});
    params.push_back(
        {"time base exponent n (2**(-n) MHz)", to_string(EXPtimbas)});
    params.push_back({"ele current range ID", to_string(EXPeconvR)});
    params.push_back(
        {"ele current conversion factor (Hz /f.s.)", to_string(EXPeconvFS)});
    params.push_back({"ion current range ID", to_string(EXPiconvR)});
    params.push_back(
        {"ion current conversion factor (Hz /f.s.)", to_string(EXPiconvFS)});
    params.push_back(
        {"counts from angular encoder divided by", to_string(EXPdecdiv)});
    params.push_back({"pulses per 10 mm", to_string(EXPdecmm)});
    params.push_back(
        {"potential of interaction region (V)", to_string(EXPwwpot)});
    params.push_back({"scaling factor spectrum 5", to_string(SHDs5scal)});

    params.push_back({"dead time (us)", to_string(SHDdeadtm)});
    params.push_back({"error dead time (us)", to_string(SHDdtmerr)});
    params.push_back({"dead time (us)", to_string(SHDdeadtm)});
    params.push_back({"error dead time (us)", to_string(SHDdtmerr)});
    params.push_back({"error kinectic factor (%)", to_string(SHDerrkf)});
    params.push_back({"error electron current (%)", to_string(SHDerrec)});
    params.push_back({"error ion current (%)", to_string(SHDerric)});
    params.push_back({"error channel width (%)", to_string(SHDerrcw)});
    params.push_back({"error detection efficiency (%)", to_string(SHDerrde)});

    params.push_back({"length of single spectrum", to_string(SHDslen)});
    params.push_back({"run time", to_string(SHDruntim)});
    params.push_back({"real time", to_string(SHDrltcnt)});
    params.push_back({"lifetime", to_string(SHDlftcnt)});
    params.push_back({"processed positions", to_string(SHDdatcnt)});
    params.push_back({"positions out of range", to_string(SHDoutcnt)});
    params.push_back({"processed counter 1 data", to_string(SHDct1cnt)});
    params.push_back({"processed counter 2 data", to_string(SHDct2cnt)});
    params.push_back({"processed counter 3 data", to_string(SHDct3cnt)});
    params.push_back({"processed counter 4 data", to_string(SHDct4cnt)});
    params.push_back({"fifo full counter", to_string(SHDfulcnt)});
    params.push_back({"rejected data", to_string(SHDrejcnt)});
    params.push_back({"sequence errors", to_string(SHDseqcnt)});
    params.push_back({"error counter", to_string(SHDerrcnt)});

    ecur_range = EcurRange(EXPeconvR);
    icur_range = IcurRange(EXPiconvR);
    time_base = 1e6 * pow(2, -round(EXPtimbas)); // in Hz
  }

  /////////////////////////////////////////////////// ESS_2021
  ////////////////////////////////
  else if (data_format == ESS_2021) {
    read_uint16(strzfile, byteorder); // 2 bytes
    read_uint16(strzfile, byteorder); // 2 bytes  for aligment
    uint32_t SHDrltcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDlftcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDdatcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDoutcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct1cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct2cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct3cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDct4cnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDseqcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDbovcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDfulcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDrejcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDerrcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint16_t SHDslen = read_uint16(strzfile, byteorder);   // 2 bytes
    read_uint16(strzfile, byteorder); // 2 bytes for aligment
    read_uint32(strzfile, byteorder); // 4 bytes for aligment
    read_float64(strzfile, byteorder); //  0+1
    float EXPionQ = read_float32(strzfile, byteorder);    //  2
    float EXPionA = read_float32(strzfile, byteorder);    //  3
    float EXPionE = read_float32(strzfile, byteorder);    //  4
    float EXPeff = read_float32(strzfile, byteorder);     //  5
    float EXPtimbas = read_float32(strzfile, byteorder);  //  6
    float EXPeconvR = read_float32(strzfile, byteorder);  //  7
    float EXPeconvFS = read_float32(strzfile, byteorder); //  8
    float EXPiconvR = read_float32(strzfile, byteorder);  //  9
    float EXPiconvFS = read_float32(strzfile, byteorder); // 10
    read_float32(strzfile, byteorder);  // 11
    read_float32(strzfile, byteorder);  // 12
    float EXPwwpot = read_float32(strzfile, byteorder);   // 13
    read_float32(strzfile, byteorder);   // 14
    read_uint32(strzfile, byteorder); // 4 bytes for aligment
    gSHDguntyp = read_uint16(strzfile, byteorder);      // 2 bytes
    read_uint16(strzfile, byteorder); // 2 bytes for aligment
    for (int i = 0; i < 10; i++) {
      gSHDgunpar[i] = read_float32(strzfile, byteorder);
    }
    float SHDdeadtm = read_float32(strzfile, byteorder);
    float SHDdtmerr = read_float32(strzfile, byteorder);
    uint32_t SHDruntim = read_uint32(strzfile, byteorder);
    char ECFid[3];
    strzfile.read(ECFid, 2);
    ECFid[2] = '\0';
    read_uint16(strzfile, byteorder); // 2 bytes (alignment, skipped)
    read_uint32(strzfile, byteorder); // 4 bytes (alignment, skipped)
    uint32_t ECFstep1 = read_uint32(strzfile, byteorder);
    uint32_t ECFstep2 = read_uint32(strzfile, byteorder);
    double ECFstpsiz1 = read_float64(strzfile, byteorder);
    double ECFstpsiz2 = read_float64(strzfile, byteorder);
    double ECFmine1 = read_float64(strzfile, byteorder);
    double ECFmine2 = read_float64(strzfile, byteorder);
    double ECFmaxe1 = read_float64(strzfile, byteorder);
    double ECFmaxe2 = read_float64(strzfile, byteorder);
    double ECFontime1 = read_float64(strzfile, byteorder);
    double ECFontime2 = read_float64(strzfile, byteorder);
    double ECFofftime1 = read_float64(strzfile, byteorder);
    double ECFofftime2 = read_float64(strzfile, byteorder);
    params.push_back({"ion charge", to_string(EXPionQ)});
    params.push_back({"ion mass (u)", to_string(EXPionA)});
    params.push_back({"ion energy (keV)", to_string(EXPionE)});
    params.push_back({"detection efficiency (%)", to_string(EXPeff)});
    params.push_back(
        {"time base exponent n (2**(-n) MHz)", to_string(EXPtimbas)});
    params.push_back({"ele current range ID", to_string(EXPeconvR)});
    params.push_back(
        {"ele current conversion factor (Hz /f.s.)", to_string(EXPeconvFS)});
    params.push_back({"ion current range ID", to_string(EXPiconvR)});
    params.push_back(
        {"ion current conversion factor (Hz /f.s.)", to_string(EXPiconvFS)});
    params.push_back({"dead time (us)", to_string(SHDdeadtm)});
    params.push_back({"error dead time (us)", to_string(SHDdtmerr)});
    params.push_back(
        {"potential of interaction region (V)", to_string(EXPwwpot)});
    params.push_back({"gun type", to_string(gSHDguntyp)});
    params.push_back({"ECF ID", ECFid});
    params.push_back(
        {"number of energy steps in range 1", to_string(ECFstep1)});
    params.push_back({"energy step width in range 1", to_string(ECFstpsiz1)});
    params.push_back({"min energy of range 1", to_string(ECFmine1)});
    params.push_back({"max energy of range 1", to_string(ECFmaxe1)});
    params.push_back({"dwell time of range 1", to_string(ECFontime1)});
    params.push_back({"wait time of range 1", to_string(ECFofftime1)});
    params.push_back(
        {"number of energy steps in range 2", to_string(ECFstep2)});
    params.push_back({"energy step width in range 2", to_string(ECFstpsiz2)});
    params.push_back({"min energy of range 2", to_string(ECFmine2)});
    params.push_back({"max energy of range 2", to_string(ECFmaxe2)});
    params.push_back({"dwell time of range 2", to_string(ECFontime2)});
    params.push_back({"wait time of range 2", to_string(ECFofftime2)});

    params.push_back({"length of single spectrum", to_string(SHDslen)});
    params.push_back({"run time", to_string(SHDruntim)});
    params.push_back({"real time", to_string(SHDrltcnt)});
    params.push_back({"lifetime", to_string(SHDlftcnt)});
    params.push_back({"processed positions", to_string(SHDdatcnt)});
    params.push_back({"positions out of range", to_string(SHDoutcnt)});
    params.push_back({"processed counter 1 data", to_string(SHDct1cnt)});
    params.push_back({"processed counter 2 data", to_string(SHDct2cnt)});
    params.push_back({"processed counter 3 data", to_string(SHDct3cnt)});
    params.push_back({"processed counter 4 data", to_string(SHDct4cnt)});
    params.push_back({"fifo full counter", to_string(SHDfulcnt)});
    params.push_back({"buffer overruns", to_string(SHDbovcnt)});
    params.push_back({"rejected data", to_string(SHDrejcnt)});
    params.push_back({"sequence errors", to_string(SHDseqcnt)});
    params.push_back({"error counter", to_string(SHDerrcnt)});

    ecur_range = EcurRange(EXPeconvR);           // in A
    icur_range = IcurRange(EXPiconvR);           // in A
    time_base = 1e6 * pow(2, -round(EXPtimbas)); // in Hz
  }

  ///////////////////////////////////// MASS_VXW   ||   MASS_VXI
  /////////////////////////////////
  else if ((data_format == MASS_VXW) || (data_format == MASS_VXI)) {
    if (data_format == MASS_VXI) {
      read_uint16(strzfile, byteorder); // 2 bytes
    }
    read_uint16(strzfile, byteorder); // 2 bytes
    uint32_t SHDrltcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDlftcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDdatcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDoutcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDioncnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDtimcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDgaucnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDseqcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDbovcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDrejcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDerrcnt = read_uint32(strzfile, byteorder); // 4 bytes
    uint32_t SHDfulcnt = read_uint32(strzfile, byteorder); // 4 bytes
    if (data_format == MASS_VXW) {
      read_uint32(strzfile, byteorder); // 4 bytes
    }
    read_uint16(strzfile, byteorder); // 2 bytes
    uint16_t SHDslen = read_uint16(strzfile, byteorder);   // 2 bytes
    float EXPblo = read_float32(strzfile, byteorder);      //  0
    float EXPbhi = read_float32(strzfile, byteorder);      //  1
    float EXPuacc = read_float32(strzfile, byteorder);     //  2
    float EXPdiah = read_float32(strzfile, byteorder);     //  3
    float EXPdiav = read_float32(strzfile, byteorder);     //  4
    float EXPfarad = read_float32(strzfile, byteorder);    //  5
    float EXPtimbas = read_float32(strzfile, byteorder);   //  6
    float EXPpress = read_float32(strzfile, byteorder);    //  7
    float EXPbconvFS = read_float32(strzfile, byteorder);  //  8
    float EXPiconvR = read_float32(strzfile, byteorder);   //  9
    float EXPiconvFS = read_float32(strzfile, byteorder);  // 10
    float EXPtstart = read_float32(strzfile, byteorder);   // 11
    float EXPtwait = read_float32(strzfile, byteorder);    // 12
    float EXPtmeas = read_float32(strzfile, byteorder);    // 13
    char SHDgastype[51];
    strzfile.read(SHDgastype, 50);
    SHDgastype[50] = '\0';
    uint32_t SHDruntim = read_uint32(strzfile, byteorder);

    params.push_back({"gas in ion source", SHDgastype});
    params.push_back(
        {"gas pressure (mbar)", to_string(EXPpress * 1E6) + "e-6"});
    params.push_back({"min B field (kG)", to_string(EXPblo)});
    params.push_back({"max B field (kG)", to_string(EXPbhi)});
    params.push_back({"accel. voltage (kV)", to_string(EXPuacc)});
    params.push_back({"diaphragm horizonztal (mm)", to_string(EXPdiah)});
    params.push_back({"diaphragm vertical (mm)", to_string(EXPdiav)});
    params.push_back({"Faraday cup no.", to_string(EXPfarad)});
    params.push_back(
        {"time base exponent n (2**(-n) MHz)", to_string(EXPtimbas)});
    params.push_back(
        {"B-field conversion factor (G /mV)", to_string(EXPbconvFS)});
    params.push_back({"ion current range ID", to_string(EXPiconvR)});
    params.push_back(
        {"ion current conversion factor (Hz /f.s.)", to_string(EXPiconvFS)});
    params.push_back({"waiting time at startup (s)", to_string(EXPtstart)});
    params.push_back(
        {"waiting time before measurement gate (ms)", to_string(EXPtwait)});
    params.push_back(
        {"duration of measurement gate (ms)", to_string(EXPtmeas)});

    params.push_back({"length of single spectrum", to_string(SHDslen)});
    params.push_back({"run time", to_string(SHDruntim)});
    params.push_back({"real time", to_string(SHDrltcnt)});
    params.push_back({"lifetime", to_string(SHDlftcnt)});
    params.push_back({"processed positions", to_string(SHDdatcnt)});
    params.push_back({"positions out of range", to_string(SHDoutcnt)});
    params.push_back({"processed ion data", to_string(SHDioncnt)});
    params.push_back({"processed time data", to_string(SHDtimcnt)});
    params.push_back({"processed Bfield data", to_string(SHDgaucnt)});
    params.push_back({"fifo full counter", to_string(SHDfulcnt)});
    params.push_back({"buffer overruns", to_string(SHDbovcnt)});
    params.push_back({"rejected data", to_string(SHDrejcnt)});
    params.push_back({"sequence errors", to_string(SHDseqcnt)});
    params.push_back({"error counter", to_string(SHDerrcnt)});

    icur_range = IcurRange(EXPiconvR);
    time_base = 1e6 * pow(2, -round(EXPtimbas)); // in Hz
  }
  // ------------------------- end of special headers
  // -------------------------------------
  /////////////////////////////////////////////////////////////////////////////////////////

  // output of electron gun parameters
  if (gSHDguntyp == 4) {
    params.push_back({"MULTI-MODE e-GUN", "3500 eV max"});
    params.push_back({"max energy (eV)", to_string(gSHDgunpar[0])});
    params.push_back({"ramp speed (eV/s)", to_string(gSHDgunpar[1])});
    params.push_back({"P0+P5 (%)", to_string(100 * gSHDgunpar[2])});
    params.push_back({"P1+P4 (%)", to_string(100 * gSHDgunpar[3])});
    params.push_back({"P2+P3 (%)", to_string(100 * gSHDgunpar[4])});
    params.push_back({"P6 (%)", to_string(100 * gSHDgunpar[5])});
    params.push_back({"collector (%)", to_string(100 * gSHDgunpar[6])});
  } else if (gSHDguntyp == 2) {
    params.push_back({"HIGH-INTENSITY e-GUN", "1000 eV max"});
  }

  // ------------------ output of parameters to screen ---------------------
  for (const Parameter &p : params) {
    cout << " " << setw(41) << p.name << ": " << p.value << endl;
  }

  // --------------------- last consistency checks ---------------------------
  if ((nrows == 0) || (ncols == 0)) {
    cout << endl << " No data in input file " << filename << "." << endl;
    return 0;
  }
  if (nbytes != 4) {
    cout << endl
         << " Number of bytes (" << nbytes << ") per datum not supported!";
    return 0;
  }

  // --------------------- output to ascii file ---------------------------
  ofstream fout(ascfilename);
  fout << "####################################################################"
          "#############"
       << endl;
  fout << "### STRZ file converted to ascii by hydrocal" << endl;
  for (const Parameter &p : params) {
    fout << "### " << setw(41) << p.name << ": " << p.value << endl;
  }
  fout << "### " << setw(43)
       << "electron current range (full scale, A): " << ecur_range << endl;
  fout << "### " << setw(43)
       << "ion current range (full scale, A): " << icur_range << endl;
  fout << "### " << setw(43) << "time base (Hz): " << time_base << endl;
  fout << "####################################################################"
          "#############"
       << endl;
  if (abs_data_flag) {
    fout << "   ionization counts [1]";
    fout << "   ion curr (counts) [2]";
    fout << "   ele curr (counts) [3]";
    fout << "   meas time (counts) [4]";
    fout << "   ioniz counts/ion curr. times 2**15 [5]";
    if (nrows > 5)
      fout << "   ele loss (counts)  [6]";
    fout << endl;
  } else if (scan_data_flag) {
    fout << "   ionization counts [1]";
    fout << "   ion curr (counts) [2]";
    fout << "   ele curr (counts) [3]";
    fout << "   meas time (counts) [4]";
    if (nrows > 4)
      fout << "   ele loss (counts)  [5]";
    fout << endl;
  } else if (mass_data_flag) {
    fout << "   ion curr (counts) [1]";
    fout << "   mag field at start (counts) [2]";
    fout << "   mag field averaged (counts) [3]";
    fout << "   meas time (counts) [4]";
    fout << endl;
  }
  fout << "###-----------------------------------------------------------------"
          "-------------"
       << endl;

  if (nbytes == 4) {
    matrix<uint32_t> data(nrows, ncols, 0.0);
    for (unsigned int i = 0; i < nrows; i++) {
      for (unsigned int j = 0; j < ncols; j++) {
        data(i, j) = read_uint32(strzfile, byteorder);
      }
    }
    for (unsigned int j = 0; j < ncols; j++) {
      fout << "  " << data(0, j);
      for (unsigned int i = 1; i < nrows; i++) {
        fout << ",   " << data(i, j);
      }
      fout << endl;
    }
  } // end if (nbytes==4);
  fout.close();
  cout << endl
       << " Data written to ascii file " << ascfilename << "." << endl
       << endl;
  return 0;
}

} // end namespace STRZ
