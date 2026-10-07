/**
 * @file storagtime.cxx
 *
 * @brief Calculation of the beam lifetime in an ion storage ring
 *
 * see M. Grieser et al., Eur. Phys. J. Special Topics 207 (2012) 1117,
section 4.8
 *
 * This code is based on Manfred Grieser's program LEBENSDAUER.
 *
// SPDX-License-Identifier: MIT
 */

#include "RRDRratecoef.h"
#include "hydroconst.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include "stdin_guard.h"

using namespace std;

/////////////////////////////////////////////////////////////////777
/**
 * @brief Calculates the space charge potential of an electron beam
 *
 * @param  Ee	electron energy in eV
 * @param  Ie	electron current in A
 * @param  rr	ratio of radii (beam tube / electron beam)
 */
double SpaceCharge(double Ee, double Ie, double rr) {
  const int maxiter = 10;
  const double eps = 1.E-5;
  const double me = hydroconst::mec2_eV;

  double pi_eps0_c =
      hydroconst::pi * hydroconst::eps0_As_Vm * hydroconst::clight_m_s;

  double factor = 0.25 * Ie * (1.0 + 2.0 * log(rr)) / pi_eps0_c / me;
  double Uk = Ee / me;
  double enew = 0.92 * Uk;
  double acc = 1.0;
  int i = 0;
  do {
    i++;
    double eold = enew;
    enew = Uk - factor * (eold + 1.0) / sqrt(eold * (eold + 2.0));
    acc = fabs(eold / enew - 1.0);
    if (i > maxiter) {
      printf("Warning:: max number of iterations exceeded!\n");
      printf("          accuracy now :: %8.2e\n", acc);
      break;
    }
  } while (acc > eps);

  return Ee - enew * me;
}

void beam_storage_time(void) {
  const int Nelem = 101; // max Z (+1) of elements considered

  double B[Nelem][Nelem] = {0}; // atomic electron binding energies (source:
                                // NIST ASD ionization energies)
  B[1][1] = 13.599;
  B[2][1] = 24.588;
  B[2][2] = 54.418;
  B[3][1] = 5.392;
  B[3][2] = 75.641;
  B[3][3] = 122.455;
  B[4][1] = 9.323;
  B[4][2] = 18.211;
  B[4][3] = 153.896;
  B[4][4] = 217.72;
  B[5][1] = 8.298;
  B[5][2] = 25.155;
  B[5][3] = 37.931;
  B[5][4] = 259.374;
  B[5][5] = 340.228;
  B[6][1] = 11.26;
  B[6][2] = 24.384;
  B[6][3] = 47.888;
  B[6][4] = 64.494;
  B[6][5] = 392.091;
  B[6][6] = 490.00;
  B[7][1] = 14.534;
  B[7][2] = 29.602;
  B[7][3] = 47.450;
  B[7][4] = 77.474;
  B[7][5] = 97.891;
  B[7][6] = 552.12;
  B[7][7] = 667.05;
  B[8][1] = 13.618;
  B[8][2] = 35.118;
  B[8][3] = 54.936;
  B[8][4] = 77.414;
  B[8][5] = 113.900;
  B[8][6] = 138.12;
  B[8][7] = 739.30;
  B[8][8] = 871.42;
  B[9][1] = 17.423;
  B[9][2] = 34.971;
  B[9][3] = 62.709;
  B[9][4] = 87.141;
  B[9][5] = 114.244;
  B[9][6] = 157.17;
  B[9][7] = 185.19;
  B[9][8] = 953.94;
  B[9][9] = 1103.13;

  B[10][1] = 21.565;
  B[10][2] = 40.963;
  B[10][3] = 63.46;
  B[10][4] = 97.12;
  B[10][5] = 126.22;
  B[10][6] = 157.93;
  B[10][7] = 207.3;
  B[10][8] = 239.09;
  B[10][9] = 1195.9;
  B[10][10] = 1362.21;
  B[11][1] = 5.139;
  B[11][2] = 47.287;
  B[11][3] = 71.621;
  B[11][4] = 98.92;
  B[11][5] = 138.39;
  B[11][6] = 172.15;
  B[11][7] = 208.48;
  B[11][8] = 264.19;
  B[11][9] = 299.88;
  B[11][10] = 1465.14;
  B[11][11] = 1648.71;
  B[12][1] = 7.646;
  B[12][2] = 15.035;
  B[12][3] = 80.144;
  B[12][4] = 109.266;
  B[12][5] = 141.27;
  B[12][6] = 186.51;
  B[12][7] = 224.95;
  B[12][8] = 265.96;
  B[12][9] = 328.24;
  B[12][10] = 367.54;
  B[12][11] = 176.86;
  B[12][12] = 1962.68;

  B[13][1] = 5.986;
  B[13][2] = 18.829;
  B[13][3] = 28.448;
  B[13][4] = 119.994;
  B[13][5] = 153.72;
  B[13][6] = 190.48;
  B[13][7] = 241.44;
  B[13][8] = 284.6;
  B[13][9] = 330.11;
  B[13][10] = 399.37;
  B[13][11] = 442.08;
  B[13][12] = 2086.05;
  B[13][13] = 2304.14;

  B[14][1] = 8.152;
  B[14][2] = 16.346;
  B[14][3] = 33.493;
  B[14][4] = 45.142;
  B[14][5] = 166.769;
  B[14][6] = 205.06;
  B[14][7] = 246.53;
  B[14][8] = 303.8;
  B[14][9] = 351.11;
  B[14][10] = 401.38;
  B[14][11] = 476.08;
  B[14][12] = 523.2;
  B[14][13] = 2437.66;
  B[14][14] = 2673.2;

  B[15][1] = 10.487;
  B[15][2] = 19.726;
  B[15][3] = 30.203;
  B[15][4] = 51.444;
  B[15][5] = 65.026;
  B[15][6] = 220.43;
  B[15][7] = 263.23;
  B[15][8] = 309.42;
  B[15][9] = 371.74;
  B[15][10] = 424.51;
  B[15][11] = 479.59;
  B[15][12] = 560.43;
  B[15][13] = 611.87;
  B[15][14] = 2817.04;
  B[15][15] = 3069.97;

  B[16][1] = 10.36;
  B[16][2] = 23.33;
  B[16][3] = 34.83;
  B[16][4] = 47.31;
  B[16][5] = 72.68;
  B[16][6] = 88.054;
  B[16][7] = 280.94;
  B[16][8] = 328.24;
  B[16][9] = 379.11;
  B[16][10] = 447.10;
  B[16][11] = 504.79;
  B[16][12] = 564.67;
  B[16][13] = 651.65;
  B[16][14] = 707.16;
  B[16][15] = 3223.95;
  B[16][16] = 3494.22;

  B[17][1] = 12.968;
  B[17][2] = 23.814;
  B[17][3] = 39.61;
  B[17][4] = 53.61;
  B[17][5] = 67.8;
  B[17][6] = 97.03;
  B[17][7] = 114.20;
  B[17][8] = 348.29;
  B[17][9] = 400.06;
  B[17][10] = 455.63;
  B[17][11] = 529.28;
  B[17][12] = 591.99;
  B[17][13] = 656.71;
  B[17][14] = 749.76;
  B[17][15] = 809.41;
  B[17][16] = 3658.55;
  B[17][17] = 3946.33;

  B[18][1] = 15.76;
  B[18][2] = 27.63;
  B[18][3] = 40.74;
  B[18][4] = 59.81;
  B[18][5] = 75.02;
  B[18][6] = 91.01;
  B[18][7] = 124.32;
  B[18][8] = 143.46;
  B[18][9] = 422.45;
  B[18][10] = 478.69;
  B[18][11] = 538.96;
  B[18][12] = 618.26;
  B[18][13] = 686.11;
  B[18][14] = 755.75;
  B[18][15] = 854.78;
  B[18][16] = 918;
  B[18][17] = 4120.9;
  B[18][18] = 4426.2;

  B[19][1] = 4.341;
  B[19][2] = 31.626;
  B[19][3] = 45.73;
  B[19][4] = 60.91;
  B[19][5] = 82.66;
  B[19][6] = 100.0;
  B[19][7] = 117.56;
  B[19][8] = 154.71;
  B[19][9] = 175.82;
  B[19][10] = 503.45;
  B[19][11] = 564.15;
  B[19][12] = 629.11;
  B[19][13] = 714.04;
  B[19][14] = 787.16;
  B[19][15] = 861.8;
  B[19][16] = 968;
  B[19][17] = 1034.8;
  B[19][18] = 4611.1;
  B[19][19] = 4934.1;

  B[20][1] = 6.113;
  B[20][2] = 11.872;
  B[20][3] = 50.914;
  B[20][4] = 67.10;
  B[20][5] = 84.41;
  B[20][6] = 108.78;
  B[20][7] = 127.7;
  B[20][8] = 147.24;
  B[20][9] = 188.35;
  B[20][10] = 211.28;
  B[20][11] = 591.27;
  B[20][12] = 656.41;
  B[20][13] = 726.06;
  B[20][14] = 816.64;
  B[20][15] = 895.15;
  B[20][16] = 974;
  B[20][17] = 1087;
  B[20][18] = 1158.1;
  B[20][19] = 5129.2;
  B[20][20] = 5469.9;

  B[21][1] = 6.54;
  B[21][2] = 12.80;
  B[21][3] = 24.757;
  B[21][4] = 73.669;
  B[21][5] = 91.66;
  B[21][6] = 111.1;
  B[21][7] = 138.0;
  B[21][8] = 158.7;
  B[21][9] = 180.03;
  B[21][10] = 225.11;
  B[21][11] = 249.84;
  B[21][12] = 685.91;
  B[21][13] = 755.49;
  B[21][14] = 829.82;
  B[21][15] = 926.03;
  B[21][16] = 1009;
  B[21][17] = 1094;
  B[21][18] = 1213;
  B[21][19] = 1288.4;
  B[21][20] = 5675;
  B[21][21] = 6034;

  B[22][1] = 6.82;
  B[22][2] = 13.58;
  B[22][3] = 27.492;
  B[22][4] = 43.267;
  B[22][5] = 99.3;
  B[22][6] = 119.5;
  B[22][7] = 140.8;
  B[22][8] = 170.4;
  B[22][9] = 192.1;
  B[22][10] = 215.92;
  B[22][11] = 265.0;
  B[22][12] = 291.5;
  B[22][13] = 787.84;
  B[22][14] = 863.1;
  B[22][15] = 941.9;
  B[22][16] = 1044;
  B[22][17] = 1131;
  B[22][18] = 1221;
  B[22][19] = 1346;
  B[22][20] = 1425.9;
  B[22][21] = 6249.4;
  B[22][22] = 6626;

  B[23][1] = 6.74;
  B[23][2] = 14.66;
  B[23][3] = 29.311;
  B[23][4] = 46.709;
  B[23][5] = 65.282;
  B[23][6] = 128.13;
  B[23][7] = 150.17;
  B[23][8] = 173.7;
  B[23][9] = 205.8;
  B[23][10] = 230.5;
  B[23][11] = 255.05;
  B[23][12] = 308.26;
  B[23][13] = 336.28;
  B[23][14] = 896.0;
  B[23][15] = 974;
  B[23][16] = 1060;
  B[23][17] = 1168;
  B[23][18] = 1260;
  B[23][19] = 1355;
  B[23][20] = 1486;
  B[23][21] = 1570.3;
  B[23][22] = 6851;
  B[23][23] = 7246;

  B[24][1] = 6.766;
  B[24][2] = 16.5;
  B[24][3] = 30.96;
  B[24][4] = 49.1;
  B[24][5] = 69.46;
  B[24][6] = 90.64;
  B[24][7] = 160.18;
  B[24][8] = 184.7;
  B[24][9] = 209.3;
  B[24][10] = 244.4;
  B[24][11] = 270.8;
  B[24][12] = 298.0;
  B[24][13] = 354.8;
  B[24][14] = 384.17;
  B[24][15] = 1010.6;
  B[24][16] = 1097;
  B[24][17] = 1185;
  B[24][18] = 1299;
  B[24][19] = 1396;
  B[24][20] = 1496;
  B[24][21] = 1634;
  B[24][22] = 1721.4;
  B[24][23] = 7482;
  B[24][24] = 7894.8;

  B[25][1] = 7.437;
  B[25][2] = 15.64;
  B[25][3] = 33.668;
  B[25][4] = 51.2;
  B[25][5] = 72.4;
  B[25][6] = 95;
  B[25][7] = 119.27;
  B[25][8] = 194.5;
  B[25][9] = 221.8;
  B[25][10] = 248.3;
  B[25][11] = 286.0;
  B[25][12] = 314.4;
  B[25][13] = 343.6;
  B[25][14] = 403.0;
  B[25][15] = 435.6;
  B[25][16] = 1133.1;
  B[25][17] = 1244;
  B[25][18] = 1317;
  B[25][19] = 1437;
  B[25][20] = 1539;
  B[25][21] = 1644;
  B[25][22] = 1788;
  B[25][23] = 1880;
  B[25][24] = 8141.4;
  B[25][25] = 8671.9;
  B[26][1] = 7.871;
  B[26][2] = 16.183;
  B[26][3] = 30.652;
  B[26][4] = 54.8;
  B[26][5] = 75.0;
  B[26][6] = 99.1;
  B[26][7] = 125;
  B[26][8] = 151.06;
  B[26][9] = 233.6;
  B[26][10] = 262.1;
  B[26][11] = 290.3;

  B[26][12] = 330.8;
  B[26][13] = 361.0;
  B[26][14] = 392.2;
  B[26][15] = 457.0;
  B[26][16] = 489.26;
  B[26][17] = 1262;
  B[26][18] = 1358;
  B[26][19] = 1456;
  B[26][20] = 1582;
  B[26][21] = 1689;
  B[26][22] = 1799;
  B[26][23] = 1950;
  B[26][24] = 2045;
  B[26][25] = 8828;
  B[26][26] = 9277.7;

  B[27][1] = 7.86;
  B[27][2] = 17.083;
  B[27][3] = 33.5;
  B[27][4] = 51.3;
  B[27][5] = 79.5;
  B[27][6] = 102;
  B[27][7] = 129;
  B[27][8] = 157;
  B[27][9] = 186.14;
  B[27][10] = 276;
  B[27][11] = 305;
  B[27][12] = 336;
  B[27][13] = 379;
  B[27][14] = 411;
  B[27][15] = 444;
  B[27][16] = 512;
  B[27][17] = 546.8;
  B[27][18] = 1397;
  B[27][19] = 1500;
  B[27][20] = 1603;
  B[27][21] = 1735;
  B[27][22] = 1846;
  B[27][23] = 1962;
  B[27][24] = 2119;
  B[27][25] = 2218;
  B[27][26] = 9544;

  B[28][1] = 7.635;
  B[28][2] = 18.169;
  B[28][3] = 35.17;
  B[28][4] = 54.9;
  B[28][5] = 75.5;
  B[28][6] = 108;
  B[28][7] = 133;
  B[28][8] = 162;
  B[28][9] = 193;
  B[28][10] = 224.5;
  B[28][11] = 321.2;
  B[28][12] = 352;
  B[28][13] = 384;
  B[28][14] = 430;
  B[28][15] = 464;
  B[28][16] = 499;
  B[28][17] = 571;
  B[28][18] = 607.2;
  B[28][19] = 1540;
  B[28][20] = 1648;
  B[28][21] = 1756;
  B[28][22] = 1894;
  B[28][23] = 2011;
  B[28][24] = 2131;
  B[28][25] = 2295;
  B[28][26] = 2399;

  B[29][1] = 7.726;
  B[29][2] = 20.293;
  B[29][3] = 36.84;
  B[29][4] = 55.2;
  B[29][5] = 79.9;
  B[29][6] = 103;
  B[29][7] = 139;
  B[29][8] = 166;
  B[29][9] = 199;
  B[29][10] = 232;
  B[29][11] = 266;
  B[29][12] = 368.9;
  B[29][13] = 401;
  B[29][14] = 435;
  B[29][15] = 484;
  B[29][16] = 520;
  B[29][17] = 557;
  B[29][18] = 633;
  B[29][19] = 671;
  B[29][20] = 1690;
  B[29][21] = 1804;
  B[29][22] = 1916;
  B[29][23] = 2060;
  B[29][24] = 218;
  B[29][25] = 2308;
  B[29][26] = 2478;

  B[30][1] = 9.394;
  B[30][2] = 17.965;
  B[30][3] = 39.724;
  B[30][4] = 59.4;
  B[30][5] = 82.6;
  B[30][6] = 108;
  B[30][7] = 134;
  B[30][8] = 174;
  B[30][9] = 203;
  B[30][10] = 238;
  B[30][11] = 274;
  B[30][12] = 310.8;
  B[30][13] = 419.7;
  B[30][14] = 454;
  B[30][15] = 490;
  B[30][16] = 542;
  B[30][17] = 579;
  B[30][18] = 619;
  B[30][19] = 698;
  B[30][20] = 737;
  B[30][21] = 1846;
  B[30][22] = 1966;
  B[30][23] = 2084;
  B[30][24] = 2234;
  B[30][25] = 2361;
  B[30][26] = 2493;

  B[31][1] = 5.999;
  B[31][2] = 20.51;
  B[31][3] = 30.71;
  B[31][4] = 64;

  B[32][1] = 7.900;
  B[32][2] = 15.935;
  B[32][3] = 34.22;
  B[32][4] = 45.71;
  B[32][5] = 93.5;

  B[33][1] = 9.81;
  B[33][2] = 18.589;
  B[33][3] = 28.352;
  B[33][4] = 50.14;
  B[33][5] = 62.63;
  B[33][6] = 127.6;

  B[34][1] = 9.752;
  B[34][2] = 21.19;
  B[34][3] = 30.821;
  B[34][4] = 42.945;
  B[34][5] = 68.3;
  B[34][6] = 81.81;
  B[34][7] = 155.4;

  B[35][1] = 11.814;
  B[35][2] = 21.8;
  B[35][3] = 36;
  B[35][4] = 47.3;
  B[35][5] = 59.7;
  B[35][6] = 88.6;
  B[35][7] = 103;
  B[35][8] = 192.8;

  B[36][1] = 14.000;
  B[36][2] = 24.360;
  B[36][3] = 36.95;
  B[36][4] = 52.5;
  B[36][5] = 64.7;
  B[36][6] = 78.5;
  B[36][7] = 111.0;
  B[36][8] = 125.94;
  B[36][9] = 230.9;
  B[36][10] = 268;
  B[36][11] = 318;
  B[36][12] = 367;
  B[36][13] = 416;
  ;
  B[36][14] = 465;
  B[36][15] = 515;
  B[36][16] = 565;
  B[36][17] = 591;
  B[36][18] = 641;
  B[36][19] = 786;
  B[36][20] = 833;
  B[36][21] = 884;
  B[36][22] = 936;
  B[36][23] = 997;
  B[36][24] = 1050;
  B[36][25] = 1152;
  B[36][26] = 1205;
  B[36][27] = 2928;
  B[36][28] = 3070;
  B[36][29] = 3227;
  B[36][30] = 3381;
  B[36][31] = 3594;
  B[36][32] = 3760;
  B[36][33] = 3966;
  B[36][34] = 4112;
  B[36][35] = 17310;
  B[36][36] = 17940;

  B[50][1] = 7.344;
  B[50][2] = 14.633;
  B[50][3] = 30.506;
  B[50][4] = 40.74;
  B[50][5] = 94.0;
  B[50][6] = 112.9;
  B[50][7] = 135.0;
  B[50][8] = 156.0;
  B[50][9] = 184.0;
  B[50][10] = 208.0;
  B[50][11] = 232.0;
  B[50][12] = 258.0, B[50][13] = 282.0;
  B[50][14] = 379.0;
  B[50][15] = 407.0;
  B[50][16] = 437.0;
  B[50][17] = 466.0;
  B[50][18] = 506.0;
  B[50][19] = 537.0;
  B[50][20] = 608.0;
  B[50][21] = 642.35;
  B[50][22] = 1127.0;
  B[50][23] = 1195.0;
  B[50][24] = 1269.0,

  B[82][1] = 7.417;
  B[82][2] = 15.03;
  B[82][3] = 31.94;
  B[82][4] = 42.33;
  B[82][5] = 68.8;
  B[82][6] = 82.9;
  B[82][7] = 100.1;
  B[82][8] = 120.0;
  B[82][9] = 138.0;
  B[82][10] = 158.0;
  B[82][11] = 182.0;
  B[82][12] = 203.0;
  B[82][13] = 224;
  B[82][14] = 245.1;
  B[82][15] = 338.1;
  B[82][16] = 374;
  B[82][17] = 401;
  B[82][18] = 427;
  B[82][19] = 478;
  B[82][20] = 507;
  B[82][21] = 570;
  B[82][22] = 610;
  B[82][23] = 650;
  B[82][24] = 690;
  B[82][25] = 750;
  B[82][26] = 810;
  B[82][27] = 870;
  B[82][28] = 930;
  B[82][29] = 990;
  B[82][30] = 1050;
  B[82][31] = 1120;
  B[82][32] = 1180;
  B[82][33] = 1240;
  B[82][34] = 1300;
  B[82][35] = 1360;
  B[82][36] = 1430;
  B[82][37] = 1704;
  B[82][38] = 1760;
  B[82][39] = 1819;
  B[82][40] = 1884;
  B[82][41] = 1945;
  B[82][42] = 2004;
  B[82][43] = 2101;
  B[82][44] = 2163;
  B[82][45] = 2230;
  B[82][46] = 2292;
  B[82][47] = 2543;
  B[82][48] = 2605;
  B[82][49] = 2671;
  B[82][50] = 2735;
  B[82][51] = 2965;
  B[82][52] = 3036;
  B[82][53] = 3211;
  B[82][54] = 3282.1;
  B[82][55] = 5414;
  B[82][56] = 5555;
  B[82][57] = 5703;
  B[82][58] = 5862;
  B[82][59] = 6015;
  B[82][60] = 6162;
  B[82][61] = 6442;
  B[82][62] = 6597;
  B[82][63] = 6767;
  B[82][64] = 6924;
  B[82][65] = 7362;
  B[82][66] = 7500;
  B[82][67] = 7650;
  B[82][68] = 7790;
  B[82][69] = 8520;
  B[82][70] = 8680;
  B[82][71] = 9000;
  B[82][72] = 9150;
  B[82][73] = 19590;
  B[82][74] = 19970;
  B[82][75] = 20380;
  B[82][76] = 20750;
  B[82][77] = 23460;
  B[82][78] = 23940;
  B[82][79] = 24550;
  B[82][80] = 24938;
  B[82][81] = 99492;
  B[82][82] = 101336;

  B[92][1] = 6.194;
  B[92][2] = 11.6;
  B[92][3] = 19.8;
  B[92][4] = 36.7;
  B[92][5] = 46.0;
  B[92][6] = 62.0;
  B[92][7] = 89.0;
  B[92][8] = 101.0;
  B[92][9] = 116.0;
  B[92][10] = 128.9;
  B[92][11] = 158.0;
  B[92][12] = 173.0;
  B[92][13] = 210;
  B[92][14] = 227;
  B[92][15] = 323;
  B[92][16] = 348;
  B[92][17] = 375;
  B[92][18] = 402;
  B[92][19] = 431;
  B[92][20] = 458;
  B[92][21] = 497;
  B[92][22] = 525;
  B[92][23] = 557;
  B[92][24] = 585;
  B[92][25] = 730;
  B[92][26] = 770;
  B[92][27] = 800;
  B[92][28] = 840;
  B[92][29] = 930;
  B[92][30] = 970;
  B[92][31] = 1070;
  B[92][32] = 1110;
  B[92][33] = 1210;
  B[92][34] = 1290;
  B[92][35] = 1370;
  B[92][36] = 1440;
  B[92][37] = 1520;
  B[92][38] = 1590;
  B[92][39] = 1670;
  B[92][40] = 1750;
  B[92][41] = 1830;
  B[92][42] = 1910;
  B[92][43] = 1990;
  B[92][44] = 2070;
  B[92][45] = 2140;
  B[92][46] = 2220;
  B[92][47] = 2578;
  B[92][48] = 2646;
  B[92][49] = 2718;
  B[92][50] = 2794;
  B[92][51] = 2867;
  B[92][52] = 2938;
  B[92][53] = 3073;
  B[92][54] = 3147;
  B[92][55] = 3228;
  B[92][56] = 3301;
  B[92][57] = 3602;
  B[92][58] = 3675;
  B[92][59] = 3753;
  B[92][60] = 3827;
  B[92][61] = 4214;
  B[92][62] = 4299;
  B[92][63] = 4513;
  B[92][64] = 4598;
  B[92][65] = 7393;
  B[92][66] = 7550;
  B[92][67] = 7730;
  B[92][68] = 7910;
  B[92][69] = 8090;
  B[92][70] = 8260;
  B[92][71] = 8650;
  B[92][72] = 8830;
  B[92][73] = 9030;
  B[92][74] = 9210;
  B[92][75] = 9720;
  B[92][76] = 9870;
  B[92][77] = 10040;
  B[92][78] = 10200;
  B[92][79] = 11410;
  B[92][80] = 11600;
  B[92][81] = 11990;
  B[92][82] = 12160;
  B[92][83] = 25260;
  B[92][84] = 25680;
  B[92][85] = 26150;
  B[92][86] = 26590;
  B[92][87] = 31060;
  B[92][88] = 31640;
  B[92][89] = 32400;
  B[92][90] = 32836.5;
  B[92][91] = 129570.3;
  B[92][92] = 131821;

  const double cLight = hydroconst::clight_m_s;
  const double mec2 = hydroconst::mec2_eV;
  const double muc2 = hydroconst::muc2_eV;

  double cooling_energy = 0; // cooling energy in eV
  double ion_A = 48;         // ion mass
  int ion_Z = 22;            // nuclear charge of the stored ion
  int ion_Q = 4;             // ion charge

  // TSR default values
  double ring_acceptance = 0.002;     // ring acceptance in mrad
  double ring_pressure = 5.0E-11;     // pressure of ring vacuum in mbar
  double ring_temperature = 25;       // ring temperature in degree Celsius
  double ring_dipole_rho = 1.15;      // bending radius of dipole magnet in m
  double electron_current = 0.04;     // cooler electron current in A
  double cathode_diameter = 0.009525; // 3/8"
  double expansion_factor = 27.56; // magnetic expansion factor of electron beam
                                   // yielding 5 cm beam diameter
  double cathode_temperature = 0.11; // temperature of cooler cathode in eV
  double longitudinal_temperature =
      0.0002; // longitudinal temperature of the electron beam in eV
  double cooler_length = 1.2;       // length of straight section in cooler in m
  double ring_circumference = 55.4; // circumference of storage ring in m
  double straight_section =
      8.3; // length of straight sections of storage ring in m
  double tube_diameter = 0.1; // inner diameter of vacuum pipe in cooler in m
  double RR_enhancement_factor =
      2.5; // accounts for RR enhancement at very low energies

  double gas_composition[Nelem] = {
      0}; // composition of the residual gas in terms of particle numbers
  gas_composition[1] = 0.9343;  //  H(Z= 1)
  gas_composition[6] = 0.0223;  //  C(Z= 6)
  gas_composition[7] = 0.0075;  //  N(Z= 7)
  gas_composition[8] = 0.0329;  //  0(Z= 8)
  gas_composition[18] = 0.0030; // Ar(Z=18)

  cout << endl
       << " Estimation of the ion-beam lifetime in a heavy-ion storage ring "
       << endl;
  cout << " based on a FORTRAN code by M. Grieser (MPIK Heidelberg)" << endl;
  int ring_selector = 1;
  cout << " Which storage ring? 1: TSR, 2: CRYRING@ESR, 3: ESR, 0: manual "
          "input of ring parameters"
       << endl;
  cout << " Make a selection (0/1/2/3) ..........: ";
  cin >> ring_selector;

  cout << " Give cooling energy in eV ...........: ";
  cin >> cooling_energy;
  cout << " Give ion mass A in amu ..............: ";
  cin >> ion_A;
  cout << " Give nuclear charge Z ...............: ";
  cin >> ion_Z;
  cout << " Give ion charge state Q .............: ";
  cin >> ion_Q;

  if (ring_selector == 1) {       // TSR default values
    ring_dipole_rho = 1.15;       // bending radius of dipole magnet in m
    ring_acceptance = 0.002;      // ring acceptance in mrad
    ring_temperature = 25;        // ring temperature in degree Celsius
    ring_pressure = 5.0E-11;      // pressure of ring vacuum in mbar
    gas_composition[1] = 0.9343;  //  H(Z= 1)
    gas_composition[6] = 0.0223;  //  C(Z= 6)
    gas_composition[7] = 0.0075;  //  N(Z= 7)
    gas_composition[8] = 0.0329;  //  0(Z= 8)
    gas_composition[18] = 0.0030; // Ar(Z=18)
    electron_current = 0.04;      // cooler electron current in A
    cathode_diameter = 0.009525;  // 3/8"
    tube_diameter = 0.1;          // inner diameter of vacuum pipe in m
    expansion_factor = 25.0;      // magnetic expansion factor of electron beam
    cathode_temperature = 0.11;   // temperature of cooler cathode in eV
    longitudinal_temperature =
        0.0002;          // longitudinal temperature of the electron beam in eV
    cooler_length = 1.2; // length of cooler (without toroids) in m
    ring_circumference = 55.4; // circumference of storage ring in m
    straight_section =
        11.8; // length of ecool straight section in m (approximate value)
    RR_enhancement_factor =
        2.5; // accounts for RR enhancement at very low energies
  } else if (ring_selector == 2) { // CRYRING@ESR default values
    ring_dipole_rho = 1.2;         // bending radius of dipole magnet
    ring_acceptance = 0.002;       // ring acceptance in mrad
    ring_temperature = 25;         // ring temperature in degree Celsius
    ring_pressure = 3.0E-11;       // pressure of ring vacuum in mbar
    gas_composition[1] = 0.9343;   //  H(Z= 1)
    gas_composition[6] = 0.0223;   //  C(Z= 6)
    gas_composition[7] = 0.0075;   //  N(Z= 7)
    gas_composition[8] = 0.0329;   //  0(Z= 8)
    gas_composition[18] = 0.0030;  // Ar(Z=18)
    electron_current = 0.01225;    // cooler electron current in A
    cathode_diameter = 0.004;      // 4 mm
    tube_diameter = 0.1;           // inner diameter of vacuum pipe in m
    expansion_factor = 33.0;       // magnetic expansion factor of electron beam
    cathode_temperature = 0.11;    // temperature of cooler cathode in eV
    longitudinal_temperature =
        0.0001;          // longitudinal temperature of the electron beam in eV
    cooler_length = 0.9; // length of cooler (without toroids) in m
    ring_circumference = 54.178; // circumference of storage ring in m
    straight_section =
        4.5; // length of ecool straight section in m (approximate value)
    RR_enhancement_factor =
        2.5; // accounts for RR enhancement at very low energies
  } else if (ring_selector == 3) { // ESR default values
    ring_dipole_rho = 6.25;        // bending radius of dipole magnet
    ring_acceptance = 0.002;       // ring acceptance in mrad
    ring_temperature = 25;         // ring temperature in degree Celsius
    ring_pressure = 3.0E-11;       // pressure of ring vacuum in mbar
    gas_composition[1] = 0.9343;   //  H(Z= 1)
    gas_composition[6] = 0.0223;   //  C(Z= 6)
    gas_composition[7] = 0.0075;   //  N(Z= 7)
    gas_composition[8] = 0.0329;   //  0(Z= 8)
    gas_composition[18] = 0.0030;  // Ar(Z=18)
    electron_current = 0.2;        // cooler electron current in A
    cathode_diameter = 0.0508;     // 50.8 mm
    tube_diameter = 0.25;          // inner diameter of vacuum pipe in m
    expansion_factor = 1.0;        // magnetic expansion factor of electron beam
    cathode_temperature = 0.11;    // temperature of cooler cathode in eV
    longitudinal_temperature =
        0.0002;          // longitudinal temperature of the electron beam in eV
    cooler_length = 2.5; // length of cooler (without toroids) in m
    ring_circumference = 108.36; // circumference of storage ring in m
    straight_section =
        21.0; // length of ecool straight section in m (approximate value)
    RR_enhancement_factor =
        2.5; // accounts for RR enhancement at very low energies
  } else {
    cout << endl;
    cout << " Give ring acceptance in mrad ........: ";
    cin >> ring_acceptance;
    ring_acceptance *= 0.001; // mrad -> rad
    cout << " Give ring circumference in m ........: ";
    cin >> ring_circumference;
    cout << " Give cooler length in m .............: ";
    cin >> cooler_length;
    cout << " Give length of straight section in m : ";
    cin >> straight_section;
    cout << " Give diameter of cathode in mm ......: ";
    cin >> cathode_diameter;
    cathode_diameter *= 0.001; // mm -> m
    cout << " Give bending radius of dipole in m ..: ";
    cin >> ring_dipole_rho;
  }

  char answer;
  cout << endl;
  cout << " electron current: " << electron_current << " A " << endl;
  cout << " electron current ok? (y/n) ..........: ";
  cin >> answer;
  if ((answer == 'N') || (answer == 'n')) {
    cout << " Give electron current in A ..........: ";
    cin >> electron_current;
  }
  cout << endl;

  cout << " expansion factor: " << expansion_factor << endl;
  cout << " expansion factor ok? (y/n) ..........: ";
  cin >> answer;
  if ((answer == 'N') || (answer == 'n')) {
    cout << " Give expansion factor ...............: ";
    cin >> expansion_factor;
  }
  double transverse_temperature = cathode_temperature / expansion_factor;

  cout << endl;
  cout << " transverse temperature: " << transverse_temperature * 1E3 << " meV"
       << endl;
  cout << " transverse temperature ok? (y/n) ....: ";
  cin >> answer;
  if ((answer == 'N') || (answer == 'n')) {
    cout << " Give transverse temperature in meV ..: ";
    cin >> transverse_temperature;
    transverse_temperature *= 1E-3; // meV -> eV
  }

  cout << endl;
  cout << " longitudinal temperature: " << longitudinal_temperature * 1E3
       << " meV" << endl;
  cout << " longitudinal temperature ok? (y/n) ..: ";
  cin >> answer;
  if ((answer == 'N') || (answer == 'n')) {
    cout << " Give longitudinal temperature in meV : ";
    cin >> longitudinal_temperature;
    longitudinal_temperature *= 1E-3; // meV -> eV
  }

  cout << endl;
  cout << " ring temperature: " << ring_temperature << " °C" << endl;
  cout << " temperature ok? (y/n) ...............: ";
  cin >> answer;
  if ((answer == 'N') || (answer == 'n')) {
    cout << " Give ring temperature in °C .........: ";
    cin >> ring_temperature;
  }

  cout << endl;
  cout << " ring vacuum: " << ring_pressure << " mbar" << endl;
  cout << " pressure ok? (y/n) ..................: ";
  cin >> answer;
  if ((answer == 'N') || (answer == 'n')) {
    cout << " Give pressure in mbar ...............: ";
    cin >> ring_pressure;
  }

  cout << endl;
  cout << " residual-gas composition " << endl;
  cout << "  Z   %" << endl;
  double sum_gas = 0.0;
  for (int n = 1; n < Nelem; n++) {
    if (gas_composition[n] < 1E-9)
      continue;
    sum_gas += gas_composition[n];
    cout << " " << setw(2) << n << "  " << setw(5) << 100 * gas_composition[n]
         << endl;
  }
  cout << " sum: " << sum_gas * 100 << endl;
  cout << " residual-gas composition ok? (y/n) ..: ";
  cin >> answer;
  if ((answer == 'N') || (answer == 'n')) {
  input_of_composition:
    cout << endl;
    for (int n = 0; n < Nelem; n++)
      gas_composition[n] = 0.0;
    sum_gas = 0.0;
    int counter = 0;
    for (;;) {
      counter++;
      int z;
      double fraction;
      cout << " Give Z of element " << setw(2) << counter
           << ": (0 quits) .....: ";
      cin >> z;
      if (z == 0)
        break;
      if (z >= Nelem) {
        cout << "ERROR: Z-value out of range" << endl;
      } else {
        cout << " Give fraction of element " << setw(2) << counter
             << " in % ....: ";
        cin >> fraction;
        gas_composition[z] = 0.01 * fraction;
        sum_gas += fraction;
      }
    }
    if (fabs(sum_gas - 100.0) > 1E-3) {
      cout << endl;
      cout << " ERORR: Sum of all components: " << sum_gas
           << "% != 100%. Start over." << endl;
      goto input_of_composition;
    }

  } // end if (answer ...)

  ring_pressure *= 100.0;     // conversion from mbar to Pascal
  ring_temperature += 273.16; // conversion from degree Celsius to Kelvin

  bool warn_flag = false;
  double ebind[Nelem] = {0};
  int ion_N = ion_Z - ion_Q; // number of electrons on stored ion
  for (int n = 1; n <= ion_N; n++) {
    ebind[n] = B[ion_Z][ion_Z + 1 - n];
    // cout << n << "  " << ebind[n] << endl;
    if (ebind[n] == 0.0)
      warn_flag = true;
  }
  if (warn_flag)
    cout << endl
         << " WARNING:: Not all binding energies found in internal table. "
            "Stripping cross section is wrong!"
         << endl;

  double atom_density =
      2.0 * ring_pressure /
      (hydroconst::kB_J_K *
       ring_temperature); // in m^-3, factor 2 because of 2 atoms per molecule
  double ion_mass = ion_A * muc2;
  double ion_energy = cooling_energy * ion_mass / mec2;
  double ion_gamma = 1.0 + cooling_energy / mec2;
  double ion_beta = sqrt(1.0 - pow(ion_gamma, -2));
  double ion_velocity = cLight * ion_beta; // in m/s
  double rigidity =
      ion_gamma * ion_beta * ion_mass / ion_Q / hydroconst::clight_m_s;
  double ring_dipole_E = ion_mass * ion_gamma * ion_beta * ion_beta / ion_Q /
                         ring_dipole_rho; // electric field in dipole in V / m
  int RR_nmax = (int)sqrt(sqrt(pow(ion_Q, 3) * 5.142E11 / (9 * ring_dipole_E)));

  cout << endl;
  cout << " max n for RR from field ionization: " << RR_nmax << endl;
  cout << " max n ok? (y/n) .....................: ";
  cin >> answer;
  if ((answer == 'N') || (answer == 'n')) {
    cout << " Give max n for RR ...................: ";
    cin >> RR_nmax;
  }

  cout << endl;
  cout << " RR enhancement factor: " << RR_enhancement_factor << endl;
  cout << " factor ok? (y/n) ....................: ";
  cin >> answer;
  if ((answer == 'N') || (answer == 'n')) {
    cout << " Give RR enhancement factor ..........: ";
    cin >> RR_enhancement_factor;
    ;
  }

  //-------------------------------------------------------------------------
  // beam lifetime due to single and multiple scattering fromresidual gas
  double sum_msc = 0.0; // multiple scattering
  double sum_ssc = 0.0; // single scattering
  for (int n = 1; n < Nelem; n++) {
    if (gas_composition[n] < 1E-9)
      continue;
    double Zgas = int(n);
    // cout << n << " " << gas_composition[n] << endl;
    double thetamin =
        2.84E-6 * pow(Zgas, 1.0 / 3.0) / (ion_A * ion_beta * ion_gamma);
    double t = hydroconst::pi * pow(1.5E-18 * Zgas * ion_Q / ion_A, 2.0) *
               14.67 * pow(ring_acceptance, -2.0) * pow(ion_beta, -4.0);
    sum_msc += t * atom_density * gas_composition[n] * ion_velocity *
               log(ring_acceptance / thetamin);
    sum_ssc += t * atom_density * gas_composition[n] * ion_velocity / 3.67;
  }
  double tmsc = 1.0 / sum_msc;
  double tssc = 1.0 / sum_ssc;

  //-------------------------------------------------------------------------
  // beam lifetime due to electron capture in collisions with residual gas
  // particles see Schlachter et al. Phys. Rev. A, 27 (1983) 3372
  double sum_c = 0.0;
  for (int n = 1; n < Nelem; n++) {
    if (gas_composition[n] < 1E-9)
      continue;
    double Zgas = int(n);
    double keV_per_nucleon = 0.001 * ion_energy / ion_A;
    double scaled_energy =
        keV_per_nucleon * pow(ion_Q, -0.7) * pow(Zgas, -1.25);
    double a1 = -0.037 * pow(scaled_energy, 2.2);
    double a2 = -2.44e-5 * pow(scaled_energy, 2.6);
    double scaled_xsec =
        1.1E-12 * pow(scaled_energy, -4.8) * (1.0 - exp(a1)) * (1.0 - exp(a2));
    double xsec = scaled_xsec * sqrt(ion_Q) * pow(Zgas, -1.8); // in m^2

    sum_c += xsec * atom_density * gas_composition[n] * ion_velocity;
  }
  double tc = 1.0 / sum_c;

  //-------------------------------------------------------------------------
  // beam lifetime due to stripping
  const double a0 = hydroconst::a0_m;         // Bohr radius
  const double e0 = hydroconst::Ryd_eV;       // Bohrenergy eV
  const double alpha = hydroconst::alpha;     // fine-structure constant
  double v0 = alpha * hydroconst::clight_m_s; // Bohr velocity in m/s
  double sum_s = 1E-99;                       // to prevent division by zero
  for (int n = 1; n < Nelem; n++) {
    if (gas_composition[n] < 1E-9)
      continue;
    double Zgas = int(n);
    double sum = 0;
    double sum1 = 0;
    for (int k = 1; k <= ion_N; k++) {
      if (ebind[k] < 1E-3)
        continue;
      sum += e0 / ebind[k];
      sum1 += sqrt(e0 / ebind[k]);
    }
    double xsec = 0.0;
    if (ion_Q >= Zgas) {
      xsec = 4.0 * hydroconst::pi * pow(a0 * v0 / ion_velocity, 2) * Zgas *
             (Zgas + 1.0) * sum;
    } else {
      xsec = hydroconst::pi * a0 * a0 * pow(Zgas, 2.0 / 3.0) * v0 /
             ion_velocity * sum1;
    }
    sum_s += xsec * atom_density * gas_composition[n] * ion_velocity;
  }
  double ts = 1.0 / sum_s;

  //-------------------------------------------------------------------------
  // beam lifetimes from collisions with the electrons in the cooler

  double electron_beam_radius = 0.5 * cathode_diameter * sqrt(expansion_factor);
  double electron_beam_area = hydroconst::pi * pow(electron_beam_radius, 2.0);
  double electron_density =
      electron_current /
      (ion_velocity * electron_beam_area * hydroconst::e_As); // in m^-3
  double length_ratio = ring_circumference / cooler_length;

  // estimate for the space charge
  // cooling energy should be electron_energy - space_charge;
  // we should first estimate the electron energy and then calculate the space
  // charge here we just say that the electron is the cooling energy
  double space_charge1 =
      SpaceCharge(cooling_energy, electron_current,
                  0.5 * tube_diameter / electron_beam_radius);
  double space_charge2 =
      SpaceCharge(cooling_energy + space_charge1, electron_current,
                  0.5 * tube_diameter / electron_beam_radius);

  /*
  // for alpha_rec and alpha_coll see EPAC 1988, page 962
  double xxx =
  log(11.32*ion_Q/sqrt(transverse_temperature))+0.14*pow(transverse_temperature/(ion_Q*ion_Q),1.0/3.0);
  double alpha_rec = 3.02E-19*ion_Q*ion_Q/sqrt(transverse_temperature)*xxx;
  double alpha_coll
  = 2.0E-39*electron_density*pow(ion_Q,3.0)*pow(transverse_temperature,-4.5);

  double trec = length_ratio/(electron_density*alpha_rec);
  double tcoll = length_ratio/(electron_density*alpha_coll);
  */

  // alternatively RR with Bethe Salpeter formula
  // first find the lowest not fully occupied subshell and the number of
  // electrons therein
  int n_principal, n_sum = 0;
  for (n_principal = 1; n_principal < ion_Z; n_principal++) {
    n_sum += 2 * n_principal * n_principal;
    if (n_sum > ion_N)
      break;
  }
  double RR_nmin = n_principal; // lowest not fully occupied shell
  n_sum -= 2 * RR_nmin * RR_nmin;
  int l_orbital, l_sum = 0;
  for (l_orbital = 0; l_orbital < RR_nmin; l_orbital++) {
    l_sum += 2 * (2 * l_orbital + 1);
    if (l_sum > (ion_N - n_sum))
      break;
  }
  double RR_lmin = l_orbital; // lowest not fully occupied subshell
  l_sum -= 2 * (2 * RR_lmin + 1);
  double RR_nele =
      ion_N - n_sum -
      l_sum; // number of electron in lowest not fully occupied subshell
  double alphaRRscl = 1.0E-6 * RR_enhancement_factor *
                      alpharrscl_cooler(1.0e-6, longitudinal_temperature,
                                        transverse_temperature, ion_Q, RR_nmin,
                                        RR_nmax, RR_lmin, RR_nele);
  double trec = length_ratio / (electron_density * alphaRRscl);

  // sums of all lifetime limiting processes
  double tecool = 1.0 / (sum_s + sum_ssc + sum_c +
                         1.0 / trec); //+1.0/tcoll); // with ecool (no multiple
                                      //scattering in cooled beam)
  double tring = 1.0 / (sum_s + sum_ssc + sum_msc + sum_c); // without ecool

  cout << endl;
  cout << " ion mass in u ...................................: " << ion_A
       << endl;
  cout << " ion nuclear charge ..............................: " << ion_Z
       << endl;
  cout << " ion charge state ................................: " << ion_Q
       << endl;
  cout << " beam energy in MeV ..............................: "
       << ion_energy * 1E-6 << endl;
  cout << " beam energy in MeV/u ............................: "
       << ion_energy * 1E-6 / ion_A << endl;
  cout << " beam beta .......................................: " << ion_beta
       << endl;
  cout << " beam rigidity in T m.............................: " << rigidity
       << endl;
  cout << " ring circumference in m .........................: "
       << ring_circumference << endl;
  cout << " ring acceptance in mrad .........................: "
       << 1000.0 * ring_acceptance << endl;
  cout << " ring temperature in °C ..........................: "
       << ring_temperature - 273.16 << endl;
  cout << " ring vacuum pressure in mbar ....................: "
       << ring_pressure * 0.01 << endl;
  cout << " residual gas density in m^-3 ....................: "
       << 0.5 * atom_density << endl; // factor 0.5 for 2-atomic molecules
  cout << " residual gas composition " << endl;
  cout << "  Z   %" << endl;
  for (int n = 1; n < Nelem; n++) {
    if (gas_composition[n] < 1E-9)
      continue;
    cout << " " << setw(2) << n << "  " << setw(5) << 100 * gas_composition[n]
         << endl;
  }
  cout << " cooler length in m ..............................: "
       << cooler_length << endl;
  cout << " electron current in A ...........................: "
       << electron_current << endl;
  cout << " magnetic expansion factor .......................: "
       << expansion_factor << endl;
  cout << " electron-beam diameter in mm ....................: "
       << 2000.0 * electron_beam_radius << endl;
  cout << " electron density in cm^-3 .......................: "
       << 1.0E-6 * electron_density << endl;
  cout << " cooler cathode temperature in meV ...............: "
       << 1000.0 * cathode_temperature << endl;
  cout << " cooling energy in eV ............................: "
       << cooling_energy << endl;
  cout << " space_charge (0th iteration) in eV ..............: "
       << space_charge1 << endl;
  cout << " space_charge (1st iteration) in eV ..............: "
       << space_charge2 << endl;

  cout << " transverse temperature in meV ...................: "
       << 1000.0 * transverse_temperature << endl;
  cout << " longitudinal temperature in meV .................: "
       << 1000.0 * longitudinal_temperature << endl;
  cout << " RR nmin .........................................: " << RR_nmin
       << endl;
  cout << " RR lmin .........................................: " << RR_lmin
       << endl;
  cout << " RR nele [no. electrons in subshell (nmin,lmin)] .: " << RR_nele
       << endl;
  cout << " RR nmax .........................................: " << RR_nmax
       << endl;
  cout << " RR enhancement factor ...........................: "
       << RR_enhancement_factor << endl;
  cout << " RR alpha in cm^3 s^-1 ...........................: "
       << alphaRRscl * 1.0E6 << endl;
  cout << endl;
  cout << " beam lifetime due to multiple scattering in s ...: " << tmsc
       << endl;
  cout << " beam lifetime due to single scattering in s .....: " << tssc
       << endl;
  cout << " beam lifetime due to charge capture in s ........: " << tc << endl;
  cout << " beam lifetime due to stripping in s .............: " << ts << endl;
  cout << " beam lifetime due to recombination in cooler in s: " << trec
       << endl;
  // cout << " beam lifetime due to collisions in cooler in s ..: " << tcoll <<
  // endl << endl;
  cout << " beam lifetime with electron cooling in s ........: " << tecool
       << endl;
  cout << " beam lifetime without electron cooling in s .....: " << tring
       << endl;
  cout << endl;

  // detector count rates
  double Nion = 1.0E6;
  /*
  cout << " Give number of stored ions.......................: ";
  cin >> Nion;
  cout << endl;
  */
  double rate_rec_cooler = Nion * electron_density * alphaRRscl *
                           (1.0 - ion_beta * ion_beta) * cooler_length /
                           ring_circumference;
  double rate_rec_gas = sum_c * Nion * straight_section / ring_circumference;
  double rate_ion_gas = sum_s * Nion * straight_section / ring_circumference;
  cout << " detector count rates for 1E6 stored ions and a residual gas "
          "pressure of "
       << ring_pressure * 0.01 << " mbar: " << endl
       << endl;
  cout << "       recombination detector from ecool (Hz) ....: "
       << rate_rec_cooler << endl;
  cout << "     recombination detector from res.gas (Hz) ....: " << rate_rec_gas
       << endl;
  cout << "    total rate on recombination detector (Hz) ....: "
       << rate_rec_gas + rate_rec_cooler << endl
       << endl;
  cout << "        ionization detector from res.gas (Hz) ....: " << rate_ion_gas
       << endl;
  cout << endl;
}
