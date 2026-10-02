// $Id: fele.h 2037 2026-07-17 15:21:56Z iamp $
// SPDX-License-Identifier: MIT
#pragma once

typedef double (*FELE)(double, double, double, double);

double fecool(double e, double eint, double ktpar, double ktperp);

double felegauss(double x, double x0, double fwhm, double dummy);

double flatmax(double e, double eint, double ktpar, double ktperp);

double trapezoid(double e, double eint, double wb, double wt);

double fmaxwell(double kT, double E, double dummy1, double dummy2);

void test_fele(void);
