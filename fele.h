// $Id: fele.h 168 2013-04-23 10:13:25Z iamp $

typedef double (*FELE)(double, double, double, double);

double fecool(double e, double eint, double ktpar, double ktperp);

double felegauss(double x, double x0, double fwhm, double dummy);

double flatmax(double e, double eint, double ktpar, double ktperp);

double trapezoid(double e, double eint, double wb, double wt);

double fmaxwell(double kT, double E, double dummy1, double dummy2);
