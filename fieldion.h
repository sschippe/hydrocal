// $Id: fieldion.h 168 2013-04-23 10:13:25Z iamp $
double const Fau = 5.1422082e9; // V/cm
double const Tau = 2.4188843e-17; // s

void testFI (void);

void calc_survival(double dt, double f, double z, int nmax, 
                   double *psurv, int *n_one, int *n_zero);

