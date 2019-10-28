// $Id: gaussint.h 168 2013-04-23 10:13:25Z iamp $
void gauss_cof(int n, double* x, double* w, int* flag);
/*
 *  DESCRIPTION:
 *  This function provides the 
 *
 *      abscissae x , x , ..... x   and weights w , w , ..... w
 *                 0   1         n-1             0   1          n-1
 *
 *  for n = 2(1)15 , n = 16(4)24 , n = 32(8)48 and n = 64(16)96
 *  related to Gauss-Legendre integration.
 *  The coefficients are tabulated with 15 decimal digits.
 */

void laguer_cof(int n, double* x, double* w, int *flag);
/*
 *  DESCRIPTION:
 *  This function provides the
 *
 *      abscissae x , x , ..... x   and weights w , w , ..... w
 *                 0   1         n-1             0   1          n-1
 *
 *  for n = 2(1)15  related to Gauss-Laguerre integration.
 */
