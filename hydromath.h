// $Id: hydromath.h 375 2016-02-05 16:40:12Z iamp $

void TestMath(void); /**< subprogram for performing function tests */


//_____________________________________________________________________________________________________________________
/**
*@brief Cubic spline interpolating function from the Numerical Recipies
*/
void spline(int n, double *x, double*y, double *y2, double yp1=1E31, double ypn=1E31);

//_____________________________________________________________________________________________________________________
/**
* @brief Cubic spline interpolation from the Numerical Recipies
*/
double splint(double x, int n, double *xa, double *ya, double *y2a);



double e1exp(double x); /**< exponential integral E1(x) */



/** elliptic integrals of the first and second kind 
 *
*/
int elliptic(double m, double* K, double* E); 


/** solver for the cubic equation a*x^3 + b*x^2 + c*x + d = 0 
 * returns the number of real roots and the real roots found 
*/
int cuberoot(double a, double b, double c, double d, 
             double *x1, double *x2, double *x3);

double derf( double x); /**< error function erf(x)   */


double derfc(double x); /**< complementary error function 1-erf(x) */


/** Numerically stable implementation of the error function
 * erf(x+h)-erf(x-h) 
 */
double daerf(double x, double h); 


/** ln(n!) up to n<=nmaxfactorial (calculated with MATHEMATICA) 
 *
*/
const int nmaxfactorial = 1000; 


double lnf(double n);  /**< ln(n!) */ 
