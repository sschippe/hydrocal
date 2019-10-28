// $Id: polari.h 168 2013-04-23 10:13:25Z iamp $
/* depolarization due to fine and hyperfine couplings 
   by Stefan Schippers, Giessen, 1996 */


double G2depol(int l, double s, double i = 0 , double gi = 4.43, 
               double q = 0.14, int n = 3, double z = 1.0,
               double v = 0.0, double lt = 0.0, double lp =0.0);

void testHF(void);  // test calculation of hyperfine splittings

void testPol(void); // test polarization coefficients to be
                    // multipliet with m-selective cross sections

void testDepol(void); // test depolarization coefficients
