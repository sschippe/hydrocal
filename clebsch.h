// $Id: clebsch.h 375 2016-02-05 16:40:12Z iamp $

#ifndef hydrocal_clebsch
#define hydrocal_clebsch

double NineJ( double ja1, double ja2, double ja3, 
              double jb1, double jb2, double jb3,
              double jc1, double jc2, double jc3); 

double SixJ(double j1, double j2, double j3, 
            double l1, double l2, double l3);

double ThreeJ(double j1, double j2, double j3, 
              double m1, double m2, double m3);

double CG(double j1, double j2, 
          double m1, double m2, 
          double j, double m);

double Strength(double l1, double l2, double ml1, double dm);

/**
 * @brief interactive calculation Clebsch Gordan coefficients for transformation into basis of Stark states 
 */
void StarkCG(void); 

/**
 * @brief interactive calculation of individual Clebsch Gordan coefficients
 */
void testCG(void);

/**
 * @brief interactive calculation of individual 3J-symbols
 */
void test3J(void);

/**
 * @brief interactive calculation of individual 6J-symbols
 */
void test6J(void);

/**
 * @brief interactive calculation of individual 9J-symbols
 */
void test9J(void);

/**
 * @brief interactive calculation of relative line strengths
 */
void testStrength(void);

#endif
