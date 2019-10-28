// $Id: osci.h 168 2013-04-23 10:13:25Z iamp $
/****************************************************************************/
/* Hydrogenic dipole oscillator strengths                                   */
/* see Kazem Omidvar and Patricia T. Guimares,                              */
/*     The Astrophysical Journal Suplement Series 73 (1990) 555-602         */
/*                                                                          */
/* Copyright: 1997 by Stefan.E.Schippers@strz.uni-giessen.de                */
/*                                                                          */
/****************************************************************************/


double fosciBB(int n1, int l1, int n2, int l2);
//bound bound oscillator strength for transition n1,l1 -> n2,l2

double fosciBC(int n, int l);

double fosciBC(double e, int n, int l);
//bound continuum oscillator strength per unit Rydberg at energy e (in Ryd)

void testOsciBB(void);

void testOsciBC(void);



