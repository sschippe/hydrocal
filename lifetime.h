// $Id: lifetime.h 168 2013-04-23 10:13:25Z iamp $
// Hydrogenic dipole transition probabilities, lifetimes and branching ratios
// see Bethe and Salpeter
// by Stefan Schippers                                             march 1996

double hydrotrans(double z, int n1, int l1, int n2, int l2);
 // hydrogenic transition probabilities for n1,l2 -> n2,l2  
 // transitions in [1/s]                                     

double hydrolife(double z, int n1, int l1, int n2min=0);
// hydrogenic lifetime of state n,l in seconds

double hydrobranch(double z, int n1, int l1, int n2, int l2);
// hydrogenic branching ratio for a n1,l1 -> n2,l2 transition

void testLifetime(void);
