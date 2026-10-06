// SPDX-License-Identifier: MIT
#pragma once
// Relativistic (Dirac) hydrogenic transition probabilities, lifetimes and
// branching ratios (counterpart of lifetime.h)

double diractrans(double z, int n1, int kappa1, int n2, int kappa2);
// Dirac E1+E2+M1 transition rate n1,kappa1 -> n2,kappa2 in [1/s]
// (0 if the transition is not allowed)

double diraclife(double z, int n, int kappa);
// Dirac lifetime of state n,kappa in seconds (10 s if no allowed decay)

double diracbranch(double z, int n1, int kappa1, int n2, int kappa2);
// Dirac branching ratio for a n1,kappa1 -> n2,kappa2 transition

void testDiracLifetime(void);
