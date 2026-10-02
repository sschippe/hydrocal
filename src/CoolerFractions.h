// $Id: CoolerFractions.h 2039 2026-07-20 07:57:32Z iamp $
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <vector>

void fractionTable(int iselect);

int readfraction(std::vector<double> &fraction, std::string &header, int nmax);

void ngamma(void);

void setup_batch(void);
