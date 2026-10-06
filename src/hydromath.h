/**
 * @file hydromath.h
 *
 * @brief some mathematical functions
 *
// SPDX-License-Identifier: MIT
 */

#include "matrix.h"
#include <vector>

#pragma once

void TestMath();

void spline(int n, std::vector<double> x, std::vector<double> y,
            std::vector<double> &y2, double yp1 = 1E31, double ypn = 1E31);

double splint(double x, int n, std::vector<double> xa, std::vector<double> ya,
              std::vector<double> y2a);

double hypconfl(double a, double c, double x);

double e1exp(double x);

int elliptic(double m, double &K, double &E);

int cuberoot(double a, double b, double c, double d, double &x1, double &x2,
             double &x3);

double derf(double x);

double derfc(double x);

double daerf(double x, double h);

double lngamma(double z);

inline constexpr int nmaxfactorial = 1000;

double lnf(double n);

matrix<double> matrixExp(const matrix<double> &A);

std::vector<double> matrixSolve(const matrix<double> &A,
                                const std::vector<double> &b);
