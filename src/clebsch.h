// SPDX-License-Identifier: MIT
#pragma once
double NineJ(double ja1, double ja2, double ja3, double jb1, double jb2,
             double jb3, double jc1, double jc2, double jc3);

double SixJ(double j1, double j2, double j3, double l1, double l2, double l3);

double ThreeJ(double j1, double j2, double j3, double m1, double m2, double m3);

double CG(double j1, double j2, double m1, double m2, double j, double m);

double Strength(double l1, double l2, double ml1, double dm);

void StarkCG();
void testCG();
void test3J();
void test6J();
void test9J();
void testStrength();
