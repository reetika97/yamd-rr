//
// Created by raor on 04.06.24.
//

#ifndef LJ_DIRECT_SUMMATION_H
#define LJ_DIRECT_SUMMATION_H

#include "atoms.h"
double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);
double ekin_direct_summation(Atoms &atoms, double mass);

#endif // LJ_DIRECT_SUMMATION_H
