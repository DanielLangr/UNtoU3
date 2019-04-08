// test_6414.cpp - a test driver for UNtoU3 class.
// 
// License: BSD 2-Clause (https://opensource.org/licenses/BSD-2-Clause)
//
// Copyright (c) 2019, Daniel Langr
// All rights reserved.
//
// Program implements the U(N) to U(3) reduction where n=5 (N=21) for the input irrep
// [f] = [2,2,2,2,2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
// which is represented by a number of twos n2=6, number of ones n1=1, and number of zeros n0=14.
// Its dimension is dim[f] = 2168999910.
// 
// After reduction, the program iterates over generated U(3) weights and calculates the sum of
// dimensions of resulting U(3) irrpes multiplied by their level dimensionalities. This sum
// should be equal to dim[f].

#include <iostream>

// #define UNTOU3_DISABLE_TCO
// #define UNTOU3_DISABLE_UNORDERED
// #define UNTOU3_DISABLE_PRECALC
#define UNTOU3_ENABLE_OPENMP
#include "UNtoU3.h"

unsigned long dim(const UNtoU3<>::U3Weight & irrep) {
   return (irrep[0] - irrep[1] + 1) * (irrep[0] - irrep[2] + 2) * (irrep[1] - irrep[2] + 1) / 2;
}

int main() {
   UNtoU3<> gen;
   // n=5 - a given HO level, N = (n+1)*(n+2)/2 = 21 
   gen.generateXYZ(5); 
   // generation of U(3) irreps in the input U(21) irrep [f]
   gen.generateU3Weights(6, 1, 14);
   // calculated sum
   unsigned long sum = 0;
   // iteration over generated U(3) weights
   for (const auto & pair : gen.multMap()) {
      // get U(3) weight lables
      const auto & weight = pair.first;
      // get its level dimensionality if its nonzero and the U(3) weight is a U(3) irrep 
      if (auto D_l = gen.getLevelDimensionality(weight)) 
         // add contribution of this U(3) irrep to the sum
         sum += D_l * dim(weight);
   }
   std::cout << sum << std::endl;
}
