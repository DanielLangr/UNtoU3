// test_input.cpp - a test driver for UNtoU3 class.
// 
// License: BSD 2-Clause (https://opensource.org/licenses/BSD-2-Clause)
//
// Copyright (c) 2019, Daniel Langr
// All rights reserved.
//
// Program implements the U(N) to U(3) the input irrep [f] specified by the HO level n, N=(n+1)*(n+2)/2, 
// and its number of twos, ones, and zeros read from the standard input.
// For instance, for the input U(21) irrep [f] = [2,2,2,2,2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
// the user should provide the following numbers: 5 6 1 14.
//
// The program performs the U(N) to U(3) reduction and calculates the sum of the dimensions
// of resulting U(3) irrpes multiplied by their level dimensionalities, and print it to the
// standard output. For instance, for the input irrep specified above, the output should read:
// U(3) irreps total dim = 2168999910
//
// This sum should be equal to dim[f], which can be calculated analytically with the support 
// of rational numbers. The program performs this calculcation as well if the Boost library 
// is available and uses its Boost.Rational sublibrary. Availabitliy of Boost is indicated by
// users by definition of HAVE_BOOST preprocessor symbol. 
// For the input irrep [f] specified above, the program should first print out:
// U(N) irrep dim = 2168999910

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>

#ifdef MEASURE_TIME
#include <chrono>
#endif

#ifdef HAVE_BOOST
#include <boost/rational.hpp>
#endif

//#define UNTOU3_DISABLE_TCE
//#define UNTOU3_DISABLE_UNORDERED
//#define UNTOU3_DISABLE_PRECALC
#define UNTOU3_ENABLE_OPENMP
#include "UNtoU3.h"

#ifdef HAVE_BOOST
// Impelments analytical formula for calculation of a dimension of a generic U(N) irrep [f].
// [f] is specified by its labels passed as an array (generaly any indexed data structure)
// argument irrep.
template <typename T>
unsigned long dim(const T & irrep) {
   const auto N = irrep.size();
   boost::rational<unsigned long> result{1};
   for (uint32_t l = 2; l <= N; l++)
      for (uint32_t k = 1; k <= l - 1; k++)
         result *= { irrep[k - 1] - irrep[l - 1] + l - k, l - k };

   assert(result.denominator() == 1);
   return result.numerator();
}
#endif

// Implements analytical formula for calculcation of a dimension of an input U(3) irrep.
// (Does not require rational arithmetics.)
unsigned long dim(const UNtoU3<>::U3Weight & irrep) {
   return (irrep[0] - irrep[1] + 1) * (irrep[0] - irrep[2] + 2) * (irrep[1] - irrep[2] + 1) / 2;
}

int main() {
   // HO level 
   unsigned long n;
   // specification of intput U(N) irrep
   unsigned short n2, n1, n0;
   std::cin >> n >> n2 >> n1 >> n0;

   if (n2 + n1 + n0 != (n + 1) * (n + 2) / 2)
      throw std::invalid_argument("Arguments mismatch!");

#ifdef HAVE_BOOST
   // analytical calculation of dim([f])
   std::vector<unsigned long> f(n2, 2);
   std::fill_n(std::back_inserter(f), n1, 1);
   std::fill_n(std::back_inserter(f), n0, 0);
   std::cout << "U(N) irrep dim = " << dim(f) << std::endl;
#endif

   UNtoU3<> gen;
   // generate HO vectors for a given n
   gen.generateXYZ(n);

#ifdef MEASURE_TIME
   auto start = std::chrono::high_resolution_clock::now();
#endif

   // generation of U(3) irreps in the input U(N) irrep [f]
   gen.generateU3Weights(n2, n1, n0);

#ifdef MEASURE_TIME
   auto end = std::chrono::high_resolution_clock::now();
   std::cout << "U3 weights generation time: "
      << std::chrono::duration<double>(end - start).count()
      << " [s]" << std::endl;
#endif

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
   std::cout << "U(3) irreps total dim = " << sum << std::endl;
}
