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
   gen.generateXYZ(5); // n
   gen.generateU3Weights(6, 1, 14);
   unsigned long sum = 0;
   for (const auto & pair : gen.multMap()) {
      const auto & irrep = pair.first;
      if (auto D_l = gen.getLevelDimensionality(irrep)) 
         sum += D_l * dim(irrep);
   }
   std::cout << sum << std::endl;
}
