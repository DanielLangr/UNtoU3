#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>

#ifdef HAVE_BOOST
#include <boost/rational.hpp>
#endif

//#define UNTOU3_DISABLE_TCE
//#define UNTOU3_DISABLE_UNORDERED
//#define UNTOU3_DISABLE_PRECALC
#define UNTOU3_ENABLE_OPENMP
#include "UNtoU3.h"

#ifdef HAVE_BOOST
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

// special case for U(3) irreps (does not require Boost rational numbers)
unsigned long dim(const UNtoU3<>::U3Weight & irrep) {
   return (irrep[0] - irrep[1] + 1) * (irrep[0] - irrep[2] + 2) * (irrep[1] - irrep[2] + 1) / 2;
}

int main() {
   unsigned long n;
   unsigned short n2, n1, n0;
   std::cin >> n >> n2 >> n1 >> n0;

   if (n2 + n1 + n0 != (n + 1) * (n + 2) / 2)
      throw std::invalid_argument("Arguments mismatch!");

#ifdef HAVE_BOOST
   std::vector<unsigned long> f(n2, 2);
   std::fill_n(std::back_inserter(f), n1, 1);
   std::fill_n(std::back_inserter(f), n0, 0);
   std::cout << "U(N) irrep dim = " << dim(f) << std::endl;
#endif

   UNtoU3<> gen;
   gen.generateXYZ(n);
   gen.generateU3Weights(n2, n1, n0);
   unsigned long sum = 0;
   for (const auto & pair : gen.multMap()) {
      const auto & irrep = pair.first;
      if (auto D_l = gen.getLevelDimensionality(irrep)) 
         sum += D_l * dim(irrep);
   }
   std::cout << "U(3) irreps total dim = " << sum << std::endl;
}
