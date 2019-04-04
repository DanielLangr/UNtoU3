#include <iostream>
#include "UNtoU3.h"

int main() {
   UNtoU3<> gen;
   gen.generateXYZ(2); // n
   gen.generateU3Weights(1, 4, 1);
   for (const auto & pair : gen.multMap()) {
      const auto & irrep = pair.first;
      if (auto D_l = gen.getLevelDimensionality(irrep))
         std::cout << "[" << irrep[0] << "," << irrep[1] << "," << irrep[2] << "] : " << D_l << "\n";
   }
}
