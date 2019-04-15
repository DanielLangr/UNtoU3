#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <stdexcept>
#include <vector>

#ifdef HAVE_BOOST
#include <boost/rational.hpp>
#endif

namespace U3 {
   using LABELS = std::array<uint32_t, 3>;
   using SPS = std::array<std::vector<uint32_t>, 3>;
   enum { NZ, NX, NY };
};

namespace UN {
   using LABELS = std::vector<uint8_t>;
   using BASIS_STATE_WEIGHT_VECTOR = std::vector<uint8_t>;
   using U3MULT_LIST = std::map<U3::LABELS, uint32_t>;
}

enum { MULT, LM, MU, S2 };

void Weight2U3Label(const UN::BASIS_STATE_WEIGHT_VECTOR& vWeights, const U3::SPS& ShellSPS,
                    U3::LABELS& vU3Labels) {
   vU3Labels[U3::NZ] =
       std::inner_product(ShellSPS[U3::NZ].begin(), ShellSPS[U3::NZ].end(), vWeights.begin(), 0);
   vU3Labels[U3::NX] =
       std::inner_product(ShellSPS[U3::NX].begin(), ShellSPS[U3::NX].end(), vWeights.begin(), 0);
   vU3Labels[U3::NY] =
       std::inner_product(ShellSPS[U3::NY].begin(), ShellSPS[U3::NY].end(), vWeights.begin(), 0);
}

void GenerateU3Labels(const UN::LABELS& vGelfandParentRow, uint32_t uSumGelfandParentRow,
                      const U3::SPS& ShellSPS, UN::BASIS_STATE_WEIGHT_VECTOR& vWeights,
                      UN::U3MULT_LIST& mU3LabelsOccurance) {
   size_t N = vGelfandParentRow.size() - 1;

   std::vector<UN::LABELS> vvAllowedLabels;
   UN::LABELS vGelfandRow(N);
   std::vector<uint32_t> vElemsPerChange(N, 1);  // Fill with 1 cause vElemsPerChange[0] = 1;
   std::vector<size_t> vNLabels(N);

   // evaluate all allowed Gelfand patterns based on a parent Gelfand row
   // (vGelfandParentRow) and store them in vvAllowedLabels
   // iNAllowedCombinations is equal to the number of allowed Gelfand patterns
   uint32_t iNAllowedCombinations = 1;
   for (size_t i = 0; i < N; i++) {
      uint32_t uLabelMin = std::min(vGelfandParentRow[i + 1], vGelfandParentRow[i]);
      uint32_t uLabelMax = std::max(vGelfandParentRow[i + 1], vGelfandParentRow[i]);
      UN::LABELS vLabels(uLabelMax - uLabelMin + 1);
      std::iota(vLabels.begin(), vLabels.end(), uLabelMin);
      vNLabels[i] = vLabels.size();
      iNAllowedCombinations *= vNLabels[i];
      vvAllowedLabels.push_back(std::move(vLabels));
   }

   for (size_t i = 1; i < N; i++) {
      vElemsPerChange[i] = vElemsPerChange[i - 1] *
                           vNLabels[i - 1];  // if i == 0 => vElemsPerChange[0]=1 since constructor
   }

   for (uint32_t index = 0; index < iNAllowedCombinations; index++) {
      for (size_t i = 0; i < N; i++) {
         size_t iElement = (index / vElemsPerChange[i]) % vNLabels[i];
         vGelfandRow[i] = vvAllowedLabels[i][iElement];
      }
      uint32_t uSumGelfandRow = std::accumulate(vGelfandRow.begin(), vGelfandRow.end(), 0);
      vWeights[N] = uSumGelfandParentRow - uSumGelfandRow;
      if (N > 1) {  // this condition is due to n = 0 case when N = 0 and hence one wants to keep
                    // vWeights[0] = uSumGelfandRow;
         GenerateU3Labels(vGelfandRow, uSumGelfandRow, ShellSPS, vWeights, mU3LabelsOccurance);
      } else {
         if (N == 1) {
            vWeights[0] = uSumGelfandRow;
         }
         U3::LABELS vU3Labels = {0, 0, 0};
         Weight2U3Label(vWeights, ShellSPS, vU3Labels);
         mU3LabelsOccurance[vU3Labels] += 1;
      }
   }
}

uint32_t GetMultiplicity(const U3::LABELS u3_labels, const UN::U3MULT_LIST& u3_mult_map) {
   auto u3_mult = u3_mult_map.find(u3_labels);
   assert(u3_mult != u3_mult_map.end());

   uint32_t f1 = u3_labels[0], f2 = u3_labels[1], f3 = u3_labels[2];
   uint32_t mult = u3_mult->second;

   u3_mult = u3_mult_map.find({f1 + 1, f2 + 1, f3 - 2});
   mult += (u3_mult == u3_mult_map.end()) ? 0 : u3_mult->second;

   u3_mult = u3_mult_map.find({f1 + 2, f2 - 1, f3 - 1});
   mult += (u3_mult == u3_mult_map.end()) ? 0 : u3_mult->second;

   u3_mult = u3_mult_map.find({f1 + 2, f2, f3 - 2});
   mult -= (u3_mult == u3_mult_map.end()) ? 0 : u3_mult->second;

   u3_mult = u3_mult_map.find({f1 + 1, f2 - 1, f3});
   mult -= (u3_mult == u3_mult_map.end()) ? 0 : u3_mult->second;

   u3_mult = u3_mult_map.find({f1, f2 + 1, f3 - 1});
   mult -= (u3_mult == u3_mult_map.end()) ? 0 : u3_mult->second;

   return mult;
}

void GenerateU3SPS(int n, U3::SPS& ShellSPS) {
   for (int k = 0; k <= n; k++) {
      uint32_t nz = n - k;
      for (int nx = k; nx >= 0; nx--) {
         ShellSPS[U3::NX].push_back(nx);
         ShellSPS[U3::NY].push_back(n - nz - nx);
         ShellSPS[U3::NZ].push_back(nz);
      }
   }
}

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
unsigned long dim(const U3::LABELS & irrep) {
   return (irrep[0] - irrep[1] + 1) * (irrep[0] - irrep[2] + 2) * (irrep[1] - irrep[2] + 1) / 2;
}


int main() {
   unsigned long n;
   unsigned short n2, n1, n0;
   std::cin >> n >> n2 >> n1 >> n0;
   unsigned long N = (n + 1) * (n + 2) / 2;

   if ((n2 + n1 + n0) != N)
      throw std::invalid_argument("Arguments mismatch!");

   U3::SPS ShellSPS;
   GenerateU3SPS(n, ShellSPS);

   UN::LABELS UNLabels(n2, 2);
   std::fill_n(std::back_inserter(UNLabels), n1, 1);
   std::fill_n(std::back_inserter(UNLabels), n0, 0);

#ifdef HAVE_BOOST
   std::cout << "U(N) irrep dim = " << dim(UNLabels) << std::endl;
#endif

   uint32_t sumUNLabels = std::accumulate(UNLabels.begin(), UNLabels.end(), 0);

   UN::U3MULT_LIST mU3_mult;
   UN::BASIS_STATE_WEIGHT_VECTOR Weight(UNLabels.size());

   GenerateU3Labels(UNLabels, sumUNLabels, ShellSPS, Weight, mU3_mult);

   unsigned long sum = 0;
   for (const auto& u3labels_mult : mU3_mult) {
      U3::LABELS U3Labels(u3labels_mult.first);

      if (U3Labels[U3::NZ] >= U3Labels[U3::NX] && U3Labels[U3::NX] >= U3Labels[U3::NY]) 
         if (uint32_t u3_mult = GetMultiplicity(U3Labels, mU3_mult))
            sum += u3_mult * dim(U3Labels);
   }
   std::cout << "U(3) irreps total dim = " << sum << std::endl;
}
