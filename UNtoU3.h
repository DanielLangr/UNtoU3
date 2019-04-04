#ifndef UNTOU3_H
#define UNTOU3_H

/*
 * Define before inclsion of this header file if needed.
 * #define UNTOU3_ENABLE_OPENMP
 * #define UNTOU3_DISABLE_TCE
 * #define UNTOU3_DISABLE_UNORDERED
 * #define UNTOU3_DISABLE_PRECALC
*/

#include <array>
#include <cassert>
#include <cstdint>
#include <vector>

#ifndef UNTOU3_DISABLE_UNORDERED
#include <unordered_map>
#else 
#include <map>
#endif

#ifndef UNTOU3_DISABLE_UNORDERED
template <typename T>
struct array_3_hasher
{
   std::size_t operator()(const std::array<T, 3>& key) const
   {
      return key[0] + (key[1] << 8) + (key[2] << 16);
   }
};
#endif /* UNTOU3_DISABLE_UNORDERED */

template <typename T = uint32_t, typename U = uint32_t>
class UNtoU3 {
   public:
      UNtoU3();
      ~UNtoU3();

      using U3Weight = std::array<T, 3>;

#ifndef UNTOU3_DISABLE_UNORDERED
      using U3MultMap = std::unordered_map<U3Weight, U, array_3_hasher<T>>;
#else 
      using U3MultMap = std::map<U3Weight, U>;
#endif /* U3MultMap */

      using GelfandRow = std::array<uint16_t, 3>;
      using GRT = GelfandRow::value_type;

      enum { NZ, NX, NY };

      void generateXYZ(int);
      void generateU3Weights(uint16_t n2, uint16_t n1, uint16_t n0);
      const U3MultMap& multMap() const { return mult_; }
      U getLevelDimensionality(const U3Weight&) const;

   private:
      std::array<std::vector<uint32_t>, 3> xyz_;
      U3MultMap mult_;

#ifndef UNTOU3_DISABLE_PRECALC
      std::array<std::array<std::array<uint8_t, 4>, 4>, 4> cnt_;
      std::array<std::array<std::array<uint16_t*, 4>, 4>, 4> ptr_;
      std::array<uint16_t, 3 * 45> conts_;
#endif

#ifdef UNTOU3_ENABLE_OPENMP
      static U3MultMap* mult_tl_;
#pragma omp threadprivate(mult_tl_)
#endif 

      void generateU3WeightsRec(GelfandRow gpr, U3Weight pp);

#ifndef UNTOU3_DISABLE_PRECALC
      void init_cnt_ptr();
      void init_conts();
#endif
};

#ifdef UNTOU3_ENABLE_OPENMP
template <typename T, typename U>
typename UNtoU3<T, U>::U3MultMap* UNtoU3<T, U>::mult_tl_;
#endif

template <typename T, typename U>
UNtoU3<T, U>::UNtoU3() 
{
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp parallel
   mult_tl_ = new U3MultMap{};
#endif 

#ifndef UNTOU3_DISABLE_PRECALC
   init_cnt_ptr();
#endif
}

template <typename T, typename U>
UNtoU3<T, U>::~UNtoU3() 
{
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp parallel
   delete mult_tl_;
#endif
}

template <typename T, typename U>
void UNtoU3<T, U>::generateXYZ(int n)
{
   for (auto & v : xyz_) v.clear();

   for (int k = 0; k <= n; k++) {
      uint32_t nz = n - k;
      for (int nx = k; nx >= 0; nx--) {
         xyz_[NX].push_back(nx);
         xyz_[NY].push_back(n - nz - nx);
         xyz_[NZ].push_back(nz);
      }
   }

#ifndef UNTOU3_DISABLE_PRECALC
   init_conts();
#endif
}

template <typename T, typename U>
U UNtoU3<T, U>::getLevelDimensionality(const U3Weight& labels) const
{
   T f1 = labels[0], f2 = labels[1], f3 = labels[2];
   if ((f1 < f2) || (f2 < f3)) return 0;

   auto u3_mult = mult_.find({f1, f2, f3});
   assert(u3_mult != mult_.end());
   auto mult = u3_mult->second;

   u3_mult = mult_.find({f1 + 1, f2 + 1, f3 - 2});
   mult += (u3_mult == mult_.end()) ? 0 : u3_mult->second;

   u3_mult = mult_.find({f1 + 2, f2 - 1, f3 - 1});
   mult += (u3_mult == mult_.end()) ? 0 : u3_mult->second;

   u3_mult = mult_.find({f1 + 2, f2, f3 - 2});
   mult -= (u3_mult == mult_.end()) ? 0 : u3_mult->second;

   u3_mult = mult_.find({f1 + 1, f2 - 1, f3});
   mult -= (u3_mult == mult_.end()) ? 0 : u3_mult->second;

   u3_mult = mult_.find({f1, f2 + 1, f3 - 1});
   mult -= (u3_mult == mult_.end()) ? 0 : u3_mult->second;

   return mult;
}

template <typename T, typename U>
void UNtoU3<T, U>::generateU3Weights(uint16_t n2, uint16_t n1, uint16_t n0)
{
   mult_.clear(); // ???
// for (auto & e : mult_) e.second = 0; // ???
   
#ifdef UNTOU3_ENABLE_OPENMP

#pragma omp parallel
   {
      mult_tl_->clear(); // ???
   // for (auto & e : *mult_tl_) e.second = 0; // ???
#pragma omp single
      generateU3WeightsRec({n2, n1, n0}, {0, 0, 0});
#pragma omp critical
      for (const auto temp : *mult_tl_)
         mult_[temp.first] += temp.second;
   }

#else  /* UNTOU3_ENABLE_OPENMP */

   generateU3WeightsRec({n2, n1, n0}, {0, 0, 0});

#endif /* UNTOU3_ENABLE_OPENMP */
}

template <typename T, typename U>
void UNtoU3<T, U>::generateU3WeightsRec(GelfandRow gpr, U3Weight pp) 
{
   size_t N = gpr[0] + gpr[1] + gpr[2] - 1;

#ifndef UNTOU3_DISABLE_TCE
   while
#else
   if
#endif
   (N > 
#ifndef UNTOU3_DISABLE_PRECALC
   2
#else 
   0
#endif 
   ) {
       if (gpr[0]) {
           if (gpr[1] || gpr[2]) {
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { (GRT)(gpr[0] - 1), gpr[1], gpr[2] }, 
                       { pp[0] + 2 * xyz_[0][N], pp[1] + 2 * xyz_[1][N], pp[2] + 2 * xyz_[2][N] });
               if (gpr[2]) 
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
                   generateU3WeightsRec( { (GRT)(gpr[0] - 1), (GRT)(gpr[1] + 1), (GRT)(gpr[2] - 1) },
                           { pp[0] + xyz_[0][N], pp[1] + xyz_[1][N], pp[2] + xyz_[2][N] });
           }
           else {
#ifndef UNTOU3_DISABLE_TCE
               gpr[0]--;
               pp[0] += 2 * xyz_[0][N]; pp[1] += 2 * xyz_[1][N]; pp[2] += 2 * xyz_[2][N];
#else /* UNTOU3_DISABLE_TCE */
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { (GRT)(gpr[0] - 1), gpr[1], gpr[2] }, 
                       { pp[0] + 2 * xyz_[0][N], pp[1] + 2 * xyz_[1][N], pp[2] + 2 * xyz_[2][N] });
#endif /* UNTOU3_DISABLE_TCE */
           }
       }

       if (gpr[1]) 
           if (gpr[2])
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { gpr[0], (GRT)(gpr[1] - 1), gpr[2] },
                       { pp[0] + xyz_[0][N], pp[1] + xyz_[1][N], pp[2] + xyz_[2][N] });
           else {
#ifndef UNTOU3_DISABLE_TCE
              gpr[1]--; 
              pp[0] += xyz_[0][N]; pp[1] += xyz_[1][N]; pp[2] += xyz_[2][N];
#else /* UNTOU3_DISABLE_TCE */
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { gpr[0], (GRT)(gpr[1] - 1), gpr[2] },
                       { pp[0] + xyz_[0][N], pp[1] + xyz_[1][N], pp[2] + xyz_[2][N] });
#endif /* UNTOU3_DISABLE_TCE */
           }

       if (gpr[2]) {
#ifndef UNTOU3_DISABLE_TCE
           gpr[2]--;
#else /* UNTOU3_DISABLE_TCE */
#ifdef UNTOU3_ENABLE_OPENMP
#pragma omp task if (N > 8) firstprivate(gpr, pp)
#endif
               generateU3WeightsRec( { gpr[0], gpr[1], (GRT)(gpr[2] - 1) }, { pp[0], pp[1], pp[2] }); 
#endif /* UNTOU3_DISABLE_TCE */
       }

#ifndef UNTOU3_DISABLE_TCE
       N--;
#endif
   }
#ifdef UNTOU3_DISABLE_TCE
   else {
#endif

#ifndef UNTOU3_DISABLE_PRECALC

   auto c = cnt_[gpr[0]][gpr[1]][gpr[2]];
   if (c > 0) {
      auto p = ptr_[gpr[0]][gpr[1]][gpr[2]];
      for (int i = 0; i < c; i++) {
         U3Weight pp_;
         pp_[0] = pp[0] + *p++;
         pp_[1] = pp[1] + *p++;
         pp_[2] = pp[2] + *p++;
#ifdef UNTOU3_ENABLE_OPENMP
         (*mult_tl_)[pp_] += 1;
#else
         mult_[pp_] += 1;
#endif
      }
   }
   else  {
#ifdef UNTOU3_ENABLE_OPENMP
      (*mult_tl_)[pp] += 1;
#else
      mult_[pp] += 1;
#endif
   }

#else /* UNTOU3_DISABLE_PRECALC */

   auto temp = 2 * gpr[0] + gpr[1];
   pp[0] += temp * xyz_[0][0];
   pp[1] += temp * xyz_[1][0];
   pp[2] += temp * xyz_[2][0];

#ifdef UNTOU3_ENABLE_OPENMP
   (*mult_tl_)[pp] += 1;
#else
   mult_[pp] += 1;
#endif

#endif /* UNTOU3_DISABLE_PRECALC */

#ifdef UNTOU3_DISABLE_TCE
   }
#endif
}

#ifndef UNTOU3_DISABLE_PRECALC

template <typename T, typename U>
void UNtoU3<T, U>::init_cnt_ptr()
{
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         for (int k = 0; k < 4; k++)
            cnt_[i][j][k] = 0;

   cnt_[3][0][0] = 1; cnt_[0][3][0] = 1;                     //  2
   cnt_[2][1][0] = 3; cnt_[2][0][1] = 6; cnt_[1][2][0] = 3;  // 12
   cnt_[0][2][1] = 3; cnt_[1][0][2] = 6; cnt_[0][1][2] = 3;  // 12
   cnt_[1][1][1] = 8;                                        //  8
   cnt_[2][0][0] = 1; cnt_[0][2][0] = 1;                     //  2
   cnt_[1][1][0] = 2; cnt_[1][0][1] = 3; cnt_[0][1][1] = 2;  //  7
   cnt_[1][0][0] = 1; cnt_[0][1][0] = 1;                     //  2

   auto p = conts_.data();
   ptr_[3][0][0] = p; p += 1 * 3; ptr_[0][3][0] = p; p += 1 * 3;                    
   ptr_[2][1][0] = p; p += 3 * 3; ptr_[2][0][1] = p; p += 6 * 3; ptr_[1][2][0] = p; p += 3 * 3;  
   ptr_[0][2][1] = p; p += 3 * 3; ptr_[1][0][2] = p; p += 6 * 3; ptr_[0][1][2] = p; p += 3 * 3;  
   ptr_[1][1][1] = p; p += 8 * 3;                                      
   ptr_[2][0][0] = p; p += 1 * 3; ptr_[0][2][0] = p; p += 1 * 3;                   
   ptr_[1][1][0] = p; p += 2 * 3; ptr_[1][0][1] = p; p += 3 * 3; ptr_[0][1][1] = p; p += 2 * 3;  
   ptr_[1][0][0] = p; p += 1 * 3; ptr_[0][1][0] = p;                    
}

template <typename T, typename U>
void UNtoU3<T, U>::init_conts() 
{
   auto p = conts_.data();
   // (3,0,0)
   *p++ = 2 * xyz_[0][2] + 2 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 2 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 2 * xyz_[2][1] + 2 * xyz_[2][0];

   // (0,3,0)
   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   // (2,1,0)
   *p++ = 1 * xyz_[0][2] + 2 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 2 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 2 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 1 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 1 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 1 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 2 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 2 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 2 * xyz_[2][1] + 1 * xyz_[2][0];

   // (2,0,1)
   *p++ = 0 * xyz_[0][2] + 2 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 2 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 2 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 2 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 2 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 2 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 0 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 0 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 0 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 2 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 2 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 2 * xyz_[2][1] + 0 * xyz_[2][0];

   // (1,2,0)
   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 2 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 2 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 2 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   // (0,2,1)
   *p++ = 0 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 0 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 0 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 0 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 0 * xyz_[2][0];

   // (1,0,2)
   *p++ = 0 * xyz_[0][2] + 0 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 0 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 0 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 0 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 0 * xyz_[0][2] + 2 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 2 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 2 * xyz_[2][1] + 0 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 0 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 0 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 0 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 0 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 0 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 0 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 0 * xyz_[2][1] + 0 * xyz_[2][0];

   // (0,1,2)
   *p++ = 0 * xyz_[0][2] + 0 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 0 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 0 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 0 * xyz_[0][2] + 1 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 1 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 1 * xyz_[2][1] + 0 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 0 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 0 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 0 * xyz_[2][1] + 0 * xyz_[2][0];

   // (1,1,1)
   *p++ = 0 * xyz_[0][2] + 1 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 1 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 1 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 0 * xyz_[0][2] + 2 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 0 * xyz_[1][2] + 2 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 0 * xyz_[2][2] + 2 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 0 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 0 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 0 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 2 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 2 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 2 * xyz_[2][1] + 0 * xyz_[2][0];

   *p++ = 1 * xyz_[0][2] + 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][2] + 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][2] + 1 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 0 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 0 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 0 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 2 * xyz_[0][2] + 1 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 2 * xyz_[1][2] + 1 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 2 * xyz_[2][2] + 1 * xyz_[2][1] + 0 * xyz_[2][0];

   // (2,0,0)
   *p++ = 2 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 2 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 2 * xyz_[2][1] + 2 * xyz_[2][0];

   // (0,2,0)
   *p++ = 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][1] + 1 * xyz_[2][0];

   // (1,1,0)
   *p++ = 1 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 1 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 1 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 2 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 2 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 2 * xyz_[2][1] + 1 * xyz_[2][0];

   // (1,0,1)
   *p++ = 0 * xyz_[0][1] + 2 * xyz_[0][0];
   *p++ = 0 * xyz_[1][1] + 2 * xyz_[1][0];
   *p++ = 0 * xyz_[2][1] + 2 * xyz_[2][0];

   *p++ = 2 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 2 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 2 * xyz_[2][1] + 0 * xyz_[2][0];

   *p++ = 1 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][1] + 1 * xyz_[2][0];

   // (0,1,1)
   *p++ = 0 * xyz_[0][1] + 1 * xyz_[0][0];
   *p++ = 0 * xyz_[1][1] + 1 * xyz_[1][0];
   *p++ = 0 * xyz_[2][1] + 1 * xyz_[2][0];

   *p++ = 1 * xyz_[0][1] + 0 * xyz_[0][0];
   *p++ = 1 * xyz_[1][1] + 0 * xyz_[1][0];
   *p++ = 1 * xyz_[2][1] + 0 * xyz_[2][0];

   // (1,0,0)
   *p++ = 2 * xyz_[0][0];
   *p++ = 2 * xyz_[1][0];
   *p++ = 2 * xyz_[2][0];

   // (0,1,0)
   *p++ = 1 * xyz_[0][0];
   *p++ = 1 * xyz_[1][0];
   *p++ = 1 * xyz_[2][0];
}

#endif /* UNTOU3_DISABLE_PRECALC */

#endif /* UNTOU3_H */
