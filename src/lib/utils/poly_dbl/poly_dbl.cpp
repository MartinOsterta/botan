/*
* (C) 2017,2018 Jack Lloyd
*
* Botan is released under the Simplified BSD License (see license.txt)
*/

#include <botan/internal/poly_dbl.h>
#include <botan/loadstor.h>
#include <botan/exceptn.h>

namespace Botan {

namespace {

/*
* The minimum weight irreducible binary polynomial of size n
*
* See http://www.hpl.hp.com/techreports/98/HPL-98-135.pdf
*/
enum class MinWeightPolynomial : uint64_t {
   P64   = 0x1B,
   P128  = 0x87,
   P192  = 0x87,
   P256  = 0x425,
   P512  = 0x125,
   P1024 = 0x80043,
};

template<size_t Size, MinWeightPolynomial P>
void poly_double(uint8_t out[], const uint8_t in[])
   {
   static_assert(Size % sizeof(word) == 0, "Valid word size");

   const size_t LIMBS = Size / sizeof(word);
   const size_t SHIFT = sizeof(word)*8 - 1;
   word W[LIMBS];
   word C[LIMBS];
   load_be(W, in, LIMBS);

   const word POLY = static_cast<word>(P);

   for(size_t i = 0; i != LIMBS; ++i)
      {
      C[i] = W[(i+1) % LIMBS] >> SHIFT;
      W[i] <<= 1;
      }

   C[LIMBS-1] *= POLY;

   for(size_t i = 0; i != LIMBS; ++i)
      W[i] ^= C[i];

   copy_out_be(out, LIMBS*sizeof(word), W);
   }

template<size_t Size, MinWeightPolynomial P>
void poly_double_le(uint8_t out[], const uint8_t in[])
   {
   typedef word limb;
   static_assert(Size % sizeof(limb) == 0, "Valid limb size");

   const size_t LIMBS = Size / sizeof(limb);
   const size_t SHIFT = sizeof(word)*8 - 1;
   limb W[LIMBS];
   limb C[LIMBS];
   load_le(W, in, LIMBS);

   for(size_t i = 0; i != LIMBS; ++i)
      {
      C[i] = W[(i-1) % LIMBS] >> SHIFT;
      W[i] <<= 1;
      }

   C[0] *= static_cast<limb>(P);

   for(size_t i = 0; i != LIMBS; ++i)
      W[i] ^= C[i];

   copy_out_le(out, LIMBS*sizeof(word), W);
   }

}

void poly_double_n(uint8_t out[], const uint8_t in[], size_t n)
   {
   switch(n)
      {
      case 8:
         return poly_double<8, MinWeightPolynomial::P64>(out, in);
      case 16:
         return poly_double<16, MinWeightPolynomial::P128>(out, in);
      case 24:
         return poly_double<24, MinWeightPolynomial::P192>(out, in);
      case 32:
         return poly_double<32, MinWeightPolynomial::P256>(out, in);
      case 64:
         return poly_double<64, MinWeightPolynomial::P512>(out, in);
      case 128:
         return poly_double<128, MinWeightPolynomial::P1024>(out, in);
      default:
         throw Invalid_Argument("Unsupported size for poly_double_n");
      }
   }

void poly_double_n_le(uint8_t out[], const uint8_t in[], size_t n)
   {
   switch(n)
      {
      case 8:
         return poly_double_le<8, MinWeightPolynomial::P64>(out, in);
      case 16:
         return poly_double_le<16, MinWeightPolynomial::P128>(out, in);
      case 24:
         return poly_double_le<24, MinWeightPolynomial::P192>(out, in);
      case 32:
         return poly_double_le<32, MinWeightPolynomial::P256>(out, in);
      case 64:
         return poly_double_le<64, MinWeightPolynomial::P512>(out, in);
      case 128:
         return poly_double_le<128, MinWeightPolynomial::P1024>(out, in);
      default:
         throw Invalid_Argument("Unsupported size for poly_double_n_le");
      }
   }

}
