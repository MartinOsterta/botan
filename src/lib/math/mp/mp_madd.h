/*
* Lowest Level MPI Algorithms
* (C) 1999-2008,2013 Jack Lloyd
*     2006 Luca Piccarreta
*     2020 Elektrobit Automotive GmbH
*
* Botan is released under the Simplified BSD License (see license.txt)
*/

#ifndef BOTAN_MP_WORD_MULADD_H_
#define BOTAN_MP_WORD_MULADD_H_

#include <botan/types.h>
#include <botan/internal/mul128.h>

namespace Botan {

#if (BOTAN_MP_WORD_BITS == 32)
  typedef uint64_t dword;
  #define BOTAN_HAS_MP_DWORD

#elif (BOTAN_MP_WORD_BITS == 64)
  #if defined(BOTAN_TARGET_HAS_NATIVE_UINT128)
    typedef uint128_t dword;
    #define BOTAN_HAS_MP_DWORD
  #else
    // No native 128 bit integer type; use mul64x64_128 instead
  #endif

#else
  #error BOTAN_MP_WORD_BITS must be 32 or 64
#endif

#if defined(BOTAN_USE_GCC_INLINE_ASM)

  #if defined(BOTAN_TARGET_ARCH_IS_X86_32) && (BOTAN_MP_WORD_BITS == 32)
    #define BOTAN_MP_USE_X86_32_ASM
  #elif defined(BOTAN_TARGET_ARCH_IS_X86_64) && (BOTAN_MP_WORD_BITS == 64)
    #define BOTAN_MP_USE_X86_64_ASM
  #elif defined(BOTAN_TARGET_ARCH_IS_X86_64) && (BOTAN_MP_WORD_BITS == 64) && defined(BOTAN_USE_GCC_INLINE_ASM)
    #define BOTAN_MP_USE_X86_64_ASM
  #elif defined(BOTAN_TARGET_ARCH_IS_ARM64) && (BOTAN_MP_WORD_BITS == 64) && defined(BOTAN_USE_GCC_INLINE_ASM)
    #define BOTAN_MP_USE_ARM_64_ASM
#endif

/*
* Word Multiply/Add
*/
inline word word_madd2(word a, word b, word* c)
   {
#if defined(BOTAN_MP_USE_X86_32_ASM)
   asm(R"(
      mull %[b]
      addl %[c],%[a]
      adcl $0,%[carry]
      )"
      : [a]"=a"(a), [b]"=rm"(b), [carry]"=&d"(*c)
      : "0"(a), "1"(b), [c]"g"(*c) : "cc");

   return a;

#elif defined(BOTAN_MP_USE_X86_64_ASM)
      asm(R"(
         mulq %[b]
         addq %[c],%[a]
         adcq $0,%[carry]
      )"
      : [a]"=a"(a), [b]"=rm"(b), [carry]"=&d"(*c)
      : "0"(a), "1"(b), [c]"g"(*c) : "cc");

   return a;

#elif defined(BOTAN_MP_USE_ARM_64_ASM)
      word a_hi;
      asm(R"(
         umulh %[a_hi], %[a], %[b]
         mul   %[a], %[a], %[b]
         adds  %[a], %[a], %[carry]
         adc   %[carry], %[a_hi], XZR
      )"
      : [a]"+r"(a), [carry]"+r"(*c), [a_hi]"=&r"(a_hi)
      : [b]"r"(b) : "cc");

   return a;

#elif defined(BOTAN_HAS_MP_DWORD)
   const dword s = static_cast<dword>(a) * b + *c;
   *c = static_cast<word>(s >> BOTAN_MP_WORD_BITS);
   return static_cast<word>(s);
#else
   static_assert(BOTAN_MP_WORD_BITS == 64, "Unexpected word size");

   word hi = 0, lo = 0;

   mul64x64_128(a, b, &lo, &hi);

   lo += *c;
   hi += (lo < *c); // carry?

   *c = hi;
   return lo;
#endif
   }

/*
* Word Multiply/Add
*/
inline word word_madd3(word a, word b, word c, word* d)
   {
#if defined(BOTAN_MP_USE_X86_32_ASM)
   asm(R"(
      mull %[b]

      addl %[c],%[a]
      adcl $0,%[carry]

      addl %[d],%[a]
      adcl $0,%[carry]
      )"
      : [a]"=a"(a), [b]"=rm"(b), [carry]"=&d"(*d)
      : "0"(a), "1"(b), [c]"g"(c), [d]"g"(*d) : "cc");

   return a;

#elif defined(BOTAN_MP_USE_X86_64_ASM)
   asm(R"(
      mulq %[b]
      addq %[c],%[a]
      adcq $0,%[carry]
      addq %[d],%[a]
      adcq $0,%[carry]
      )"
      : [a]"=a"(a), [b]"=rm"(b), [carry]"=&d"(*d)
      : "0"(a), "1"(b), [c]"g"(c), [d]"g"(*d) : "cc");

   return a;

#elif defined(BOTAN_MP_USE_ARM_64_ASM)
   word a_hi;
   asm(R"(
      umulh %[a_hi], %[a], %[b]
      mul  %[a], %[a], %[b]
      adds %[a], %[a], %[c]
      adc  %[a_hi], %[a_hi], XZR
      adds %[a], %[a], %[d]
      adc  %[d], %[a_hi], XZR
   )"
   : [a]"+r"(a), [a_hi]"+r"(a_hi), [d]"+r"(*d)
   : [b]"r"(b), [c]"r"(c) : "cc");

   return a;

#elif defined(BOTAN_HAS_MP_DWORD)
   const dword s = static_cast<dword>(a) * b + c + *d;
   *d = static_cast<word>(s >> BOTAN_MP_WORD_BITS);
   return static_cast<word>(s);
#else
   static_assert(BOTAN_MP_WORD_BITS == 64, "Unexpected word size");

   word hi = 0, lo = 0;

   mul64x64_128(a, b, &lo, &hi);

   lo += c;
   hi += (lo < c); // carry?

   lo += *d;
   hi += (lo < *d); // carry?

   *d = hi;
   return lo;
#endif
   }

}

#endif
