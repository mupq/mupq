// SPDX-License-Identifier: Apache-2.0 or CC0-1.0
#ifndef COMPAT_H
#define COMPAT_H


#define PQCLEAN_VLA(__t,__x,__s) __t __x[__s]

/*********************************
 * Prevent branching on variable *
 *********************************/

#if defined(__GNUC__) || defined(__clang__)
  // Prevent the compiler from
  //    1) inferring that b is 0/1-valued, and
  //    2) handling the two cases with a branch.
  // This is not necessary when verify.c and kem.c are separate translation
  // units, but we expect that downstream consumers will copy this code and/or
  // change how it is built.
# define PQCLEAN_PREVENT_BRANCH_HACK(b)  __asm__("" : "+r"(b) : /* no inputs */);
#else
# define PQCLEAN_PREVENT_BRANCH_HACK(b)
#endif

#endif
