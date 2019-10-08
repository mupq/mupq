#include "ffi.h"
#include "ffi_elt.h"

/**
 * \fn void ffi_elt_reduce(ffi_elt o, const ffi_elt_ur e)
 * \brief This function reduces a finite field element
 *
 * This function uses f = X^97 + X^6 + 1
 *
 * \param[out] o Finite field element equal to \f$ e \pmod f \f$
 * \param[in] e Finite field element
 */
void ffi_elt_reduce(ffi_elt o, const ffi_elt_ur e) {
// Reduce e[3] e[2]
  o[0] = 0;
  o[1] = 0;
  uint64_t v2 = 0;
  uint64_t tmp = e[3] << 15;
  uint64_t tmp2 = e[3] >> 40;
  v2 ^= e[2] ^ tmp2 ^ (tmp2 >> 9);
  o[1] ^= e[1] ^ tmp ^ (tmp << 9);
  tmp = v2 << 15;
  tmp2 = v2 >> 40;
  o[1] ^= tmp2 ^ (tmp2 >> 9);
  o[0] ^= e[0] ^ tmp ^ (tmp << 9);

  // Reduce bits 113-127 of e[1]
  tmp = (o[1] >> 49);
  o[0] ^= tmp ^ (tmp << 9);

  // Clear reduced bits in e[1]
  o[1] &= 0x1FFFFFFFFFFFF;
}
