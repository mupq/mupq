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
  uint64_t tmp = e[3] << 27;
  uint64_t tmp2 = e[3] >> 30;
  v2 ^= e[2] ^ (tmp2 >> 6) ^ tmp2 ^ (tmp2 >> 1) ^ (tmp2 >> 7);
  o[1] ^= e[1] ^ (tmp) ^ (tmp << 1) ^ (tmp << 6) ^ (tmp << 7);
  tmp = v2 << 27;
  tmp2 = v2 >> 30;
  o[1] ^= tmp2 ^ (tmp2 >> 1) ^ (tmp2 >> 6) ^ (tmp2 >> 7);
  o[0] ^= e[0] ^ tmp ^ (tmp << 1) ^ (tmp << 6) ^ (tmp << 7);

  // Reduce bits 101-127 of e[1]
  tmp = (o[1] >> 37);
  o[0] ^= tmp ^ (tmp << 1) ^ (tmp << 6) ^ (tmp << 7) ;

  // Clear reduced bits in e[1]
  o[1] &= 0x1FFFFFFFFF;
}
