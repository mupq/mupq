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
  // Reduce e[2]
  o[0] = 0;
  o[1] = 0;
  uint64_t tmp = e[2] << 45;
  uint64_t tmp2 = e[2] >> 12;
  o[0] ^= e[0] ^ (tmp << 2) ^ tmp ^ (tmp << 4) ^ (tmp << 7);
  o[1] ^= e[1] ^ tmp2 ^ (tmp2 >> 3) ^ (tmp2 >> 5)  ^ (tmp2 >> 7);

  // Reduce bits 83-127 of e[1]
  tmp = (o[1] >> 19);
  o[0] ^= tmp ^ (tmp << 2) ^ (tmp << 4) ^ (tmp << 7);

  // Clear reduced bits in e[1]
  o[1] &= 0x7FFFF;
}
