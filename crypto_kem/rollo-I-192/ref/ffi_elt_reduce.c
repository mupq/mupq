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
  o[0] ^= e[0] ^ (e[2] << 39);
  o[1] ^= e[1] ^ (e[2] >> 25) ^ (e[2] << 13);

   uint64_t tmp = e[2] >> 52;
   o[1] ^= (tmp << 13);
   o[0] ^= (tmp << 39);

   // Reduce bits 89-127 of e[1]
   tmp = (o[1] >> 25);
   o[0] ^= tmp ^ (tmp << 38);
   o[1] ^= (tmp >> 26);
   o[1] &= 0x1FFFFFF;
}
