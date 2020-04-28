/**
 * \file api.h
 * \brief NIST KEM API used by the ROLLOI_KEM IND-CPA scheme
 */


#ifndef ROLLOI_KEM_H
#define ROLLOI_KEM_H

int rolloI_kem_keygen(uint8_t* pk, uint8_t* sk);
int rolloI_kem_encaps(uint8_t* ct, uint8_t* ss, const uint8_t* pk);
int rolloI_kem_decaps(uint8_t* ss, const uint8_t* ct, const uint8_t* sk);

#endif
