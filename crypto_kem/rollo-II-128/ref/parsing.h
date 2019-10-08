/** 
 * \file parsing.h
 * \brief Functions to parse secret key, public key and ciphertext of the ROLLO scheme
 */

#ifndef PARSING_H
#define PARSING_H

#include "rolloII_types.h"
#include "ffi_poly.h"
#include "ffi_vspace.h"

void rolloII_secret_key_to_string(unsigned char* skString, const unsigned char* seed);
void rolloII_secret_key_from_string(secretKey* sk, unsigned char* skString);


void rolloII_public_key_to_string(unsigned char* pkString, publicKey* pk);
void rolloII_public_key_from_string(publicKey* pk, const unsigned char* pkString);


void rolloII_ciphertext_to_string(unsigned char* ctString, ciphertext* ct);
void rolloII_ciphertext_from_string(ciphertext* ct, const unsigned char* ctString);

#endif
