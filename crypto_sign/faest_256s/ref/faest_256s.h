/*
 *  SPDX-License-Identifier: MIT
 */

#include "faest_defines.h"

#include <stddef.h>
#include <stdint.h>

#ifndef FAEST_256S_H
#define FAEST_256S_H

#ifdef __cplusplus
extern "C" {
#endif

#define FAEST_256S_PUBLIC_KEY_SIZE 64
#define FAEST_256S_PRIVATE_KEY_SIZE 32
#define FAEST_256S_SIGNATURE_SIZE 22100

/* Signature API */

/**
 * Key generation function.
 * Generates a public and private key pair, for the specified parameter set.
 *
 * @param[out] pk         The new public key.
 * @param[out] sk         The new private key.
 *
 * @return Returns 0 for success, or a nonzero value indicating an error.
 */
FAEST_EXPORT int FAEST_CALLING_CONVENTION faest_256s_keygen(uint8_t* pk, uint8_t* sk);

/**
 * Signature function.
 * Signs a message with the given keypair.
 *
 * @param[in] sk      The signer's private key.
 * @param[in] pk      The signer's public key.
 * @param[in] message The message to be signed.
 * @param[in] message_len The length of the message, in bytes.
 * @param[in] rho     Additonal randomness; providing randomness renders the signature non-determinstic
 * @param[in] rho_len Length of rho, in bytes.
 * @param[out] signature A buffer to hold the signature. The specific max number of
 * bytes required for a parameter set is given by 256s_signature_size().
 * @param[in,out] signature_len The length of the provided signature buffer.
 * On success, this is set to the number of bytes written to the signature buffer.
 *
 * @return Returns 0 for success, or a nonzero value indicating an error.
 *
 * @see faest_verify(), faest_keygen(), faest_signature_size()
 */
FAEST_EXPORT int FAEST_CALLING_CONVENTION faest_256s_sign(const uint8_t* sk, const uint8_t* pk, const uint8_t* message, size_t message_len, const uint8_t* rho, size_t rho_len, uint8_t* signature, size_t* signature_len);

/**
 * Verification function.
 * Verifies a signature is valid with respect to a public key and message.
 *
 * @param[in] pk      The signer's public key.
 * @param[in] message The message the signature purpotedly signs.
 * @param[in] message_len The length of the message, in bytes.
 * @param[in] signature The signature to verify.
 * @param[in] signature_len The length of the signature.
 *
 * @return Returns 0 for success, indicating a valid signature, or a nonzero
 * value indicating an error or an invalid signature.
 *
 * @see faest_sign(), faest_keygen()
 */
FAEST_EXPORT int FAEST_CALLING_CONVENTION faest_256s_verify(const uint8_t* pk, const uint8_t* message, size_t message_len, const uint8_t* signature, size_t signature_len);

/**
 * Check that a key pair is valid.
 *
 * @param[in] sk The private key to check
 * @param[in] pk The public key to check
 *
 * @return Returns 0 if the key pair is valid, or a nonzero value indicating an error
 */
FAEST_EXPORT int FAEST_CALLING_CONVENTION faest_256s_validate_keypair(const uint8_t* pk, const uint8_t* sk);

/**
 * Clear data of a private key.
 *
 * @param[out] key The private key to clear
 */
FAEST_EXPORT void FAEST_CALLING_CONVENTION faest_256s_clear_private_key(uint8_t* key);

#ifdef __cplusplus
}
#endif

#endif

// vim: ft=c
