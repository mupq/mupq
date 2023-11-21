/*
 * Copyright 2023 Carlo Sanna, Javier Verbel, and Floyd Zweydinger.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef SIGN_H
#define SIGN_H

#include <stddef.h>
#include <stdint.h>
#include "config.h"
#include "matrix.h"

/* Generate a random pair of public key 'pk' and secret key 'sk'. */
#define crypto_sign_keypair MIRITH_NAMESPACE(keypair)
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk);

/* Sign the message 'msg' using the secret key 'sk'.
 * 
 * Write the signed message over 'sig_msg', which must be a pointer to
 * a pre-allocated array of 'CRYPTO_BYTES + msg_len' bytes.
 * Write the length of the signed message over 'sig_msg_len'.
 *
 * The signed message is formatted as '[signature][message]'. */
#define crypto_sign MIRITH_NAMESPACETOP
int crypto_sign(uint8_t *sig_msg, size_t *sig_msg_len,
    uint8_t *msg, size_t msg_len,
    uint8_t *sk);

/* Verify the signed message 'sig_msg' using the public key 'pk'.
 * If the verification succeeds, return '0', write the message over 'msg',
 * and write the length of the message over 'msg_len'.
 * If the verification fails, return '-1'. */
#define crypto_sign_open MIRITH_NAMESPACE(open)
int crypto_sign_open(uint8_t *msg, size_t *msg_len,
    uint8_t *sig_msg, size_t sig_msg_len,
    uint8_t *pk);

/* Sign the message 'msg' using the (uncompressed) public key 'M' and
 * secret key 'a', 'K', 'E'.
 * 
 * Write the signed message over 'sig_msg', which must be a pointer to
 * a pre-allocated array of 'CRYPTO_BYTES + msg_len' bytes.
 * Write the length of the signed message over 'sig_msg_len'.
 *
 * The signed message is formatted as '[signature][message]'. */
int crypto_sign_unpacked_keys(uint8_t *sig_msg, size_t *sig_msg_len,
    uint8_t *msg, size_t msg_len,
    ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)],
    ff_t a[matrix_bytes_size(PAR_K, 1)],
    ff_t K[matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    ff_t E[matrix_bytes_size(PAR_M, PAR_N)]);
    
#endif
