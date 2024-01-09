#ifndef MQOM_KEYGEN_H
#define MQOM_KEYGEN_H

#include <stdint.h>
#include <stddef.h>
#include "parameters.h"
#include "witness.h"

#define PARAM_PUBLICKEYBYTES (PARAM_INSTANCE_SIZE)
#ifndef SHORT_SECRET_KEY
#define PARAM_SECRETKEYBYTES (PARAM_INSTANCE_SIZE + PARAM_SOL_SIZE)
#else
#define PARAM_SECRETKEYBYTES (PARAM_SEED_SIZE)
#endif

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk);
int crypto_sign_valid_keys(unsigned char *pk, unsigned char *sk);

int mqom_keygen(unsigned char *pk, unsigned char *sk);
int mqom_validate_keys(const unsigned char *pk, const unsigned char *sk);

// MQOM Internal
typedef struct mqom_public_key_t {
    instance_t* inst;
} mqom_public_key_t;

typedef struct mqom_secret_key_t {
    instance_t* inst;
    solution_t* wit;
} mqom_secret_key_t;

int mqom_keygen_internal(mqom_public_key_t* pk, mqom_secret_key_t* sk);
int mqom_validate_keys_internal(const mqom_public_key_t* pk, const mqom_secret_key_t* sk);

int deserialize_public_key(mqom_public_key_t* key, const uint8_t* buf, size_t buflen);
int serialize_public_key(uint8_t* buf, const mqom_public_key_t* key, size_t buflen);

int deserialize_secret_key(mqom_secret_key_t* key, const uint8_t* buf, size_t buflen);
int serialize_secret_key(uint8_t* buf, const mqom_secret_key_t* key, size_t buflen);

int mqom_free_keys_internal(mqom_public_key_t* pk, mqom_secret_key_t* sk);

#endif /* MQOM_KEYGEN_H */
