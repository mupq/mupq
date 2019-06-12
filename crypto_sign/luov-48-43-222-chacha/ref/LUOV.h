#ifndef LUOV_H
#define LUOV_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "parameters.h"
#include "Column.h"
#include "LinearAlgebra.h"
#include "buffer.h"
#include "prng.h"
#include "api.h"
#include "ColumnGenerator.h"

#define PK_SEED(pk) pk
#define PK_Q2(pk) (pk + 32)

#define SIG_SOL(sig) sig
#define SIG_SALT(sig) (sig+ VARS*(FIELD_SIZE/8) )

int luov_keygen(unsigned char *pk, unsigned char *sk);
int luov_sign(unsigned char *sig, size_t *smlen, const unsigned char* document, size_t len, const unsigned char *sk);
int luov_verify(unsigned char* m, size_t *mlen, const unsigned char* sm, size_t smlen, const unsigned char *pk);

#endif
