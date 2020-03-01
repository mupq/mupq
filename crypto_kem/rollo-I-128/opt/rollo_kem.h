#ifndef LAKE_H
#define LAKE_H

#include "typedef.h"

typedef struct _SecretKey
{
  pu4 pu4Sk_a;
  pu4 pu4Sk_b;
  pu4 pu4F;
} _SecretKey;

typedef struct _PublicKey
{
  pu4 pu4Pk;
  u1  u1Pk_degree;
} _PublicKey;

#endif
