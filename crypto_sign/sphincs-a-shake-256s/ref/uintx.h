#ifndef UINTX_H
#define UINTX_H

#include <stdint.h>
#include <stdbool.h>

typedef uint64_t uint128_t[2];

bool less_u128(const uint128_t *x, const uint128_t *y);

bool eq_u128(const uint128_t *x, const uint128_t *y);

void add_u128(uint128_t *z, const uint128_t *x, const uint128_t *y);

void sub_u128(uint128_t *z, const uint128_t *x, const uint128_t *y);

void set0_u128(uint128_t *z);

void set1_u128(uint128_t *z);

static const uint128_t UINT128_MAX = {~0ULL, ~0ULL};

typedef uint64_t uint256_t[4];

bool less_u256(const uint256_t *x, const uint256_t *y);

void set0_u256(uint256_t *z);
void set1_u256(uint256_t *z);

void add_u256(uint256_t *z, const uint256_t *x, const uint256_t *y);
void sub_u256(uint256_t *z, const uint256_t *x, const uint256_t *y);

#endif
