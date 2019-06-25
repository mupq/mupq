//  r5_addsub.h
//  2019-03-10  Markku-Juhani O. Saarinen <mjos@pqshield.com>
//  Low-level functions for handling ternary vectors.

//  Copyright (c) 2019, PQShield Ltd.

#ifndef _R5_ADDSUB_H_
#define _R5_ADDSUB_H_

#include "r5_parameter_sets.h"

//  generic versions

void r5_modq_addsub_d(modq_t *dst,
    const modq_t *p_add, const modq_t *p_sub);

void r5_modq_addsub3_d(modq_t *dst,
    const modq_t *p_add1, const modq_t *p_sub1,
    const modq_t *p_add2, const modq_t *p_sub2,
    const modq_t *p_add3, const modq_t *p_sub3);

void r5_modq_addsub_perm_nbar_d(modq_t *dst, const uint16_t *perm,
    const modq_t *p_add, const modq_t *p_sub);

void r5_modq_addsub3_perm_nbar_d(modq_t *dst, const uint16_t *perm,
    const modq_t *p_add1, const modq_t *p_sub1,
    const modq_t *p_add2, const modq_t *p_sub2,
    const modq_t *p_add3, const modq_t *p_sub3);

void r5_modp_addsub_mu(modp_t *dst,
    const modp_t *p_add, const modp_t *p_sub);

#endif /* _R5_ADDSUB_H_ */

