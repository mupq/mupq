/** \file tuov_keypair.h
 *  \brief Formats of key pairs.
 */

#ifndef _TUOV_KEYPAIR_H_
#define _TUOV_KEYPAIR_H_


#include "params.h"

#ifdef  __cplusplus
extern  "C" {
#endif

/* alignment 1 for sturct */
#pragma pack(push,1)


/** \brief public key for classic tuov
 * 
 *  when generating the keys, each Pi are stored in two 
 *  parts Pi1, Pi2, correspond to 
 *  the first m1 equations and the last m - m1 equations.
 *  Pi1 is an m1-batched row-major matrix
 *  Pi2 is an (m-m1)-batched row-major matrix, 
 *  but in tuov_publicmap() we want to work with m-batched matrix, 
 *  so in the last step of generate_keypair/expand_pk I 
 *  will use combine_P() to make each Pi an m-batched matrix, 
 *  and we combine 
 *  P1 P2 P3                new P1  newP2
 *     P5 P6   ------->             newP3
 *        P9
 *  new P1 is m-batched V * V upper triangular matrix, correspond to P1
 *  new P2 is m-batched V * O matrix, correspond to P2, P3
 *  new P3 is m-batched O * O upper triangular matrix, correspond to P5, P6, P9
 */
typedef
struct {
    unsigned char pk[_PK_BYTE];         /* store P1, P2, P3, P5, P6, P9 in order */
} pk_t;


/** \brief secret key for classic tuov
 *  when generating the keys, each Fi are stored in two parts Fi1, Fi2
 *  each Fi and T and S are listed below
 *  (if we are working with gf16, T3 will use m/2 bytes)
 *  S       - col-major matrix of size m1 * m1
 *  T1      - a col-major               (n - m) * 1         matrix
 *  T3      - a col-major               1 * (m - 1)         matrix
 *  T4      - a col-major               (n - m) * (m - 1)   matrix
 *  F1      - an m-batched row-major    (n - m) * (n - m)   upper triangular matrix
 *  F2      - an m-batched row-major    (n - m) * 1         matrix
 *  F3      - an m-batched row-major    (n - m) * (m - 1)   matrix
 *  F5      - an m1-batched row-major   1 * 1               matrix
 *  F6      - an m1-batched row-major   1 * (m - 1)         matrix
 *  but in the last step of generate_keypair/expand_sk I 
 *  will use combine_F() to make the following change
 *  F1      change to a m-batched V * V upper trangular matrix
 *  F2 & F3 change to m V * O row-major matrices, the first O-1 cols of each matrix is F3, the last col is F2
 */
typedef
struct {
    unsigned char sk_seed[LEN_SKSEED];

    unsigned char ST[_SK_T_BYTE+_SK_S_BYTE];        /* store S, T1, T3, T4 in order */
    unsigned char F[_SK_F_BYTE];                    /* store F1, F2, F3, F5, F6 in order */
} sk_t;


/** \brief compressed public key
 * 
 *  P51     - an m1-batched row-major   1 * 1 matrix
 *  P61     - an m1-batched row-major   1 * (m - 1) matrix
 *  P9      - an m-batched row-major    (m - 1) * (m - 1) matrix
 */
typedef
struct {
    unsigned char pk_seed[LEN_PKSEED];
    unsigned char P51[_PK_P5_BYTE >> 1];
    unsigned char P61[_PK_P6_BYTE >> 1];
    unsigned char P9[_PK_P9_BYTE];
} cpk_t;


/** \brief compressed secret key
 * 
 *  compressed secret key
 */
typedef
struct {
    unsigned char sk_seed[LEN_SKSEED];
    unsigned char pk_seed[LEN_PKSEED];
} csk_t;




/* restores alignment */
#pragma pack(pop)

#ifdef  __cplusplus
}
#endif


#endif
