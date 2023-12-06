/** \file   tuov_keypair_computation.c
 *  \brief  Implementations for functions in tuov_keypair.c
 */

#ifndef _TUOV_KEYPAIR_COMPUTATION_H_
#define _TUOV_KEYPAIR_COMPUTATION_H_

/** \brief compute matrix product M = S * M
 *  
 *  \param[in, out] M       - an col-major m1 * D matrix
 *  \param[in] S            - an col-major m1 * m1 matrix 
 *  \param[in] D            - number of cols of M
*/
void multiple_S(unsigned char *M, const unsigned char *S, unsigned D);

/** \brief 
 *  \param[in, out] F        - part of the secret key
*/
void combine_F(unsigned char *F);


/** \brief in the computation of keygen, we store P1, P2, P3... in two parts, i.e. P11, P12, ...
 *          we want to combine two m1-batched matrix Pi1, Pi2 to a m-batched matrix, the 
 *          output P contains 6 m-batched matrix P1, P2, P3, ..., P6 in order.
 *  \param[in, out] P        - part of the public key
*/
void combine_P(unsigned char *P);

/** \brief compute Q51 from T1, Q11, Q21, 
 * 
 * \param[out]  Q51         - an m1-batched row-major 1 * 1 matrix
 * \param[in]   T1          - a col-major (n - m) * 1 matrix
 * \param[in]   Q11         - an m1-batched row-major (n - m) * (n - m) upper triangular matrix
 * \param[in]   Q21         - an m1-batched row-major (n - m) * 1 matrix
*/
void calculate_Q51(unsigned char *Q51, const unsigned char *T1, const unsigned char *Q11, const unsigned char *Q21);

/** \brief compute Q61 from T1, T3, T4, Q11, Q21, Q31, Q51
 *  \param[out] Q61         - an m1-batched row-major 1 * (m - 1) matrix
 *  \param[in] T1           - a col-major (n - m) * 1 matrix
 *  \param[in] T3           - a col-major 1 * (m - 1) matrix
 *  \param[in] T4           - a col-major (n - m) * (m - 1) matrix
 *  \param[in] Q11          - an m1-batched row-major (n - m) * (n - m) matrix
 *  \param[in] Q21          - an m1-batched row-major (n - m) * 1 matrix
 *  \param[in] Q31          - an m1-batched row-major (n - m) * (m - 1) matrix
 *  \param[in] Q51          - an m1-batched row-major 1 * 1 matrix
*/
void calculate_Q61(unsigned char *Q61, const unsigned char *T1, const unsigned char *T3, const unsigned char *T4, 
                    const unsigned char *Q11, const unsigned char *Q21, const unsigned char *Q31, const unsigned char *Q51);

/** \brief compute Q91 from T3, T4, Q11, Q21, Q31, Q51, Q61
 *  \param[out] Q91         - an m1-batched row-major (m - 1) * (m - 1) upper triangular matrix
 *  \param[in] T3           - a col-major 1 * (m - 1) matrix
 *  \param[in] T4           - a col-major (n - m) * (m - 1) matrix
 *  \param[in] Q11          - an m1-batched row-major (n - m) * (n - m) matrix
 *  \param[in] Q21          - an m1-batched row-major (n - m) * 1 matrix
 *  \param[in] Q31          - an m1-batched row-major (n - m) * (m - 1) matrix
 *  \param[in] Q51          - an m1-batched row-major 1 * 1 matrix 
 *  \param[in] Q61          - an m1-batched row-major 1 * (m - 1) matrix
 */ 
void calculate_Q91(unsigned char *Q91, const unsigned char *T3, const unsigned char *T4, const unsigned char *Q11,
                    const unsigned char *Q21, const unsigned char *Q31, const unsigned char *Q51, const unsigned char *Q61);


/** \brief compute F21 from T1, Q11, Q21
 *  \param F21              - an m1-batched row-major (n - m) * 1 matrix
 *  \param[in] T1           - a col-major (n - m) * 1 matrix
 *  \param Q11              - an m1-batched row-major (n - m) * (n - m) matrix
 *  \param Q21              - an m1-batched row-major (n - m) * 1 matrix
 */ 
void calculate_F21(unsigned char *F21, const unsigned char *T1, const unsigned char *Q11, const unsigned char *Q21);

/** \brief compute F31 from T3, T4, Q11, Q21, Q31
 *  \param F31 
 *  \param[in] T3           - a col-major 1 * (m - 1) matrix
 *  \param[in] T4           - a col-major (n - m) * (m - 1) matrix
 *  \param Q11              - an m1-batched row-major (n - m) * (n - m) matrix
 *  \param Q21              - an m1-batched row-major (n - m) * 1 matrix 
 *  \param Q31              - an m1-batched row-major (n - m) * (m - 1) matrix
 */
void calculate_F31(unsigned char *F31, const unsigned char *T3, const unsigned char *T4, const unsigned char *Q11, 
                    const unsigned char *Q21, const unsigned char *Q31);

/** \brief compute F52 from T1, Q12, Q22, Q52
 *  \param F52 
 *  \param[in] T1           - a col-major (n - m) * 1 matrix
 *  \param Q12              - an m1-batched row-major (n - m) * (n - m) matrix
 *  \param Q22              - an m1-batched row-major (n - m) * 1 matrix 
 *  \param Q52              - an m1-batched row-major 1 * 1 matrix 
 */
void calculate_F52(unsigned char *F52, const unsigned char *T1, const unsigned char *Q12, const unsigned char *Q22, 
                    const unsigned char *Q52);

/** \brief compute F62 from T1, T3, T4, Q12, Q22, Q32, Q52, Q62
 *  \param F62 
 *  \param[in] T1           - a col-major (n - m) * 1 matrix
 *  \param[in] T3           - a col-major 1 * (m - 1) matrix
 *  \param[in] T4           - a col-major (n - m) * (m - 1) matrix
 *  \param Q12              - an m1-batched row-major (n - m) * (n - m) matrix
 *  \param Q22              - an m1-batched row-major (n - m) * 1 matrix
 *  \param Q32              - an m1-batched row-major (n - m) * (m - 1) matrix
 *  \param Q52              - an m1-batched row-major 1 * 1 matrix 
 *  \param Q62              - an m1-batched row-major 1 * (m - 1) matrix
 */
void calculate_F62(unsigned char *F62, const unsigned char *T1, const unsigned char *T3, const unsigned char *T4, 
                    const unsigned char *Q12, const unsigned char *Q22, const unsigned char *Q32, const unsigned char *Q52, 
                    const unsigned char *Q62);

#endif
