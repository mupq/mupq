#ifndef _PARAMS_BISCUIT_H_
#define _PARAMS_BISCUIT_H_

/* #define COMPACT_SK */

#define CDIV8(x) (((x) + 7) / 8)

#define NBITS1(x) ((x) >= (1 << 0) ? 1 : 0)
#define NBITS2(x) ((x) >= (1 << 1) ? 1 + NBITS1 ((x) >> 1) : NBITS1 (x))
#define NBITS4(x) ((x) >= (1 << 2) ? 2 + NBITS2 ((x) >> 2) : NBITS2 (x))
#define NBITS8(x) ((x) >= (1 << 4) ? 4 + NBITS4 ((x) >> 4) : NBITS4 (x))
#define NBITS16(x) ((x) >= (1 << 8) ? 8 + NBITS8 ((x) >> 8) : NBITS8 (x))
#define NBITS32(x) ((x) >= (1 << 16) ? 16 + NBITS16 ((x) >> 16) : NBITS16 (x))
#define LOG2(x) NBITS32 ((x) - 1)

#define S_BITS(q, n) ((n) * LOG2 (q))
#define T_BITS(q, m) ((m) * LOG2 (q))
#define YZ_BITS(q, m, d) ((2 * (d) - 3) * (m) * LOG2 (q))
#define A_BITS(q, m, d) (((d) - 1) * (m) * LOG2 (q))

#ifndef COMPACT_SK
#define GET_PRIVKEY_BYTES(lambda, tau, N, q, n, m, d)                   \
  (((lambda) >> 3) +                                                    \
   CDIV8 (S_BITS((q), (n)) + T_BITS((q), (m)) + YZ_BITS((q), (m), (d))))
#else
#define GET_PRIVKEY_BYTES(lambda, tau, N, q, n, m, d) (2 * ((lambda) >> 3))
#endif
#define GET_PUBKEY_BYTES(lambda, tau, N, q, n, m, d)    \
  (((lambda) >> 3) + CDIV8 (T_BITS ((q), (m))))
#define GET_SIGNATURE_BYTES(lambda, tau, N, q, n, m, d)                 \
  (3 * ((lambda) >> 2)                                                  \
   + (tau) * (LOG2 (N) * ((lambda) >> 3) + ((lambda) >> 2))             \
   + CDIV8 ((tau) *                                                     \
            (S_BITS((q), (n)) + YZ_BITS((q), (m), (d))                  \
            + A_BITS((q), (m), (d)))))

#endif
