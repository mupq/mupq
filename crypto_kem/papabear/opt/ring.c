/** Optimized ring arithmetic */
#include "ring.h"

/** Karatsuba inner loop */
static inline void __attribute__((always_inline)) triplemac(
    dslimb_t *ac0, dslimb_t *ac1, dslimb_t *ac2,
    limb_t a0, limb_t a1, limb_t b0, limb_t b1
) {
    /* Multiply and accumulate three times into dslimbs, for Karatsuba */
    *ac0 += (dlimb_t)b0 * a1;
    *ac1 += (dlimb_t)b1 * (a1+a0);
    *ac2 += (dslimb_t)(slimb_t)(b0-b1) * (slimb_t)a0;
}

/** Multiply and accumulate c += a*b */
void mac(gf_t c, const gf_t a, const gf_t b) {
    /* One-level Karatsuba */
    /* PERF: for some reason whether this is signed or unsigned makes
     * a huge difference, depending on platform.
     */
    unsigned i,j;
    dslimb_t accum0 = 0, accum1 = 0;
    for (i=0; i<NLIMBS/2; i++) {
        dslimb_t accum2 = 0;

        for (j=0; j<=i; j++) {
            limb_t b0 = b[j], b1 = b[j+NLIMBS/2];
            limb_t a1 = a[i-j+NLIMBS/2], a0 = a[i-j];
            triplemac(&accum0,&accum1,&accum2,a0,a1,b0,b1);
        }

        accum0 = accum0 - accum2 - accum2;

        for (; j<NLIMBS/2; j++) {
            limb_t b0 = b[j], b1 = b[j+NLIMBS/2];
            limb_t a1 = a[i+NLIMBS-j], a0 = a[i+NLIMBS/2-j];
            triplemac(&accum1,&accum2,&accum0,a0,a1,b0,b1);
        }

        accum0 += accum2 + c[i];
        accum1 += accum2 + c[i+NLIMBS/2];
        c[i]          = accum0 & LMASK;
        c[i+NLIMBS/2] = accum1 & LMASK;

        accum0 >>= LBITS;
        accum1 >>= LBITS;
    }

    accum0 += accum1;
    accum0 += c[NLIMBS/2];
    accum1 += c[0];
    c[NLIMBS/2] = accum0 & LMASK;
    c[0] = accum1 & LMASK;

    accum0 >>= LBITS;
    accum1 >>= LBITS;

    c[NLIMBS/2+1] += accum0;
    c[1] += accum1;
}

typedef struct { limb_t x; } __attribute__((packed)) unaligned_limb_t;
static inline void write_limb(uint8_t out[sizeof(limb_t)], limb_t x) {
    ((unaligned_limb_t *)out)->x = x;
}

static inline limb_t read_limb(const uint8_t *in) {
    return ((const unaligned_limb_t *)in)->x;
}

/** Serialize a gf_t to bytes */
void contract(uint8_t ch[GF_BYTES], gf_t ll) {
    canon(ll);
    signed bbits=0;
    unsigned j=0;
    limb_t buffer=0;
    for (unsigned i=0;i<NLIMBS;i++) {
        if ((unsigned)bbits+LBITS >= 8*sizeof(limb_t)) {
            limb_t tmp = ll[i];
            buffer |= tmp<<bbits;
            write_limb(&ch[j], buffer);
            j += sizeof(limb_t);
            buffer = tmp>>(8*sizeof(limb_t)-bbits);
            bbits += LBITS - 8*sizeof(limb_t);
        } else {
            buffer |= ll[i]<<bbits;
            bbits += LBITS;
        }
    }
    for (; bbits>0; bbits-=8,buffer>>=8) {
        ch[j++] = buffer;
    }
}

/** Deserialize a gf_t from bytes */
void expand(gf_t ll, const uint8_t ch[GF_BYTES]) {
    limb_t buffer=0;
    for (unsigned i=0,j=0,bbits=0;;i+=sizeof(limb_t)) {
        if (i+sizeof(limb_t) >= GF_BYTES) {
            unsigned target = GF_BYTES-sizeof(limb_t);
            limb_t tmp = read_limb(&ch[target]) >> 8*(i-target);
            ll[j++] = (buffer | tmp<<bbits) & LMASK;
            return;
        } else {
            limb_t tmp = read_limb(&ch[i]);
            ll[j++] = (buffer | tmp<<bbits) & LMASK;
            buffer = tmp>>(LBITS-bbits);
            bbits += 8*sizeof(limb_t) - LBITS;
            while (bbits >= LBITS) {
                ll[j++] = buffer & LMASK;
                buffer >>= LBITS;
                bbits -= LBITS;
            }
        }
    }
}

/** Reduce a gf_t to canonical form, i.e. strictly less than N. */
void canon(gf_t c) {
    const limb_t DELTA = (limb_t)1<<(LBITS-1);

    /* Reduce to 0..2p */
    slimb_t hi = c[NLIMBS-1] - DELTA;
    c[NLIMBS-1] = (hi & LMASK) + DELTA;
    c[NLIMBS/2] += hi>>LBITS;

    /* Strong reduce.  First subtract modulus */
    dslimb_t scarry = hi>>LBITS;
    for (unsigned i=0; i<NLIMBS; i++) {
        scarry = scarry + (slimb_t)c[i] - modulus(i);
        c[i] = scarry & LMASK;
        scarry >>= LBITS;
    }

    /* add it back */
    dlimb_t carry = 0;
    for (unsigned i=0; i<NLIMBS; i++) {
        carry = carry + c[i] + (scarry & modulus(i));
        c[i] = carry & LMASK;
        carry >>= LBITS;
    }
    assert(carry+scarry==0);
}
