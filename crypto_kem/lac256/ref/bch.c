/*
 * This library provides encoding/decoding of binary
 * Bose-Chaudhuri-Hocquenghem (BCH) codes.
 *
 * Call init_bch to get a pointer to a newly allocated bch_control structure for
 * the given m (Galois field order), t (error correction capability) and
 * (optional) primitive polynomial parameters.
 *
 * Call encode_bch to compute and store ecc parity bytes to a given buffer.
 * Call decode_bch to detect and locate errors in received data.
 * Algorithmic details:
 *
 * Encoding is performed by processing 32 input bits in parallel, using 4
 * remainder lookup tables.
 *
 * The final stage of decoding involves the following internal steps:
 * a. Syndrome computation
 * b. Error locator polynomial computation using Berlekamp-Massey algorithm
 * c. Error locator root finding (by far the most expensive step)
 *

 */
#include "lac_param.h"
# include <stdint.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "bch.h"
# define EINVAL -75
# define EBADMSG -76
# define MODULE_LICENSE(s)
# define MODULE_AUTHOR(s)
# define MODULE_DESCRIPTION(s)
# define EXPORT_SYMBOL_GPL(x)
# define GFP_KERNEL
# define kmalloc(size, flags) malloc(size)
# define kzalloc(size, flags) memset(malloc(size), 0, size)
# define kfree(ptr) free(ptr)
# define DIV_ROUND_UP(a, b) ((a + b - 1) / b)
# define ARRAY_SIZE(arr) (sizeof(arr) / sizeof(*(arr)))
# define fls(x) ({ \
	unsigned int __tmp = x; \
	unsigned int __count = 0; \
	while (__tmp >>= 1) \
		__count++; \
	__count + 1; \
})

#define GF_M(_p)               ((_p)->m)
#define GF_T(_p)               ((_p)->t)
#define GF_N(_p)               ((_p)->n)

#define BCH_ECC_WORDS(_p)      DIV_ROUND_UP(GF_M(_p)*GF_T(_p), 32)
#define BCH_ECC_BYTES(_p)      DIV_ROUND_UP(GF_M(_p)*GF_T(_p), 8)

/*
 * represent a polynomial over GF(2^m)
 */
struct gf_poly {
	unsigned int deg;    /* polynomial degree */
	unsigned int c[0];   /* polynomial terms */
};

/* given its degree, compute a polynomial size in bytes */
#define GF_POLY_SZ(_d) (sizeof(struct gf_poly)+((_d)+1)*sizeof(unsigned int))

/* polynomial of degree 1 */
struct gf_poly_deg1 {
	struct gf_poly poly;
	unsigned int   c[2];
};

static int32_t cpu_to_be32(int32_t x)
{
	return (x<<24)|((x<<8)&0xff0000)|((x>>8)&0xff00)|((x>>24)&0xff);
}

/*
 * same as encode_bch(), but process input data one byte at a time
 */
static void encode_bch_unaligned(struct bch_control *bch,
				 const unsigned char *data, unsigned int len,
				 uint32_t *ecc)
{
	int i;
	const uint32_t *p;
	const int l = BCH_ECC_WORDS(bch)-1;

	while (len--) {
		p = bch->mod8_tab + (l+1)*(((ecc[0] >> 24)^(*data++)) & 0xff);

		for (i = 0; i < l; i++)
			ecc[i] = ((ecc[i] << 8)|(ecc[i+1] >> 24))^(*p++);

		ecc[l] = (ecc[l] << 8)^(*p);
	}
}

/*
 * convert ecc bytes to aligned, zero-padded 32-bit ecc words
 */
static void load_ecc8(struct bch_control *bch, uint32_t *dst,
		      const uint8_t *src)
{
	uint8_t pad[4] = {0, 0, 0, 0};
	unsigned int i, nwords = BCH_ECC_WORDS(bch)-1;

	for (i = 0; i < nwords; i++, src += 4)
		dst[i] = (src[0] << 24)|(src[1] << 16)|(src[2] << 8)|src[3];

	memcpy(pad, src, BCH_ECC_BYTES(bch)-4*nwords);
	dst[nwords] = (pad[0] << 24)|(pad[1] << 16)|(pad[2] << 8)|pad[3];
}

/*
 * convert 32-bit ecc words to ecc bytes
 */
static void store_ecc8(struct bch_control *bch, uint8_t *dst,
		       const uint32_t *src)
{
	uint8_t pad[4];
	unsigned int i, nwords = BCH_ECC_WORDS(bch)-1;

	for (i = 0; i < nwords; i++) {
		*dst++ = (src[i] >> 24);
		*dst++ = (src[i] >> 16) & 0xff;
		*dst++ = (src[i] >>  8) & 0xff;
		*dst++ = (src[i] >>  0) & 0xff;
	}
	pad[0] = (src[nwords] >> 24);
	pad[1] = (src[nwords] >> 16) & 0xff;
	pad[2] = (src[nwords] >>  8) & 0xff;
	pad[3] = (src[nwords] >>  0) & 0xff;
	memcpy(dst, pad, BCH_ECC_BYTES(bch)-4*nwords);
}

/**
 * encode_bch - calculate BCH ecc parity of data
 * @bch:   BCH control structure
 * @data:  data to encode
 * @len:   data length in bytes
 * @ecc:   ecc parity data, must be initialized by caller
 *
 * The @ecc parity array is used both as input and output parameter, in order to
 * allow incremental computations. It should be of the size indicated by member
 * @ecc_bytes of @bch, and should be initialized to 0 before the first call.
 *
 * The exact number of computed ecc parity bits is given by member @ecc_bits of
 * @bch; it may be less than m*t for large values of t.
 */
void encode_bch(struct bch_control *bch, const uint8_t *data,
		unsigned int len, uint8_t *ecc)
{
	const unsigned int l = BCH_ECC_WORDS(bch)-1;
	unsigned int i, mlen;
	unsigned long m;
	uint32_t w, r[l+1];
	const uint32_t * const tab0 = bch->mod8_tab;
	const uint32_t * const tab1 = tab0 + 256*(l+1);
	const uint32_t * const tab2 = tab1 + 256*(l+1);
	const uint32_t * const tab3 = tab2 + 256*(l+1);
	const uint32_t *pdata, *p0, *p1, *p2, *p3;

	if (ecc) {
		/* load ecc parity bytes into internal 32-bit buffer */
		load_ecc8(bch, bch->ecc_buf, ecc);
	} else {
		memset(bch->ecc_buf, 0, sizeof(r));
	}

	/* process first unaligned data bytes */
	m = ((unsigned long)data) & 3;
	if (m) {
		mlen = (len < (4-m)) ? len : 4-m;
		encode_bch_unaligned(bch, data, mlen, bch->ecc_buf);
		data += mlen;
		len  -= mlen;
	}

	/* process 32-bit aligned data words */
	pdata = (uint32_t *)data;
	mlen  = len/4;
	data += 4*mlen;
	len  -= 4*mlen;
	memcpy(r, bch->ecc_buf, sizeof(r));

	/*
	 * split each 32-bit word into 4 polynomials of weight 8 as follows:
	 *
	 * 31 ...24  23 ...16  15 ... 8  7 ... 0
	 * xxxxxxxx  yyyyyyyy  zzzzzzzz  tttttttt
	 *                               tttttttt  mod g = r0 (precomputed)
	 *                     zzzzzzzz  00000000  mod g = r1 (precomputed)
	 *           yyyyyyyy  00000000  00000000  mod g = r2 (precomputed)
	 * xxxxxxxx  00000000  00000000  00000000  mod g = r3 (precomputed)
	 * xxxxxxxx  yyyyyyyy  zzzzzzzz  tttttttt  mod g = r0^r1^r2^r3
	 */
	while (mlen--) {
		/* input data is read in big-endian format */
		w = r[0]^cpu_to_be32(*pdata++);
		p0 = tab0 + (l+1)*((w >>  0) & 0xff);
		p1 = tab1 + (l+1)*((w >>  8) & 0xff);
		p2 = tab2 + (l+1)*((w >> 16) & 0xff);
		p3 = tab3 + (l+1)*((w >> 24) & 0xff);

		for (i = 0; i < l; i++)
			r[i] = r[i+1]^p0[i]^p1[i]^p2[i]^p3[i];

		r[l] = p0[l]^p1[l]^p2[l]^p3[l];
	}
	memcpy(bch->ecc_buf, r, sizeof(r));

	/* process last unaligned bytes */
	if (len)
		encode_bch_unaligned(bch, data, len, bch->ecc_buf);

	/* store ecc parity bytes into original parity buffer */
	if (ecc)
		store_ecc8(bch, ecc, bch->ecc_buf);
}

static inline int modulo(struct bch_control *bch, unsigned int v)
{
	const unsigned int n = GF_N(bch);
	while (v >= n) {
		v -= n;
		v = (v & n) + (v >> GF_M(bch));
	}
	return v;
}

/*
 * shorter and faster modulo function, only works when v < 2N.
 */
static inline int mod_s(struct bch_control *bch, unsigned int v)
{
	const unsigned int n = GF_N(bch);
	return (v < n) ? v : v-n;
}

static inline int deg(unsigned int poly)
{
	/* polynomial degree is the most-significant bit index */
	return fls(poly)-1;
}

/* Galois field basic operations: multiply, divide, inverse, etc. */

static inline unsigned int gf_mul(struct bch_control *bch, unsigned int a,
				  unsigned int b)
{
	return (a && b) ? bch->a_pow_tab[mod_s(bch, bch->a_log_tab[a]+
					       bch->a_log_tab[b])] : 0;
}

static inline unsigned int gf_sqr(struct bch_control *bch, unsigned int a)
{
	return a ? bch->a_pow_tab[mod_s(bch, 2*bch->a_log_tab[a])] : 0;
}

static inline unsigned int gf_div(struct bch_control *bch, unsigned int a,
				  unsigned int b)
{
	return a ? bch->a_pow_tab[mod_s(bch, bch->a_log_tab[a]+
					GF_N(bch)-bch->a_log_tab[b])] : 0;
}

static inline unsigned int gf_inv(struct bch_control *bch, unsigned int a)
{
	return bch->a_pow_tab[GF_N(bch)-bch->a_log_tab[a]];
}

static inline unsigned int a_pow(struct bch_control *bch, int i)
{
	return bch->a_pow_tab[modulo(bch, i)];
}

static inline int a_log(struct bch_control *bch, unsigned int x)
{
	return bch->a_log_tab[x];
}

static inline int a_ilog(struct bch_control *bch, unsigned int x)
{
	return mod_s(bch, GF_N(bch)-bch->a_log_tab[x]);
}
static inline int parity(unsigned int x)
{
	/*
	 * public domain code snippet, lifted from
	 * http://www-graphics.stanford.edu/~seander/bithacks.html
	 */
	x ^= x >> 1;
	x ^= x >> 2;
	x = (x & 0x11111111U) * 0x11111111U;
	return (x >> 28) & 1;
}



#ifdef BCH_CONSTANT_TIME

/*
 * compute 2t syndromes of ecc polynomial, i.e. ecc(a^j) for j=1..2t
 */
static void compute_syndromes(struct bch_control *bch, uint32_t *ecc,
			      unsigned int *syn)
{
	int i, j, s;
	unsigned int m;
	uint32_t poly;
	const int t = GF_T(bch);

	s = bch->ecc_bits;

	/* make sure extra bits in last ecc word are cleared */
	m = ((unsigned int)s) & 31;
	if (m)
		ecc[s/32] &= ~((1u << (32-m))-1);
	memset(syn, 0, 2*t*sizeof(*syn));

	/* compute v(a^j) for j=1 .. 2t-1 */
	do {
		poly = *ecc++;
		s -= 32;

		for(i=31;i>=0;i--)
		{
			for (j = 0; j < 2*t; j += 2)
			{
				syn[j] ^= (a_pow(bch, (j+1)*(i+s))&(-((poly>>i)&0x1)));
			//	printf("index:%d\n",(j+1)*(i+s));
			}
		}
	} while (s > 0);

	/* v(a^(2j)) = v(a^j)^2 */
	for (j = 0; j < t; j++)
		syn[2*j+1] = gf_sqr(bch, syn[j]);
}

static void gf_poly_copy(struct gf_poly *dst, struct gf_poly *src)
{
	memcpy(dst, src, GF_POLY_SZ(src->deg));
}


static int compute_error_locator_polynomial(struct bch_control *bch,
					    const unsigned int *syn)
{
	const unsigned int t = GF_T(bch);
	const unsigned int n = GF_N(bch);
	unsigned int i, j, tmp, l, pd = 1, d = syn[0];
	struct gf_poly *elp = bch->elp;
	struct gf_poly *pelp = bch->poly_2t[0];
	struct gf_poly *elp_copy = bch->poly_2t[1];
	int k, pp = -1;
	uint16_t mask_d,mask_pelp;

	memset(pelp, 0, GF_POLY_SZ(2*t));
	memset(elp, 0, GF_POLY_SZ(2*t));

	pelp->deg = 0;
	pelp->c[0] = 1;
	elp->deg = 0;
	elp->c[0] = 1;

	/* use simplified binary Berlekamp-Massey algorithm */
	for (i = 0; (i < t) && (elp->deg <= t); i++) {
		//if (d) {
			mask_d=(d?0xffff:0x0);
			k = 2*i-pp;

			gf_poly_copy(elp_copy, elp);
			/* e[i+1](X) = e[i](X)+di*dp^-1*X^2(i-p)*e[p](X) */
			tmp = a_log(bch, d)+n-a_log(bch, pd);
			for (j = 0; j <= pelp->deg; j++)
			{
				mask_pelp=(pelp->c[j]? 0xffff:0x0);
				l = a_log(bch, pelp->c[j]);
				elp->c[j+k] ^= (a_pow(bch, tmp+l)&mask_d&mask_pelp);
			}
			/* compute l[i+1] = max(l[i]->c[l[p]+2*(i-p]) */
			tmp = pelp->deg+(k&mask_d);
			if (tmp > elp->deg) {
				elp->deg = tmp;
				gf_poly_copy(pelp, elp_copy);
				pd = d;
				pp = 2*i;
			}
			else //just for constant time
			{
				tmp=elp->deg;
				gf_poly_copy(elp_copy, pelp);
				tmp=pd;
				tmp=2*i;
			}
	//	}
		/* di+1 = S(2i+3)+elp[i+1].1*S(2i+2)+...+elp[i+1].lS(2i+3-l) */
		if (i < t-1) {
			d = syn[2*i+2];
			for (j = 1; j <= elp->deg; j++)
				d ^= gf_mul(bch, elp->c[j], syn[2*i+2-j]);
		}
	}

	return (elp->deg > t) ? -1 : (int)elp->deg;
}

/*
 * build monic, log-based representation of a polynomial
 */
static void init_rep(struct bch_control *bch,
			   const struct gf_poly *a, uint16_t *rep, uint16_t *syn_mask, uint16_t *syn_a, unsigned int pow_start)
{
	unsigned int i, d = a->deg;
    int l = GF_N(bch)-a_log(bch, a->c[a->deg]);

	for (i = 0; i < d; i++)
	{
		rep[i] =  modulo(bch, a_log(bch, a->c[i])+l + pow_start*i);
	    syn_mask[i]= a->c[i] ? 0xFFFF : 0x0;
		syn_a[i]=i;
	}
	rep[d]=modulo(bch, pow_start*d);
	syn_mask[d]=0xFFFF;
	syn_a[d]=d;
	for(i=d+1;i<=bch->t;i++)
	{
		rep[i]=0;
		syn_mask[i]=0;
		syn_a[i]=i;
	}
}

/*
 * exhaustive root search (Chien) implementation - not used, included only for
 * reference/comparison tests
 */
static int chien_search(struct bch_control *bch, unsigned int len,
			struct gf_poly *p, unsigned int *roots)
{
	unsigned int i, j, syn, syn0, count = 0;
	const unsigned int k = 8*len+bch->ecc_bits,bound=GF_N(bch)-bch->ecc_bits;
	uint16_t syn_mask[bch->t+1],syn_rep[bch->t+1],syn_a[bch->t+1],syn_tmp;

	/* use a log-based representation of polynomial */
	init_rep(bch,p,syn_rep,syn_mask,syn_a,GF_N(bch)-k);

	syn0 = gf_div(bch, p->c[0], p->c[p->deg]);

	for (i = GF_N(bch)-k+1; i <= bound; i++) {
		/* compute elp(a^i) */

		for (j = 1, syn = syn0; j <= bch->t; j++) {
			syn_rep[j]+=syn_a[j];
			syn_tmp=syn_rep[j];
		    syn_rep[j]=(syn_tmp<bch->n ? syn_tmp: syn_tmp-bch->n);
			syn ^= (bch->a_pow_tab[syn_rep[j]]&syn_mask[j]);

		}
		roots[count] = GF_N(bch)-i;
		count+=(syn==0);
	}

	return count;
}

/**
 * decode_bch - decode received codeword and find bit error locations
 * @bch:      BCH control structure
 * @data:     received data, ignored if @calc_ecc is provided
 * @len:      data length in bytes, must always be provided
 * @recv_ecc: received ecc, if NULL then assume it was XORed in @calc_ecc
 *
 * Returns:
 *  The number of errors found, or -EBADMSG if decoding failed, or -EINVAL if
 *  invalid parameters were provided
 */
int decode_bch(struct bch_control *bch, const uint8_t *data, unsigned int len,
	       const uint8_t *recv_ecc, const uint8_t *calc_ecc,
	       const unsigned int *syn, unsigned int *errloc)
{
    (void)syn;
    (void)calc_ecc;
	const unsigned int ecc_words = BCH_ECC_WORDS(bch);
	unsigned int nbits;
	int i, err;

	/* sanity check: make sure data length can be handled */
	if (8*len > (bch->n-bch->ecc_bits))
		return -EINVAL;
	//check data and ecc pointer
    if (!data || !recv_ecc)
		return -EINVAL;

	/* compute received data ecc into an internal buffer */
	encode_bch(bch, data, len, NULL);
	/* load received ecc  */
	load_ecc8(bch, bch->ecc_buf2, recv_ecc);
	/* XOR received and calculated ecc */
	for (i = 0; i < (int)ecc_words; i++)
	{
		bch->ecc_buf[i] ^= bch->ecc_buf2[i];
	}

	//compute syndromes
	compute_syndromes(bch, bch->ecc_buf, bch->syn);
    //compute error locator polynomial
	compute_error_locator_polynomial(bch, bch->syn);
	//find roots
	err=chien_search(bch, len, bch->elp, errloc);

	//post process error location
	nbits = (len*8)+bch->ecc_bits;

	for (i = 0; i < err; i++)
	{
		errloc[i] = nbits-1-errloc[i];
		errloc[i] = (errloc[i] & ~7)|(7-(errloc[i] & 7));
	}

	return err;
}

#else

/*
 * compute 2t syndromes of ecc polynomial, i.e. ecc(a^j) for j=1..2t
 */
static void compute_syndromes(struct bch_control *bch, uint32_t *ecc,
			      unsigned int *syn)
{
	int i, j, s;
	unsigned int m;
	uint32_t poly;
	const int t = GF_T(bch);

	s = bch->ecc_bits;

	/* make sure extra bits in last ecc word are cleared */
	m = ((unsigned int)s) & 31;
	if (m)
		ecc[s/32] &= ~((1u << (32-m))-1);
	memset(syn, 0, 2*t*sizeof(*syn));

	/* compute v(a^j) for j=1 .. 2t-1 */
	do {
		poly = *ecc++;
		s -= 32;
		while (poly) {
			i = deg(poly);
			for (j = 0; j < 2*t; j += 2)
				syn[j] ^= a_pow(bch, (j+1)*(i+s));

			poly ^= (1 << i);
		}
	} while (s > 0);

	/* v(a^(2j)) = v(a^j)^2 */
	for (j = 0; j < t; j++)
		syn[2*j+1] = gf_sqr(bch, syn[j]);
}

static void gf_poly_copy(struct gf_poly *dst, struct gf_poly *src)
{
	memcpy(dst, src, GF_POLY_SZ(src->deg));
}

static int compute_error_locator_polynomial(struct bch_control *bch,
					    const unsigned int *syn)
{
	const unsigned int t = GF_T(bch);
	const unsigned int n = GF_N(bch);
	unsigned int i, j, tmp, l, pd = 1, d = syn[0];
	struct gf_poly *elp = bch->elp;
	struct gf_poly *pelp = bch->poly_2t[0];
	struct gf_poly *elp_copy = bch->poly_2t[1];
	int k, pp = -1;

	memset(pelp, 0, GF_POLY_SZ(2*t));
	memset(elp, 0, GF_POLY_SZ(2*t));

	pelp->deg = 0;
	pelp->c[0] = 1;
	elp->deg = 0;
	elp->c[0] = 1;

	/* use simplified binary Berlekamp-Massey algorithm */
	for (i = 0; (i < t) && (elp->deg <= t); i++) {
		if (d) {
			k = 2*i-pp;
			gf_poly_copy(elp_copy, elp);
			/* e[i+1](X) = e[i](X)+di*dp^-1*X^2(i-p)*e[p](X) */
			tmp = a_log(bch, d)+n-a_log(bch, pd);
			for (j = 0; j <= pelp->deg; j++) {
				if (pelp->c[j]) {
					l = a_log(bch, pelp->c[j]);
					elp->c[j+k] ^= a_pow(bch, tmp+l);
				}
			}
			/* compute l[i+1] = max(l[i]->c[l[p]+2*(i-p]) */
			tmp = pelp->deg+k;
			if (tmp > elp->deg) {
				elp->deg = tmp;
				gf_poly_copy(pelp, elp_copy);
				pd = d;
				pp = 2*i;
			}
		}
		/* di+1 = S(2i+3)+elp[i+1].1*S(2i+2)+...+elp[i+1].lS(2i+3-l) */
		if (i < t-1) {
			d = syn[2*i+2];
			for (j = 1; j <= elp->deg; j++)
				d ^= gf_mul(bch, elp->c[j], syn[2*i+2-j]);
		}
	}
//	dbg("elp=%s\n", gf_poly_str(elp));
	return (elp->deg > t) ? -1 : (int)elp->deg;
}

/*
 * solve a m x m linear system in GF(2) with an expected number of solutions,
 * and return the number of found solutions
 */
static int solve_linear_system(struct bch_control *bch, unsigned int *rows,
			       unsigned int *sol, int nsol)
{
	const int m = GF_M(bch);
	unsigned int tmp, mask;
	int rem, c, r, p, k, param[m];

	k = 0;
	mask = 1 << m;

	/* Gaussian elimination */
	for (c = 0; c < m; c++) {
		rem = 0;
		p = c-k;
		/* find suitable row for elimination */
		for (r = p; r < m; r++) {
			if (rows[r] & mask) {
				if (r != p) {
					tmp = rows[r];
					rows[r] = rows[p];
					rows[p] = tmp;
				}
				rem = r+1;
				break;
			}
		}
		if (rem) {
			/* perform elimination on remaining rows */
			tmp = rows[p];
			for (r = rem; r < m; r++) {
				if (rows[r] & mask)
					rows[r] ^= tmp;
			}
		} else {
			/* elimination not needed, store defective row index */
			param[k++] = c;
		}
		mask >>= 1;
	}
	/* rewrite system, inserting fake parameter rows */
	if (k > 0) {
		p = k;
		for (r = m-1; r >= 0; r--) {
			if ((r > m-1-k) && rows[r])
				/* system has no solution */
				return 0;

			rows[r] = (p && (r == param[p-1])) ?
				p--, 1u << (m-r) : rows[r-p];
		}
	}

	if (nsol != (1 << k))
		/* unexpected number of solutions */
		return 0;

	for (p = 0; p < nsol; p++) {
		/* set parameters for p-th solution */
		for (c = 0; c < k; c++)
			rows[param[c]] = (rows[param[c]] & ~1)|((p >> c) & 1);

		/* compute unique solution */
		tmp = 0;
		for (r = m-1; r >= 0; r--) {
			mask = rows[r] & (tmp|1);
			tmp |= parity(mask) << (m-r);
		}
		sol[p] = tmp >> 1;
	}
	return nsol;
}

/*
 * this function builds and solves a linear system for finding roots of a degree
 * 4 affine monic polynomial X^4+aX^2+bX+c over GF(2^m).
 */
static int find_affine4_roots(struct bch_control *bch, unsigned int a,
			      unsigned int b, unsigned int c,
			      unsigned int *roots)
{
	int i, j, k;
	const int m = GF_M(bch);
	unsigned int mask = 0xff, t, rows[16] = {0,};

	j = a_log(bch, b);
	k = a_log(bch, a);
	rows[0] = c;

	/* buid linear system to solve X^4+aX^2+bX+c = 0 */
	for (i = 0; i < m; i++) {
		rows[i+1] = bch->a_pow_tab[4*i]^
			(a ? bch->a_pow_tab[mod_s(bch, k)] : 0)^
			(b ? bch->a_pow_tab[mod_s(bch, j)] : 0);
		j++;
		k += 2;
	}
	/*
	 * transpose 16x16 matrix before passing it to linear solver
	 * warning: this code assumes m < 16
	 */
	for (j = 8; j != 0; j >>= 1, mask ^= (mask << j)) {
		for (k = 0; k < 16; k = (k+j+1) & ~j) {
			t = ((rows[k] >> j)^rows[k+j]) & mask;
			rows[k] ^= (t << j);
			rows[k+j] ^= t;
		}
	}
	return solve_linear_system(bch, rows, roots, 4);
}

/*
 * compute root r of a degree 1 polynomial over GF(2^m) (returned as log(1/r))
 */
static int find_poly_deg1_roots(struct bch_control *bch, struct gf_poly *poly,
				unsigned int *roots)
{
	int n = 0;

	if (poly->c[0])
		/* poly[X] = bX+c with c!=0, root=c/b */
		roots[n++] = mod_s(bch, GF_N(bch)-bch->a_log_tab[poly->c[0]]+
				   bch->a_log_tab[poly->c[1]]);
	return n;
}

/*
 * compute roots of a degree 2 polynomial over GF(2^m)
 */
static int find_poly_deg2_roots(struct bch_control *bch, struct gf_poly *poly,
				unsigned int *roots)
{
	int n = 0, i, l0, l1, l2;
	unsigned int u, v, r;

	if (poly->c[0] && poly->c[1]) {

		l0 = bch->a_log_tab[poly->c[0]];
		l1 = bch->a_log_tab[poly->c[1]];
		l2 = bch->a_log_tab[poly->c[2]];

		/* using z=a/bX, transform aX^2+bX+c into z^2+z+u (u=ac/b^2) */
		u = a_pow(bch, l0+l2+2*(GF_N(bch)-l1));
		/*
		 * let u = sum(li.a^i) i=0..m-1; then compute r = sum(li.xi):
		 * r^2+r = sum(li.(xi^2+xi)) = sum(li.(a^i+Tr(a^i).a^k)) =
		 * u + sum(li.Tr(a^i).a^k) = u+a^k.Tr(sum(li.a^i)) = u+a^k.Tr(u)
		 * i.e. r and r+1 are roots iff Tr(u)=0
		 */
		r = 0;
		v = u;
		while (v) {
			i = deg(v);
			r ^= bch->xi_tab[i];
			v ^= (1 << i);
		}
		/* verify root */
		if ((gf_sqr(bch, r)^r) == u) {
			/* reverse z=a/bX transformation and compute log(1/r) */
			roots[n++] = modulo(bch, 2*GF_N(bch)-l1-
					    bch->a_log_tab[r]+l2);
			roots[n++] = modulo(bch, 2*GF_N(bch)-l1-
					    bch->a_log_tab[r^1]+l2);
		}
	}
	return n;
}

/*
 * compute roots of a degree 3 polynomial over GF(2^m)
 */
static int find_poly_deg3_roots(struct bch_control *bch, struct gf_poly *poly,
				unsigned int *roots)
{
	int i, n = 0;
	unsigned int a, b, c, a2, b2, c2, e3, tmp[4];

	if (poly->c[0]) {
		/* transform polynomial into monic X^3 + a2X^2 + b2X + c2 */
		e3 = poly->c[3];
		c2 = gf_div(bch, poly->c[0], e3);
		b2 = gf_div(bch, poly->c[1], e3);
		a2 = gf_div(bch, poly->c[2], e3);

		/* (X+a2)(X^3+a2X^2+b2X+c2) = X^4+aX^2+bX+c (affine) */
		c = gf_mul(bch, a2, c2);           /* c = a2c2      */
		b = gf_mul(bch, a2, b2)^c2;        /* b = a2b2 + c2 */
		a = gf_sqr(bch, a2)^b2;            /* a = a2^2 + b2 */

		/* find the 4 roots of this affine polynomial */
		if (find_affine4_roots(bch, a, b, c, tmp) == 4) {
			/* remove a2 from final list of roots */
			for (i = 0; i < 4; i++) {
				if (tmp[i] != a2)
					roots[n++] = a_ilog(bch, tmp[i]);
			}
		}
	}
	return n;
}

/*
 * compute roots of a degree 4 polynomial over GF(2^m)
 */
static int find_poly_deg4_roots(struct bch_control *bch, struct gf_poly *poly,
				unsigned int *roots)
{
	int i, l, n = 0;
	unsigned int a, b, c, d, e = 0, f, a2, b2, c2, e4;

	if (poly->c[0] == 0)
		return 0;

	/* transform polynomial into monic X^4 + aX^3 + bX^2 + cX + d */
	e4 = poly->c[4];
	d = gf_div(bch, poly->c[0], e4);
	c = gf_div(bch, poly->c[1], e4);
	b = gf_div(bch, poly->c[2], e4);
	a = gf_div(bch, poly->c[3], e4);

	/* use Y=1/X transformation to get an affine polynomial */
	if (a) {
		/* first, eliminate cX by using z=X+e with ae^2+c=0 */
		if (c) {
			/* compute e such that e^2 = c/a */
			f = gf_div(bch, c, a);
			l = a_log(bch, f);
			l += (l & 1) ? GF_N(bch) : 0;
			e = a_pow(bch, l/2);
			/*
			 * use transformation z=X+e:
			 * z^4+e^4 + a(z^3+ez^2+e^2z+e^3) + b(z^2+e^2) +cz+ce+d
			 * z^4 + az^3 + (ae+b)z^2 + (ae^2+c)z+e^4+be^2+ae^3+ce+d
			 * z^4 + az^3 + (ae+b)z^2 + e^4+be^2+d
			 * z^4 + az^3 +     b'z^2 + d'
			 */
			d = a_pow(bch, 2*l)^gf_mul(bch, b, f)^d;
			b = gf_mul(bch, a, e)^b;
		}
		/* now, use Y=1/X to get Y^4 + b/dY^2 + a/dY + 1/d */
		if (d == 0)
			/* assume all roots have multiplicity 1 */
			return 0;

		c2 = gf_inv(bch, d);
		b2 = gf_div(bch, a, d);
		a2 = gf_div(bch, b, d);
	} else {
		/* polynomial is already affine */
		c2 = d;
		b2 = c;
		a2 = b;
	}
	/* find the 4 roots of this affine polynomial */
	if (find_affine4_roots(bch, a2, b2, c2, roots) == 4) {
		for (i = 0; i < 4; i++) {
			/* post-process roots (reverse transformations) */
			f = a ? gf_inv(bch, roots[i]) : roots[i];
			roots[i] = a_ilog(bch, f^e);
		}
		n = 4;
	}
	return n;
}

/*
 * build monic, log-based representation of a polynomial
 */
static void gf_poly_logrep(struct bch_control *bch,
			   const struct gf_poly *a, int *rep)
{
	int i, d = a->deg, l = GF_N(bch)-a_log(bch, a->c[a->deg]);

	/* represent 0 values with -1; warning, rep[d] is not set to 1 */
	for (i = 0; i < d; i++)
		rep[i] = a->c[i] ? mod_s(bch, a_log(bch, a->c[i])+l) : -1;
}

/*
 * compute polynomial Euclidean division remainder in GF(2^m)[X]
 */
static void gf_poly_mod(struct bch_control *bch, struct gf_poly *a,
			const struct gf_poly *b, int *rep)
{
	int la, p, m;
	unsigned int i, j, *c = a->c;
	const unsigned int d = b->deg;

	if (a->deg < d)
		return;

	/* reuse or compute log representation of denominator */
	if (!rep) {
		rep = bch->cache;
		gf_poly_logrep(bch, b, rep);
	}

	for (j = a->deg; j >= d; j--) {
		if (c[j]) {
			la = a_log(bch, c[j]);
			p = j-d;
			for (i = 0; i < d; i++, p++) {
				m = rep[i];
				if (m >= 0)
					c[p] ^= bch->a_pow_tab[mod_s(bch,
								     m+la)];
			}
		}
	}
	a->deg = d-1;
	while (!c[a->deg] && a->deg)
		a->deg--;
}

/*
 * compute polynomial Euclidean division quotient in GF(2^m)[X]
 */
static void gf_poly_div(struct bch_control *bch, struct gf_poly *a,
			const struct gf_poly *b, struct gf_poly *q)
{
	if (a->deg >= b->deg) {
		q->deg = a->deg-b->deg;
		/* compute a mod b (modifies a) */
		gf_poly_mod(bch, a, b, NULL);
		/* quotient is stored in upper part of polynomial a */
		memcpy(q->c, &a->c[b->deg], (1+q->deg)*sizeof(unsigned int));
	} else {
		q->deg = 0;
		q->c[0] = 0;
	}
}

/*
 * compute polynomial GCD (Greatest Common Divisor) in GF(2^m)[X]
 */
static struct gf_poly *gf_poly_gcd(struct bch_control *bch, struct gf_poly *a,
				   struct gf_poly *b)
{
	struct gf_poly *tmp;

//	dbg("gcd(%s,%s)=", gf_poly_str(a), gf_poly_str(b));

	if (a->deg < b->deg) {
		tmp = b;
		b = a;
		a = tmp;
	}

	while (b->deg > 0) {
		gf_poly_mod(bch, a, b, NULL);
		tmp = b;
		b = a;
		a = tmp;
	}

//	dbg("%s\n", gf_poly_str(a));

	return a;
}

/*
 * Given a polynomial f and an integer k, compute Tr(a^kX) mod f
 * This is used in Berlekamp Trace algorithm for splitting polynomials
 */
static void compute_trace_bk_mod(struct bch_control *bch, int k,
				 const struct gf_poly *f, struct gf_poly *z,
				 struct gf_poly *out)
{
	const int m = GF_M(bch);
	int i, j;

	/* z contains z^2j mod f */
	z->deg = 1;
	z->c[0] = 0;
	z->c[1] = bch->a_pow_tab[k];

	out->deg = 0;
	memset(out, 0, GF_POLY_SZ(f->deg));

	/* compute f log representation only once */
	gf_poly_logrep(bch, f, bch->cache);

	for (i = 0; i < m; i++) {
		/* add a^(k*2^i)(z^(2^i) mod f) and compute (z^(2^i) mod f)^2 */
		for (j = z->deg; j >= 0; j--) {
			out->c[j] ^= z->c[j];
			z->c[2*j] = gf_sqr(bch, z->c[j]);
			z->c[2*j+1] = 0;
		}
		if (z->deg > out->deg)
			out->deg = z->deg;

		if (i < m-1) {
			z->deg *= 2;
			/* z^(2(i+1)) mod f = (z^(2^i) mod f)^2 mod f */
			gf_poly_mod(bch, z, f, bch->cache);
		}
	}
	while (!out->c[out->deg] && out->deg)
		out->deg--;

//	dbg("Tr(a^%d.X) mod f = %s\n", k, gf_poly_str(out));
}

/*
 * factor a polynomial using Berlekamp Trace algorithm (BTA)
 */
static void factor_polynomial(struct bch_control *bch, int k, struct gf_poly *f,
			      struct gf_poly **g, struct gf_poly **h)
{
	struct gf_poly *f2 = bch->poly_2t[0];
	struct gf_poly *q  = bch->poly_2t[1];
	struct gf_poly *tk = bch->poly_2t[2];
	struct gf_poly *z  = bch->poly_2t[3];
	struct gf_poly *gcd;

//	dbg("factoring %s...\n", gf_poly_str(f));

	*g = f;
	*h = NULL;

	/* tk = Tr(a^k.X) mod f */
	compute_trace_bk_mod(bch, k, f, z, tk);

	if (tk->deg > 0) {
		/* compute g = gcd(f, tk) (destructive operation) */
		gf_poly_copy(f2, f);
		gcd = gf_poly_gcd(bch, f2, tk);
		if (gcd->deg < f->deg) {
			/* compute h=f/gcd(f,tk); this will modify f and q */
			gf_poly_div(bch, f, gcd, q);
			/* store g and h in-place (clobbering f) */
			*h = &((struct gf_poly_deg1 *)f)[gcd->deg].poly;
			gf_poly_copy(*g, gcd);
			gf_poly_copy(*h, q);
		}
	}
}

/*
 * find roots of a polynomial, using BTZ algorithm; see the beginning of this
 * file for details
 */
static int find_poly_roots(struct bch_control *bch, unsigned int k,
			   struct gf_poly *poly, unsigned int *roots)
{
	int cnt;
	struct gf_poly *f1, *f2;

	switch (poly->deg) {
		/* handle low degree polynomials with ad hoc techniques */
	case 1:
		cnt = find_poly_deg1_roots(bch, poly, roots);
		break;
	case 2:
		cnt = find_poly_deg2_roots(bch, poly, roots);
		break;
	case 3:
		cnt = find_poly_deg3_roots(bch, poly, roots);
		break;
	case 4:
		cnt = find_poly_deg4_roots(bch, poly, roots);
		break;
	default:
		/* factor polynomial using Berlekamp Trace Algorithm (BTA) */
		cnt = 0;
		if (poly->deg && (k <= GF_M(bch))) {
			factor_polynomial(bch, k, poly, &f1, &f2);
			if (f1)
				cnt += find_poly_roots(bch, k+1, f1, roots);
			if (f2)
				cnt += find_poly_roots(bch, k+1, f2, roots+cnt);
		}
		break;
	}
	return cnt;
}


/**
 * decode_bch - decode received codeword and find bit error locations
 * @bch:      BCH control structure
 * @data:     received data, ignored if @calc_ecc is provided
 * @len:      data length in bytes, must always be provided
 * @recv_ecc: received ecc, if NULL then assume it was XORed in @calc_ecc
 * @calc_ecc: calculated ecc, if NULL then calc_ecc is computed from @data
 * @syn:      hw computed syndrome data (if NULL, syndrome is calculated)
 * @errloc:   output array of error locations
 *
 * Returns:
 *  The number of errors found, or -EBADMSG if decoding failed, or -EINVAL if
 *  invalid parameters were provided
 *
 * Depending on the available hw BCH support and the need to compute @calc_ecc
 * separately (using encode_bch()), this function should be called with one of
 * the following parameter configurations -
 *
 * by providing @data and @recv_ecc only:
 *   decode_bch(@bch, @data, @len, @recv_ecc, NULL, NULL, @errloc)
 *
 * by providing @recv_ecc and @calc_ecc:
 *   decode_bch(@bch, NULL, @len, @recv_ecc, @calc_ecc, NULL, @errloc)
 *
 * by providing ecc = recv_ecc XOR calc_ecc:
 *   decode_bch(@bch, NULL, @len, NULL, ecc, NULL, @errloc)
 *
 * by providing syndrome results @syn:
 *   decode_bch(@bch, NULL, @len, NULL, NULL, @syn, @errloc)
 *
 * Once decode_bch() has successfully returned with a positive value, error
 * locations returned in array @errloc should be interpreted as follows -
 *
 * if (errloc[n] >= 8*len), then n-th error is located in ecc (no need for
 * data correction)
 *
 * if (errloc[n] < 8*len), then n-th error is located in data and can be
 * corrected with statement data[errloc[n]/8] ^= 1 << (errloc[n] % 8);
 *
 * Note that this function does not perform any data correction by itself, it
 * merely indicates error locations.
 */
int decode_bch(struct bch_control *bch, const uint8_t *data, unsigned int len,
	       const uint8_t *recv_ecc, const uint8_t *calc_ecc,
	       const unsigned int *syn, unsigned int *errloc)
{
	const unsigned int ecc_words = BCH_ECC_WORDS(bch);
	unsigned int nbits;
	int i, err, nroots;
	uint32_t sum;

	/* sanity check: make sure data length can be handled */
	if (8*len > (bch->n-bch->ecc_bits))
		return -EINVAL;

	/* if caller does not provide syndromes, compute them */
	if (!syn) {
		if (!calc_ecc) {
			/* compute received data ecc into an internal buffer */
			if (!data || !recv_ecc)
				return -EINVAL;
			encode_bch(bch, data, len, NULL);
		} else {
			/* load provided calculated ecc */
			load_ecc8(bch, bch->ecc_buf, calc_ecc);
		}
		/* load received ecc or assume it was XORed in calc_ecc */
		if (recv_ecc) {
			load_ecc8(bch, bch->ecc_buf2, recv_ecc);
			/* XOR received and calculated ecc */
			for (i = 0, sum = 0; i < (int)ecc_words; i++) {
				bch->ecc_buf[i] ^= bch->ecc_buf2[i];
				sum |= bch->ecc_buf[i];
			}
			if (!sum)
				/* no error found */
				return 0;
		}
		compute_syndromes(bch, bch->ecc_buf, bch->syn);
		syn = bch->syn;
	}

	err = compute_error_locator_polynomial(bch, syn);
	if (err > 0) {
		nroots = find_poly_roots(bch, 1, bch->elp, errloc);
		if (err != nroots)
			err = -1;
	}
	if (err > 0) {
		/* post-process raw error locations for easier correction */
		nbits = (len*8)+bch->ecc_bits;
		for (i = 0; i < err; i++) {
			if (errloc[i] >= nbits) {
				err = -1;
				break;
			}
			errloc[i] = nbits-1-errloc[i];
			errloc[i] = (errloc[i] & ~7)|(7-(errloc[i] & 7));
		}
	}
	return (err >= 0) ? err : -EBADMSG;
}
EXPORT_SYMBOL_GPL(decode_bch);

#endif
