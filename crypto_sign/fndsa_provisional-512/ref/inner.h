#ifndef FNDSA_INNER_H__
#define FNDSA_INNER_H__

/* ==================================================================== */

#include "fndsa.h"

#include "archflags.h"

#include <stddef.h>
#include <stdint.h>
#include <string.h>

/*
 * Naming conventions
 * ==================
 *
 * All identifiers with external linkage are prefixed with "fndsa_".
 * This header includes convenient macros that map the prefix-less names
 * to the actual identifiers.
 *
 * When a function has an alternate implementation, using special CPU
 * feature "fff" and which is to be used conditionally to a runtime test
 * on the support by the CPU, then the prefix for that alternate
 * implementation is "fndsa_fff_"; again, a convenient macro allows using
 * the name without the "fndsa_" prefix.
 */

/* Define FNDSA_AVX2 to 1 in order to add AVX2 support, 0 otherwise. */
#ifndef FNDSA_AVX2
#if defined __x86_64__ || defined _M_X64 || defined __i386__ || defined _M_IX86
#define FNDSA_AVX2   1
#else
#define FNDSA_AVX2   0
#endif
#endif

/* TARGET_AVX2 is applied to a function definition and allows use of AVX2
   intrinsics in that function. */
#if FNDSA_AVX2
#if defined __GNUC__ || defined __clang__
#define TARGET_AVX2    __attribute__((target("avx2,lzcnt")))
#else
#define TARGET_AVX2
#endif
#else
#define TARGET_AVX2
#endif

/* ALIGN32 is applied to a declarator and will try to make the declared
   object aligned at a 32-byte boundary in memory. */
#if defined __GNUC__ || defined __clang__
#define ALIGN32   __attribute__((aligned(32)))
#elif defined _MSC_VER
#define ALIGN32   __declspec(align(32))
#else
#define ALIGN32
#endif

/* FNDSA_SSE2 is set to 1 if SSE2 is supported (x86 only). This test
   is static (no runtime detection). On 64-bit x86, SSE2 is part of the
   ABI, so FNDSA_SSE2 will always be set. On 32-bit x86, whether SSE2
   is enabled or not depends on the compiler (recent MSVC versions set
   it by default, GCC and Clang on Linux do not). */
#ifndef FNDSA_SSE2
#if ((defined __GNUC__ || defined __clang__) && defined __SSE2__) \
	|| (defined _MSC_VER && defined _M_X64) \
	|| (defined _MSC_VER && defined _M_IX86_FP && _M_IX86_FP >= 2)
#define FNDSA_SSE2   1
#else
#define FNDSA_SSE2   0
#endif
#endif

/* TARGET_SSE2 is applied to a function definition and allows use of SSE2
   intrinsics in that function. */
#if FNDSA_SSE2
#if defined __GNUC__ || defined __clang__
#define TARGET_SSE2    __attribute__((target("sse2")))
#else
#define TARGET_SSE2
#endif
#else
#define TARGET_SSE2
#endif

#if FNDSA_AVX2 || FNDSA_SSE2
#include <immintrin.h>
#if defined __GNUC__ || defined __clang__
#include <x86intrin.h>
#endif
#endif

/* FNDSA_NEON is set to 1 if NEON is supported by the hardware, including
   vectors of 64-bit (double precision) floating-point values (aarch64 only).
   This test is static (no runtime detection). Since NEON is normally part
   of the aarch64 ABI, this should always be set on this architecture. */
#ifndef FNDSA_NEON
#if defined __aarch64__ \
	&& ((defined __ARM_NEON_FP && ((__ARM_NEON_FP & 0x0C) != 0)) \
	 || (!defined __ARM_NEON_FP && ((__ARM_FP & 0x0C) != 0)))
#define FNDSA_NEON   1
#else
#define FNDSA_NEON   0
#endif
#endif

/* TARGET_NEON is applied to a function definition and allows use of NEON
   intrinsics in that function. */
#if FNDSA_NEON
#if defined __GNUC__ || defined __clang__
#define TARGET_NEON    __attribute__((target("neon")))
#else
#define TARGET_NEON
#endif
#else
#define TARGET_NEON
#endif

#if FNDSA_NEON
#include <arm_neon.h>
#endif

/* If FNDSA_NEON_SHA3 is non-zero, then NEON opcodes will be used to
   make two SHAKE256 parallel evaluations during keygen and signing.
   This happens to be slower than the plain code on ARM Cortex-A55
   and Cortex-A76, which is why it is not enabled by default. */
#ifndef FNDSA_NEON_SHA3
#define FNDSA_NEON_SHA3   0
#endif

/* FNDSA_RV64D is set to 1 if the hardware uses the riscv64 architecture
   with the D (64-bit floating-point) extension; this is the case of most
   "general purpose" RISC-V hardware (which is "RV64GC", with "G" being
   a shorthand for I, M, A, F and D). */
#ifndef FNDSA_RV64D
#if defined __riscv && defined __riscv_xlen && __riscv_xlen >= 64 \
	&& defined __riscv_d && defined __riscv_flen && __riscv_flen >= 64
#define FNDSA_RV64D   1
#else
#define FNDSA_RV64D   0
#endif
#endif

/* If FNDSA_DIV_EMU is 1, then floating-point division will use the
   emulated code with only integer operations, instead of the hardware
   abilities. This option impacts only architecture which otherwise use
   the hardware native abilities for divisions; moreover, for platforms
   where SIMD opcodes are used (e.g. SSE2 or NEON), the divisions
   will still be performed with the SIMD unit. The main target for this
   option is RISC-V cores; in particular, the SiFive U74 offers only
   non-constant-time floating-point divisions in hardware. Enabling this
   option (and also FNDSA_SQRT_EMU) increases signing cost by about 25%. */
#ifndef FNDSA_DIV_EMU
#define FNDSA_DIV_EMU   0
#endif

/* If FNDSA_SQRT_EMU is 1, then floating-point square roots will use the
   emulated code with only integer operations, instead of the hardware
   abilities. As with FNDSA_DIV_EMU, this impacts only architecture for
   which the hardware abitilies as used for square roots, but not platforms
   with explicit SIMD opcode support (e.g. SSE2 or NEON). */
#ifndef FNDSA_SQRT_EMU
#define FNDSA_SQRT_EMU   0
#endif

/* If FNDSA_ASM_CORTEXM4 is 1, then the code will use optimized assembly
   routines, under the assumption that it runs on an ARM Cortex-M4
   (specifically an M4F: the hardware floating-point is not used since
   it's only single-precision, but the FP registers are leveraged as
   temporary storage).
   This is not set by default because some of the routines are in
   external assembly-only source files which the Makefile must then
   include explicitly. That Makefile may thus set FNDSA_ASM_CORTEXM4 too. */
#ifndef FNDSA_ASM_CORTEXM4
#define FNDSA_ASM_CORTEXM4   0
#endif

/* Failsafe check: if ARMv7M assembly code is detected, the compiler should
   know about it too. */
#if FNDSA_ASM_CORTEXM4
#if !(defined __ARM_ARCH_7EM__ && defined __ARM_FEATURE_DSP)
#error ARMv7M+DSP assembly code selected but supported.
#endif
#endif

/* Automatically recognize some architectures as being "64-bit", which
   mostly means that we assume that 64-bit shifts are constant-time
   with regard to the shift count. */
#ifndef FNDSA_64
#if defined __x86_64__ || defined _M_X64 \
	|| defined __ia64 || defined __itanium__ || defined _M_IA64 \
	|| defined __powerpc64__ || defined __ppc64__ || defined __PPC64__ \
	|| defined __64BIT__ || defined _LP64 || defined __LP64__ \
	|| defined __sparc64__ \
	|| defined __aarch64__ || defined _M_ARM64 \
	|| defined __mips64
#define FNDSA_64   1
#else
#define FNDSA_64   0
#endif
#endif

/* Recognize little-endian architectures. A little-endian system can
   extract bytes from a SHAKE context more efficiently. If
   FNDSA_LITTLE_ENDIAN is 0, then the platform is not assumed to be
   little-endian (but it still can be). */
#ifndef FNDSA_LITTLE_ENDIAN
#if defined __LITTLE_ENDIAN__ \
	|| (defined __BYTE_ORDER__ && defined __ORDER_LITTLE_ENDIAN__ \
		&& __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__) \
	|| defined _M_IX86 || defined _M_X64 || defined _M_ARM64
#define FNDSA_LITTLE_ENDIAN   1
#else
#define FNDSA_LITTLE_ENDIAN   0
#endif
#endif

/* Some architectures tolerate well unaligned accesses to 64-bit words. */
#ifndef FNDSA_UNALIGNED_64
#if defined __x86_64__ || defined _M_X64 \
	|| defined __i386__ || defined _M_IX86 \
	|| defined __aarch64__ || defined _M_ARM64
#define FNDSA_UNALIGNED_64   1
#else
#define FNDSA_UNALIGNED_64   0
#endif
#endif

/* Some architectures tolerate well unaligned accesses to 16-bit words. */
#ifndef FNDSA_UNALIGNED_16
#if defined __x86_64__ || defined _M_X64 \
	|| defined __i386__ || defined _M_IX86 \
	|| FNDSA_ASM_CORTEXM4 \
	|| defined __aarch64__ || defined _M_ARM64
#define FNDSA_UNALIGNED_16   1
#else
#define FNDSA_UNALIGNED_16   0
#endif
#endif

/* Some MSVC adjustments. */
#if defined _MSC_VER
/* Disable some warnings which are about valid and well-defined operations
   that are often unwanted (but here we want them, e.g. applying the negation
   operator on unsigned types). */
#pragma warning( disable : 4146 )
#pragma warning( disable : 4204 )
#pragma warning( disable : 4244 )
#pragma warning( disable : 4267 )
#pragma warning( disable : 4334 )
/* MSVC does not like the 'restrict' keyword. */
#define restrict
#endif

/* ==================================================================== */
/*
 * SHAKE implementation.
 *
 * A SHAKE context is (re)initialized with shake_init(). After
 * initialization, the context is in input mode, and accepts data with
 * shake_inject() calls. shake_flip() switches the context to output mode,
 * allowing data to be extracted with shake_extract() calls.
 *
 * Callers MUST NOT inject data in a context which has been flipped to
 * output mode. Callers MUST NOT extract data from a context which has
 * not been flipped to output mode. Reinitializing the context also sets
 * it back to input mode.
 *
 * A context contains no external or internal pointers, and can thus be
 * moved, cloned, or abandoned.
 */

typedef struct {
	uint64_t A[25];
	unsigned dptr, rate;
} shake_context;

#define shake_init      fndsa_shake_init
#define shake_inject    fndsa_shake_inject
#define shake_flip      fndsa_shake_flip
#define shake_extract   fndsa_shake_extract

/* Initialize context, size = 128 or 256 (for SHAKE128 or SHAKE256). */
void shake_init(shake_context *sc, unsigned size);
/* Inject some bytes in context. */
void shake_inject(shake_context *sc, const void *in, size_t len);
/* Flip context from input to output mode. */
void shake_flip(shake_context *sc);
/* Extract some bytes from context. If out is NULL, then len bytes are
   still virtually extracted, but discarded.
   In systems with little-endian encoding, the discarded bytes
   can still be obtained from the context; this is used for saving some
   RAM (especially stack space) on embedded systems. */
void shake_extract(shake_context *sc, void *out, size_t len);

/* Get the next byte from a SHAKE context. */
static inline uint8_t
shake_next_u8(shake_context *sc)
{
	if (sc->dptr == sc->rate) {
		uint8_t x;
		shake_extract(sc, &x, 1);
		return x;
	}
#if FNDSA_LITTLE_ENDIAN
	uint8_t *d = (uint8_t *)(void *)sc;
	return d[sc->dptr ++];
#else
	uint8_t x = (uint8_t)(sc->A[sc->dptr >> 3] >> ((sc->dptr & 7) << 3));
	sc->dptr ++;
	return x;
#endif
}

/* Get the next 16-bit word from SHAKE. */
static inline unsigned
shake_next_u16(shake_context *sc)
{
	if (sc->dptr + 1 >= sc->rate) {
		uint8_t x[2];
		shake_extract(sc, x, 2);
		return (unsigned)x[0] | ((unsigned)x[1] << 8);
	}
#if FNDSA_LITTLE_ENDIAN
	uint8_t *d = (uint8_t *)(void *)sc;
#if FNDSA_UNALIGNED_16
	unsigned v = *(uint16_t *)(d + sc->dptr);
#else
	unsigned v = (unsigned)d[sc->dptr] | ((unsigned)d[sc->dptr + 1] << 8);
#endif
	sc->dptr += 2;
	return v;
#else
	unsigned x0 = (uint8_t)(sc->A[sc->dptr >> 3] >> ((sc->dptr & 7) << 3));
	sc->dptr ++;
	unsigned x1 = (uint8_t)(sc->A[sc->dptr >> 3] >> ((sc->dptr & 7) << 3));
	sc->dptr ++;
	return x0 | (x1 << 8);
#endif
}

/* Get the next 64-bit word from SHAKE. */
static inline uint64_t
shake_next_u64(shake_context *sc)
{
	if ((sc->dptr + 7) >= sc->rate) {
#if FNDSA_LITTLE_ENDIAN
		uint64_t v;
		shake_extract(sc, &v, 8);
		return v;
#else
		uint8_t x[8];
		shake_extract(sc, x, 8);
		return (uint64_t)x[0]
			| ((uint64_t)x[1] << 8)
			| ((uint64_t)x[2] << 16)
			| ((uint64_t)x[3] << 24)
			| ((uint64_t)x[4] << 32)
			| ((uint64_t)x[5] << 40)
			| ((uint64_t)x[6] << 48)
			| ((uint64_t)x[7] << 56);
#endif
	}
	uint64_t x;
#if FNDSA_LITTLE_ENDIAN && FNDSA_UNALIGNED_64
	x = *(uint64_t *)((uint8_t *)(void *)sc + sc->dptr);
#else
	size_t j = sc->dptr >> 3;
	unsigned n = sc->dptr & 7;
	if (n == 0) {
		x = sc->A[j];
	} else {
		x = sc->A[j] >> (n << 3);
		x |= sc->A[j + 1] << (64 - (n << 3));
	}
#endif
	sc->dptr += 8;
	return x;
}

/* By default, we use a simple SHAKE256 for internal PRNG needs
   (in keygen to generate (f,g), in signing for the Gaussian sampling).
   If FNDSA_SHAKE256X4 is non-zero, then SHAKE256x4 is used: it is a
   PRNG consisting of four SHAKE256 running in parallel, with interleaved
   outputs. This has two main effects:
     - On x86 with AVX2 support, this makes signing faster (by about 20%).
     - It increases stack usage by about 1.1 kB, which can be a concern
       for small embedded systems (e.g. ARM Cortex-M4).
   It otherwise has no real perceivable effect, except that (of course)
   it changes the exact key pairs and signature values obtained from a
   given seed. */
#ifndef FNDSA_SHAKE256X4
#define FNDSA_SHAKE256X4   0
#endif

#if FNDSA_SHAKE256X4
/*
 * SHAKE256x4 is a PRNG based on SHAKE256; it runs four SHAKE256 instances
 * in parallel, interleaving their outputs with 64-bit granularity. The
 * four instances are initialized with a common seed, followed by a single
 * byte of value 0x00, 0x01, 0x02 or 0x03, depending on the SHAKE instance.
 */

typedef struct {
	uint64_t state[100];
	uint8_t buf[4 * 136];
	unsigned ptr;
#if FNDSA_AVX2
	int use_avx2;
#endif
} shake256x4_context;

/* Initialize a SHAKE256x4 context from a given seed.
   WARNING: seed length MUST NOT exceed 134 bytes. */
#define shake256x4_init   fndsa_shake256x4_init
void shake256x4_init(shake256x4_context *sc, const void *seed, size_t seed_len);
/* Refill the SHAKE256x4 output buffer. */
#define shake256x4_refill   fndsa_shake256x4_refill
void shake256x4_refill(shake256x4_context *sc);

/* Get the next byte of pseudorandom output. */
static inline uint8_t
shake256x4_next_u8(shake256x4_context *sc)
{
	if (sc->ptr >= sizeof sc->buf) {
		shake256x4_refill(sc);
	}
	return sc->buf[sc->ptr ++];
}

/* Get the next 16-bit word of pseudorandom output. */
static inline unsigned
shake256x4_next_u16(shake256x4_context *sc)
{
	if (sc->ptr >= (sizeof sc->buf) - 1) {
		shake256x4_refill(sc);
	}
	unsigned x = (unsigned)sc->buf[sc->ptr]
		| ((unsigned)sc->buf[sc->ptr + 1] << 8);
	sc->ptr += 2;
	return x;
}

/* Get the next 64-bit word of pseudorandom output. */
static inline uint64_t
shake256x4_next_u64(shake256x4_context *sc)
{
	if (sc->ptr >= (sizeof sc->buf) - 7) {
		shake256x4_refill(sc);
	}
	uint64_t x = (uint64_t)sc->buf[sc->ptr]
		| ((uint64_t)sc->buf[sc->ptr + 1] << 8)
		| ((uint64_t)sc->buf[sc->ptr + 2] << 16)
		| ((uint64_t)sc->buf[sc->ptr + 3] << 24)
		| ((uint64_t)sc->buf[sc->ptr + 4] << 32)
		| ((uint64_t)sc->buf[sc->ptr + 5] << 40)
		| ((uint64_t)sc->buf[sc->ptr + 6] << 48)
		| ((uint64_t)sc->buf[sc->ptr + 7] << 56);
	sc->ptr += 8;
	return x;
}
#endif

/*
 * SHA-3 implementation.
 * This is a variation on the SHAKE implementation, which implements
 * SHA-3 with outputs of 224, 256, 384 and 512 bits. The output size
 * (in bits) is provided as parameter to sha3_init(). Input is provided
 * with sha3_update(). The output is computed with sha3_close(); this
 * function also reinitializes the context.
 */
typedef shake_context sha3_context;

#define sha3_init      fndsa_sha3_init
#define sha3_update    fndsa_sha3_update
#define sha3_close     fndsa_sha3_close

/* Initialize context, size = 224, 256, 384 or 512 */
void sha3_init(sha3_context *sc, unsigned size);
/* Inject some bytes in context. */
void sha3_update(sha3_context *sc, const void *in, size_t len);
/* Compute the output and reinitialize the context. */
void sha3_close(sha3_context *sc, void *out);

/* ==================================================================== */
/*
 * Encoding/decoding primitives.
 */

#define trim_i8_encode   fndsa_trim_i8_encode
#define trim_i8_decode   fndsa_trim_i8_decode
#define mqpoly_encode    fndsa_mqpoly_encode
#define mqpoly_decode    fndsa_mqpoly_decode
#define comp_encode      fndsa_comp_encode
#define comp_decode      fndsa_comp_decode

/* Encode small polynomial f into output buffer d, nbits by element.
   Total encoded size (in bytes) is returned. The encoding MUST use
   an integral number of bytes (i.e. nbits*2^logn must be a multiple of 8). */
size_t trim_i8_encode(unsigned logn,
	const int8_t *f, unsigned nbits, uint8_t *d);

/* Decode small polynomial f from input buffer d, nbits by element.
   Returned value is the number of read bytes on success, 0 on error (an
   error is reported if an element has value -2^(nbits-1), which is not
   allowed). */
size_t trim_i8_decode(unsigned logn,
	const uint8_t *d, int8_t *f, unsigned nbits);

/* Encode polynomial h (external representation, values are in [0,q-1])
   into bytes, 14 bits per value. Total encoded size (in bytes) is
   returned. */
size_t mqpoly_encode(unsigned logn, const uint16_t *h, uint8_t *d);

/* Decode polynomial h (external representation, values are in [0,q-1])
   from bytes, 14 bits per value. Total encoded size (in bytes) is
   returned. On error (a value is out-of-range), 0 is returned. */
size_t mqpoly_decode(unsigned logn, const uint8_t *d, uint16_t *h);

/* Encode polynomial s into destination buffer d (of size dlen bytes),
   using compressed (Golomb-Rice) format. If any of the source values is
   outside of [-2047,+2047], this function fails and returns 0. If the
   destination buffer is not large enough to receive all values, this
   function fails and returns 0. Otherwise, this function succeeds and
   returns 1. On success, the entirety of the d buffer is written
   (padding bits/bytes of values 0 are appended if necessary). */
int comp_encode(unsigned logn, const int16_t *s, uint8_t *d, size_t dlen);

/* Decode polynomial s from buffer d (of size dlen bytes), using
   compressed (Golomb-Rice) format. Returned value is 1 on success, 0 on
   failure. A failure is reported in any of the following cases:
     - Source does not contain enough bytes for all the polynomial values.
     - An invalid value encoding is encountered.
     - Unused bits/bytes in the source buffer are not all zeros.
   Validly encoded values are in [-2047,+2047]. Each such value has a
   single valid (canonical) encoding. */
int comp_decode(unsigned logn, const uint8_t *d, size_t dlen, int16_t *s);

/* ==================================================================== */
/*
 * Computations modulo q = 12289.
 *
 * Polynomials are held in arrays of integers. The following
 * representations are used:
 *   name       type     range
 *  --------------------------------
 *   small     int8_t    [-127,+127]
 *   signed   uint16_t   [-2047,+2047] (int16_t cast to uint16_t)
 *   ext      uint16_t   [0,q-1]
 *   int      uint16_t   (internal)
 *   ntt      uint16_t   (internal)
 *
 * The internal convention may change depending on the implementation.
 * For most architectures, internal representation uses [1,q] (i.e. zero
 * is represented by q instead of 0). The 32-bit ARMv7 code (for ARM
 * Cortex-M4) uses [0,q] (i.e. zero can be represented by either 0 or q).
 */

#define mqpoly_small_to_int           fndsa_mqpoly_small_to_int
#define mqpoly_signed_to_int          fndsa_mqpoly_signed_to_int
#define mqpoly_int_to_small           fndsa_mqpoly_int_to_small
#define mqpoly_ext_to_int             fndsa_mqpoly_ext_to_int
#define mqpoly_int_to_ext             fndsa_mqpoly_int_to_ext
#define mqpoly_int_to_ntt             fndsa_mqpoly_int_to_ntt
#define mqpoly_ntt_to_int             fndsa_mqpoly_ntt_to_int
#define mqpoly_mul_ntt                fndsa_mqpoly_mul_ntt
#define mqpoly_div_ntt                fndsa_mqpoly_div_ntt
#define mqpoly_sub                    fndsa_mqpoly_sub
#define mqpoly_is_invertible          fndsa_mqpoly_is_invertible
#define mqpoly_div_small              fndsa_mqpoly_div_small
#define mqpoly_sqnorm_ext             fndsa_mqpoly_sqnorm_ext
#define mqpoly_sqnorm_signed          fndsa_mqpoly_sqnorm_signed
#define mqpoly_sqnorm_is_acceptable   fndsa_mqpoly_sqnorm_is_acceptable
#define mq_GM                         fndsa_mq_GM
#define mq_iGM                        fndsa_mq_iGM

/* Convert a polynomial from small (f) to int (d). */
void mqpoly_small_to_int(unsigned logn, const int8_t *f, uint16_t *d);

/* Convert a polynomial from signed to int (in-place). */
void mqpoly_signed_to_int(unsigned logn, uint16_t *d);

/* Convert a polynomial from int (d) to small (f).
   Returns 1 if all values are in [-127,+127], 0 otherwise. */
int mqpoly_int_to_small(unsigned logn, const uint16_t *d, int8_t *f);

#if FNDSA_ASM_CORTEXM4
static inline void
mqpoly_ext_to_int(unsigned logn, uint16_t *d)
{
	(void)logn;
	(void)d;
}
#else
/* Convert a polynomial from ext to int (in-place). */
void mqpoly_ext_to_int(unsigned logn, uint16_t *d);
#endif

/* Convert a polynomial from int to ext (in-place). */
void mqpoly_int_to_ext(unsigned logn, uint16_t *d);

/* Convert a polynomial from int to ntt (in-place). */
void mqpoly_int_to_ntt(unsigned logn, uint16_t *d);

/* Convert a polynomial from ntt to int (in-place). */
void mqpoly_ntt_to_int(unsigned logn, uint16_t *d);

/* Multiply polynomial a by polynomial b (both in ntt representation). */
void mqpoly_mul_ntt(unsigned logn, uint16_t *a, const uint16_t *b);

/* Divide polynomial a by polynomial b (both in ntt representation).
   If b is not invertible, then the corresponding coefficients are set
   to zero in the output ntt representation.
   Return value is 1 on success (b is invertible), 0 otherwise. */
int mqpoly_div_ntt(unsigned logn, uint16_t *a, const uint16_t *b);

/* Subtract polynomial b from polynomial a (both must be in int
   representation, or both must be in ntt representation). */
void mqpoly_sub(unsigned logn, uint16_t *a, const uint16_t *b);

/* Check whether the small polynomial f is invertible. tmp[] must
   have n elements. Returned value is 1 if invertible, 0 otherwise. */
int mqpoly_is_invertible(unsigned logn, const int8_t *f, uint16_t *tmp);

/* Compute h = g/f. Output is in external representation (coefficients
   in [0,q-1]). This function assumes that f is invertible. tmp[] must
   have n elements. */
void mqpoly_div_small(unsigned logn, const int8_t *g, const int8_t *f,
	uint16_t *h, uint16_t *tmp);

/* Compute the squared norm of a polynomial (in external representation).
   The squared norm includes an implicit normalization to [-q/2,+q/2].
   If the value exceeds 2^31-1 then 2^32-1 is returned. */
uint32_t mqpoly_sqnorm_ext(unsigned logn, const uint16_t *a);

/* Compute the squared norm of a polynomial (in signed representation).
   Since the signed representation assumes that all coefficients are
   in [-2047,+2047], and the degree is at most 2^10 = 1024, the largest
   possible squared norm is 4290774016, which cannot overflow the 32-bit
   return type. */
uint32_t mqpoly_sqnorm_signed(unsigned logn, const uint16_t *a);

/* Check whether a given squared norm (for a signature vector (s1,s2))
   is acceptable. */
int mqpoly_sqnorm_is_acceptable(unsigned logn, uint32_t norm);

/* Tables of constants for NTT and inverse NTT modulo q. */
extern const uint16_t mq_GM[];
extern const uint16_t mq_iGM[];

#if FNDSA_AVX2
#define avx2_mqpoly_small_to_int          fndsa_avx2_mqpoly_small_to_int
#define avx2_mqpoly_signed_to_int         fndsa_avx2_mqpoly_signed_to_int
#define avx2_mqpoly_int_to_small          fndsa_avx2_mqpoly_int_to_small
#define avx2_mqpoly_ext_to_int            fndsa_avx2_mqpoly_ext_to_int
#define avx2_mqpoly_int_to_ext            fndsa_avx2_mqpoly_int_to_ext
#define avx2_mqpoly_int_to_ntt            fndsa_avx2_mqpoly_int_to_ntt
#define avx2_mqpoly_ntt_to_int            fndsa_avx2_mqpoly_ntt_to_int
#define avx2_mqpoly_mul_ntt               fndsa_avx2_mqpoly_mul_ntt
#define avx2_mqpoly_div_ntt               fndsa_avx2_mqpoly_div_ntt
#define avx2_mqpoly_sub                   fndsa_avx2_mqpoly_sub
#define avx2_mqpoly_is_invertible         fndsa_avx2_mqpoly_is_invertible
#define avx2_mqpoly_div_small             fndsa_avx2_mqpoly_div_small
#define avx2_mqpoly_sqnorm_ext            fndsa_avx2_mqpoly_sqnorm_ext
#define avx2_mqpoly_sqnorm_signed         fndsa_avx2_mqpoly_sqnorm_signed
void avx2_mqpoly_small_to_int(unsigned logn, const int8_t *f, uint16_t *d);
void avx2_mqpoly_signed_to_int(unsigned logn, uint16_t *d);
int avx2_mqpoly_int_to_small(unsigned logn, const uint16_t *d, int8_t *f);
void avx2_mqpoly_ext_to_int(unsigned logn, uint16_t *d);
void avx2_mqpoly_int_to_ext(unsigned logn, uint16_t *d);
void avx2_mqpoly_int_to_ntt(unsigned logn, uint16_t *d);
void avx2_mqpoly_ntt_to_int(unsigned logn, uint16_t *d);
void avx2_mqpoly_mul_ntt(unsigned logn, uint16_t *a, const uint16_t *b);
int avx2_mqpoly_div_ntt(unsigned logn, uint16_t *a, const uint16_t *b);
void avx2_mqpoly_sub(unsigned logn, uint16_t *a, const uint16_t *b);
int avx2_mqpoly_is_invertible(unsigned logn, const int8_t *f, uint16_t *tmp);
void avx2_mqpoly_div_small(unsigned logn, const int8_t *f, const int8_t *g,
	uint16_t *h, uint16_t *tmp);
uint32_t avx2_mqpoly_sqnorm_ext(unsigned logn, const uint16_t *a);
uint32_t avx2_mqpoly_sqnorm_signed(unsigned logn, const uint16_t *a);
#endif

/* ==================================================================== */
/*
 * Utility functions.
 */

/* Hash a (pre-hashed) message into a polynomial.
    logn              degree (logarithmic, 2 to 10)
    nonce             40-byte random nonce
    hashed_vrfy_key   64-byte hashed public key (with SHAKE256)
    ctx, ctx_len      domain separation context (at most 255 bytes)
    hash_id           hash identifier
    hv, hv_len        pre-hashed message (raw if hash_id = FNDSA_HASH_ID_RAW)
    c                 output polynomial  */
#define hash_to_point   fndsa_hash_to_point
void hash_to_point(unsigned logn,
	const uint8_t *nonce, const uint8_t *hashed_vrfy_key,
	const void *ctx, size_t ctx_len,
	const char *hash_id, const void *hv, size_t hv_len,
	uint16_t *c);

#if FNDSA_AVX2
#define has_avx2   fndsa_has_avx2
/* Check for AVX2 support by the current CPU. */
int has_avx2(void);
#endif

/* Expand the top bit of a 32-bit word into a full 32-bit mask (i.e. return
   0xFFFFFFFF if x >= 0x80000000, or 0x00000000 otherwise). */
static inline uint32_t
tbmask(uint32_t x)
{
	return (uint32_t)(*(int32_t *)&x >> 31);
}

/* Get the number of leading zeros in a 32-bit value. */
static inline unsigned
lzcnt(uint32_t x)
{
#if FNDSA_ASM_CORTEXM4
	unsigned r;
	__asm__ ("clz %0, %1" : "=r" (r) : "r" (x));
	return r;
#else
	/* TODO: optimize on x86 and ARMv8 */
	uint32_t m = tbmask((x >> 16) - 1);
	uint32_t s = m & 16;
	x = (x >> 16) ^ (m & (x ^ (x >> 16)));
	m = tbmask((x >>  8) - 1);
	s |= m & 8;
	x = (x >>  8) ^ (m & (x ^ (x >>  8)));
	m = tbmask((x >>  4) - 1);
	s |= m & 4;
	x = (x >>  4) ^ (m & (x ^ (x >>  4)));
	m = tbmask((x >>  2) - 1);
	s |= m & 2;
	x = (x >>  2) ^ (m & (x ^ (x >>  2)));
	return (unsigned)(s + ((2 - x) & tbmask(x - 3)));
#endif
}

/* Same as lzcnt(), but the caller ensures that the operand is non-zero.
   TODO: on x86, this can be implemented with the bsr opcode. */
#define lzcnt_nonzero   lzcnt

/* Obtain fresh randomness from the operating system. This function shall
   ensure that the requested entropy is achieved. If the operating system
   does not have a secure random source, or if that source fails, then
   this function will fail and return 0. On success it returns 1. */
#define sysrng   fndsa_sysrng
int sysrng(void *dst, size_t len);

/* ==================================================================== */

#endif
