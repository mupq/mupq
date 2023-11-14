/*
 * Reference ISO-C11 Implementation of CROSS.
 *
 * @version 1.1 (May 2023)
 *
 * @author: Patrick Karl <patrick.karl@tum.de>
 *
 *
 * SHA2 self-contained implementation derived from :
 *
 * LibTomCrypt, modular cryptographic library
 *
 * @author: Tom St Denis
 *
 * originally placed in the public domain by the author.
 *
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "sha2.h"

/******************************* Various logical functions *******************************/
#define MIN(x, y) ( ((x)<(y)) ? (x) : (y) )
#define CONST64(n) n ## ULL

#define RORc(x, y) ( ((((uint32_t)(x)&0xFFFFFFFFUL)>>(uint32_t)((y)&31)) | \
            ((uint32_t)(x)<<(uint32_t)((32-((y)&31))&31))) & 0xFFFFFFFFUL)
#define ROR64c(x, y) ( ((((x)&CONST64(0xFFFFFFFFFFFFFFFF))>>((uint64_t)(y)&CONST64(63))) | \
            ((x)<<(((uint64_t)64-((y)&63))&63))) & CONST64(0xFFFFFFFFFFFFFFFF))

#define STORE64H(x, y)                                                                     \
do { (y)[0] = (unsigned char)(((x)>>56)&255); (y)[1] = (unsigned char)(((x)>>48)&255);     \
    (y)[2] = (unsigned char)(((x)>>40)&255); (y)[3] = (unsigned char)(((x)>>32)&255);     \
    (y)[4] = (unsigned char)(((x)>>24)&255); (y)[5] = (unsigned char)(((x)>>16)&255);     \
    (y)[6] = (unsigned char)(((x)>>8)&255); (y)[7] = (unsigned char)((x)&255); } while(0)

#define LOAD64H(x, y)                                                      \
do { x = (((uint64_t)((y)[0] & 255))<<56)|(((uint64_t)((y)[1] & 255))<<48) | \
        (((uint64_t)((y)[2] & 255))<<40)|(((uint64_t)((y)[3] & 255))<<32) | \
        (((uint64_t)((y)[4] & 255))<<24)|(((uint64_t)((y)[5] & 255))<<16) | \
        (((uint64_t)((y)[6] & 255))<<8)|(((uint64_t)((y)[7] & 255))); } while(0)

#define STORE32H(x, y)                                                                     \
do { (y)[0] = (unsigned char)(((x)>>24)&255); (y)[1] = (unsigned char)(((x)>>16)&255);   \
    (y)[2] = (unsigned char)(((x)>>8)&255); (y)[3] = (unsigned char)((x)&255); } while(0)

#define LOAD32H(x, y)                            \
do { x = ((uint32_t)((y)[0] & 255)<<24) | \
        ((uint32_t)((y)[1] & 255)<<16) | \
        ((uint32_t)((y)[2] & 255)<<8)  | \
        ((uint32_t)((y)[3] & 255)); } while(0)

#define Ch(x,y,z)       (z ^ (x & (y ^ z)))
#define Maj(x,y,z)      (((x | y) & z) | (x & y))
#define S_256(x, n)     RORc((x),(n))
#define R_256(x, n)     (((x)&0xFFFFFFFFUL)>>(n))
#define Sigma0_256(x)   (S_256(x, 2) ^ S_256(x, 13) ^ S_256(x, 22))
#define Sigma1_256(x)   (S_256(x, 6) ^ S_256(x, 11) ^ S_256(x, 25))
#define Gamma0_256(x)   (S_256(x, 7) ^ S_256(x, 18) ^ R_256(x, 3))
#define Gamma1_256(x)   (S_256(x, 17) ^ S_256(x, 19) ^ R_256(x, 10))

#define S_512(x, n)     ROR64c(x, n)
#define R_512(x, n)     (((x)&CONST64(0xFFFFFFFFFFFFFFFF))>>((uint64_t)n))
#define Sigma0_512(x)   (S_512(x, 28) ^ S_512(x, 34) ^ S_512(x, 39))
#define Sigma1_512(x)   (S_512(x, 14) ^ S_512(x, 18) ^ S_512(x, 41))
#define Gamma0_512(x)   (S_512(x, 1) ^ S_512(x, 8) ^ R_512(x, 7))
#define Gamma1_512(x)   (S_512(x, 19) ^ S_512(x, 61) ^ R_512(x, 6))


/********************************* SHA-256 ************************************************/
/**
   Initialize the hash state
   @param state     The hash state to be initialized
   @return 0 if successful
*/
int sha256_init(sha256_state *state)
{
    state->curlen = 0;
    state->length = 0;
    state->state[0] = 0x6A09E667UL;
    state->state[1] = 0xBB67AE85UL;
    state->state[2] = 0x3C6EF372UL;
    state->state[3] = 0xA54FF53AUL;
    state->state[4] = 0x510E527FUL;
    state->state[5] = 0x9B05688CUL;
    state->state[6] = 0x1F83D9ABUL;
    state->state[7] = 0x5BE0CD19UL;
    return 0;
}


/* compress 512-bits */
int s_sha256_compress(sha256_state *state, const unsigned char *buf)
{
    uint32_t S[8], W[64], t0, t1;
    int i;

    /* copy state into S */
    for (i = 0; i < 8; i++) {
        S[i] = state->state[i];
    }

    /* copy the state into 512-bits into W[0..15] */
    for (i = 0; i < 16; i++) {
        LOAD32H(W[i], buf + (4*i));
    }

    /* fill W[16..63] */
    for (i = 16; i < 64; i++) {
        W[i] = Gamma1_256(W[i - 2]) + W[i - 7] + Gamma0_256(W[i - 15]) + W[i - 16];
    }

    /* Compress */
#define RND(a,b,c,d,e,f,g,h,i,ki)                    \
     t0 = h + Sigma1_256(e) + Ch(e, f, g) + ki + W[i];   \
     t1 = Sigma0_256(a) + Maj(a, b, c);                  \
     d += t0;                                        \
     h  = t0 + t1;

    RND(S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],0,0x428a2f98);
    RND(S[7],S[0],S[1],S[2],S[3],S[4],S[5],S[6],1,0x71374491);
    RND(S[6],S[7],S[0],S[1],S[2],S[3],S[4],S[5],2,0xb5c0fbcf);
    RND(S[5],S[6],S[7],S[0],S[1],S[2],S[3],S[4],3,0xe9b5dba5);
    RND(S[4],S[5],S[6],S[7],S[0],S[1],S[2],S[3],4,0x3956c25b);
    RND(S[3],S[4],S[5],S[6],S[7],S[0],S[1],S[2],5,0x59f111f1);
    RND(S[2],S[3],S[4],S[5],S[6],S[7],S[0],S[1],6,0x923f82a4);
    RND(S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[0],7,0xab1c5ed5);
    RND(S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],8,0xd807aa98);
    RND(S[7],S[0],S[1],S[2],S[3],S[4],S[5],S[6],9,0x12835b01);
    RND(S[6],S[7],S[0],S[1],S[2],S[3],S[4],S[5],10,0x243185be);
    RND(S[5],S[6],S[7],S[0],S[1],S[2],S[3],S[4],11,0x550c7dc3);
    RND(S[4],S[5],S[6],S[7],S[0],S[1],S[2],S[3],12,0x72be5d74);
    RND(S[3],S[4],S[5],S[6],S[7],S[0],S[1],S[2],13,0x80deb1fe);
    RND(S[2],S[3],S[4],S[5],S[6],S[7],S[0],S[1],14,0x9bdc06a7);
    RND(S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[0],15,0xc19bf174);
    RND(S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],16,0xe49b69c1);
    RND(S[7],S[0],S[1],S[2],S[3],S[4],S[5],S[6],17,0xefbe4786);
    RND(S[6],S[7],S[0],S[1],S[2],S[3],S[4],S[5],18,0x0fc19dc6);
    RND(S[5],S[6],S[7],S[0],S[1],S[2],S[3],S[4],19,0x240ca1cc);
    RND(S[4],S[5],S[6],S[7],S[0],S[1],S[2],S[3],20,0x2de92c6f);
    RND(S[3],S[4],S[5],S[6],S[7],S[0],S[1],S[2],21,0x4a7484aa);
    RND(S[2],S[3],S[4],S[5],S[6],S[7],S[0],S[1],22,0x5cb0a9dc);
    RND(S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[0],23,0x76f988da);
    RND(S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],24,0x983e5152);
    RND(S[7],S[0],S[1],S[2],S[3],S[4],S[5],S[6],25,0xa831c66d);
    RND(S[6],S[7],S[0],S[1],S[2],S[3],S[4],S[5],26,0xb00327c8);
    RND(S[5],S[6],S[7],S[0],S[1],S[2],S[3],S[4],27,0xbf597fc7);
    RND(S[4],S[5],S[6],S[7],S[0],S[1],S[2],S[3],28,0xc6e00bf3);
    RND(S[3],S[4],S[5],S[6],S[7],S[0],S[1],S[2],29,0xd5a79147);
    RND(S[2],S[3],S[4],S[5],S[6],S[7],S[0],S[1],30,0x06ca6351);
    RND(S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[0],31,0x14292967);
    RND(S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],32,0x27b70a85);
    RND(S[7],S[0],S[1],S[2],S[3],S[4],S[5],S[6],33,0x2e1b2138);
    RND(S[6],S[7],S[0],S[1],S[2],S[3],S[4],S[5],34,0x4d2c6dfc);
    RND(S[5],S[6],S[7],S[0],S[1],S[2],S[3],S[4],35,0x53380d13);
    RND(S[4],S[5],S[6],S[7],S[0],S[1],S[2],S[3],36,0x650a7354);
    RND(S[3],S[4],S[5],S[6],S[7],S[0],S[1],S[2],37,0x766a0abb);
    RND(S[2],S[3],S[4],S[5],S[6],S[7],S[0],S[1],38,0x81c2c92e);
    RND(S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[0],39,0x92722c85);
    RND(S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],40,0xa2bfe8a1);
    RND(S[7],S[0],S[1],S[2],S[3],S[4],S[5],S[6],41,0xa81a664b);
    RND(S[6],S[7],S[0],S[1],S[2],S[3],S[4],S[5],42,0xc24b8b70);
    RND(S[5],S[6],S[7],S[0],S[1],S[2],S[3],S[4],43,0xc76c51a3);
    RND(S[4],S[5],S[6],S[7],S[0],S[1],S[2],S[3],44,0xd192e819);
    RND(S[3],S[4],S[5],S[6],S[7],S[0],S[1],S[2],45,0xd6990624);
    RND(S[2],S[3],S[4],S[5],S[6],S[7],S[0],S[1],46,0xf40e3585);
    RND(S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[0],47,0x106aa070);
    RND(S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],48,0x19a4c116);
    RND(S[7],S[0],S[1],S[2],S[3],S[4],S[5],S[6],49,0x1e376c08);
    RND(S[6],S[7],S[0],S[1],S[2],S[3],S[4],S[5],50,0x2748774c);
    RND(S[5],S[6],S[7],S[0],S[1],S[2],S[3],S[4],51,0x34b0bcb5);
    RND(S[4],S[5],S[6],S[7],S[0],S[1],S[2],S[3],52,0x391c0cb3);
    RND(S[3],S[4],S[5],S[6],S[7],S[0],S[1],S[2],53,0x4ed8aa4a);
    RND(S[2],S[3],S[4],S[5],S[6],S[7],S[0],S[1],54,0x5b9cca4f);
    RND(S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[0],55,0x682e6ff3);
    RND(S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],56,0x748f82ee);
    RND(S[7],S[0],S[1],S[2],S[3],S[4],S[5],S[6],57,0x78a5636f);
    RND(S[6],S[7],S[0],S[1],S[2],S[3],S[4],S[5],58,0x84c87814);
    RND(S[5],S[6],S[7],S[0],S[1],S[2],S[3],S[4],59,0x8cc70208);
    RND(S[4],S[5],S[6],S[7],S[0],S[1],S[2],S[3],60,0x90befffa);
    RND(S[3],S[4],S[5],S[6],S[7],S[0],S[1],S[2],61,0xa4506ceb);
    RND(S[2],S[3],S[4],S[5],S[6],S[7],S[0],S[1],62,0xbef9a3f7);
    RND(S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[0],63,0xc67178f2);
#undef RND

    /* feedback */
    for (i = 0; i < 8; i++) {
        state->state[i] = state->state[i] + S[i];
    }
    return 0;
}


/**
   Process a block of memory though the hash
   @param state     The hash state
   @param in        The data to hash
   @param inlen     The length of the data (bytes)
   @return 0 if successful
*/
int sha256_process(sha256_state *state, const unsigned char *in, unsigned long inlen)
{
    unsigned long n;

    if (state->curlen > sizeof(state->buf)) {
       return 1;
    }
    if (((state->length + inlen * 8) < state->length) || ((inlen * 8) < inlen)) {
      return 1;
    }
    while (inlen > 0) {
        if (state->curlen == 0 && inlen >= 64) {
           if (s_sha256_compress(state, in) != 0) {
              return 1;
           }
           state->length += 64 * 8;
           in             += 64;
           inlen          -= 64;
        } else {
           n = MIN(inlen, (64 - state->curlen));
           memcpy(state->buf + state->curlen, in, (size_t)n);
           state->curlen += n;
           in             += n;
           inlen          -= n;
           if (state->curlen == 64) {
              if (s_sha256_compress(state, state->buf) != 0) {
                 return 1;
              }
               state->length += 8*64;
               state->curlen = 0;
           }
       }
    }
    return 0;
}


/**
   Terminate the hash to get the digest
   @param state     The hash state
   @param out [out] The destination of the hash (32 bytes)
   @return 0 if successful
*/
int sha256_done(sha256_state *state, unsigned char *out)
{
    int i;

    if (state->curlen >= sizeof(state->buf)) {
       return 1;
    }

    /* increase the length of the message */
    state->length += state->curlen * 8;

    /* append the '1' bit */
    state->buf[state->curlen++] = (unsigned char)0x80;

    /* if the length is currently above 56 bytes we append zeros
     * then compress.  Then we can fall back to padding zeros and length
     * encoding like normal.
     */
    if (state->curlen > 56) {
        while (state->curlen < 64) {
            state->buf[state->curlen++] = (unsigned char)0;
        }
        s_sha256_compress(state, state->buf);
        state->curlen = 0;
    }

    /* pad upto 56 bytes of zeroes */
    while (state->curlen < 56) {
        state->buf[state->curlen++] = (unsigned char)0;
    }
    
    /* store length */
    STORE64H(state->length, state->buf+56);
    s_sha256_compress(state, state->buf);
  
    /* copy output */
    for (i = 0; i < 8; i++) {
        STORE32H(state->state[i], out+(4*i));
    }

    return 0;
}



/********************************* SHA-512 ************************************************/

/**
   Initialize the hash state
   @param state     The hash state to be initialized
   @return 0 if successful
*/
int sha512_init(sha512_state *state)
{
    state->curlen = 0;
    state->length = 0;
    state->state[0] = CONST64(0x6a09e667f3bcc908);
    state->state[1] = CONST64(0xbb67ae8584caa73b);
    state->state[2] = CONST64(0x3c6ef372fe94f82b);
    state->state[3] = CONST64(0xa54ff53a5f1d36f1);
    state->state[4] = CONST64(0x510e527fade682d1);
    state->state[5] = CONST64(0x9b05688c2b3e6c1f);
    state->state[6] = CONST64(0x1f83d9abfb41bd6b);
    state->state[7] = CONST64(0x5be0cd19137e2179);
    return 0;
}


/* the K array */
static const uint64_t K_512[80] = {
    CONST64(0x428a2f98d728ae22), CONST64(0x7137449123ef65cd),
    CONST64(0xb5c0fbcfec4d3b2f), CONST64(0xe9b5dba58189dbbc),
    CONST64(0x3956c25bf348b538), CONST64(0x59f111f1b605d019),
    CONST64(0x923f82a4af194f9b), CONST64(0xab1c5ed5da6d8118),
    CONST64(0xd807aa98a3030242), CONST64(0x12835b0145706fbe),
    CONST64(0x243185be4ee4b28c), CONST64(0x550c7dc3d5ffb4e2),
    CONST64(0x72be5d74f27b896f), CONST64(0x80deb1fe3b1696b1),
    CONST64(0x9bdc06a725c71235), CONST64(0xc19bf174cf692694),
    CONST64(0xe49b69c19ef14ad2), CONST64(0xefbe4786384f25e3),
    CONST64(0x0fc19dc68b8cd5b5), CONST64(0x240ca1cc77ac9c65),
    CONST64(0x2de92c6f592b0275), CONST64(0x4a7484aa6ea6e483),
    CONST64(0x5cb0a9dcbd41fbd4), CONST64(0x76f988da831153b5),
    CONST64(0x983e5152ee66dfab), CONST64(0xa831c66d2db43210),
    CONST64(0xb00327c898fb213f), CONST64(0xbf597fc7beef0ee4),
    CONST64(0xc6e00bf33da88fc2), CONST64(0xd5a79147930aa725),
    CONST64(0x06ca6351e003826f), CONST64(0x142929670a0e6e70),
    CONST64(0x27b70a8546d22ffc), CONST64(0x2e1b21385c26c926),
    CONST64(0x4d2c6dfc5ac42aed), CONST64(0x53380d139d95b3df),
    CONST64(0x650a73548baf63de), CONST64(0x766a0abb3c77b2a8),
    CONST64(0x81c2c92e47edaee6), CONST64(0x92722c851482353b),
    CONST64(0xa2bfe8a14cf10364), CONST64(0xa81a664bbc423001),
    CONST64(0xc24b8b70d0f89791), CONST64(0xc76c51a30654be30),
    CONST64(0xd192e819d6ef5218), CONST64(0xd69906245565a910),
    CONST64(0xf40e35855771202a), CONST64(0x106aa07032bbd1b8),
    CONST64(0x19a4c116b8d2d0c8), CONST64(0x1e376c085141ab53),
    CONST64(0x2748774cdf8eeb99), CONST64(0x34b0bcb5e19b48a8),
    CONST64(0x391c0cb3c5c95a63), CONST64(0x4ed8aa4ae3418acb),
    CONST64(0x5b9cca4f7763e373), CONST64(0x682e6ff3d6b2b8a3),
    CONST64(0x748f82ee5defb2fc), CONST64(0x78a5636f43172f60),
    CONST64(0x84c87814a1f0ab72), CONST64(0x8cc702081a6439ec),
    CONST64(0x90befffa23631e28), CONST64(0xa4506cebde82bde9),
    CONST64(0xbef9a3f7b2c67915), CONST64(0xc67178f2e372532b),
    CONST64(0xca273eceea26619c), CONST64(0xd186b8c721c0c207),
    CONST64(0xeada7dd6cde0eb1e), CONST64(0xf57d4f7fee6ed178),
    CONST64(0x06f067aa72176fba), CONST64(0x0a637dc5a2c898a6),
    CONST64(0x113f9804bef90dae), CONST64(0x1b710b35131c471b),
    CONST64(0x28db77f523047d84), CONST64(0x32caab7b40c72493),
    CONST64(0x3c9ebe0a15c9bebc), CONST64(0x431d67c49c100d4c),
    CONST64(0x4cc5d4becb3e42b6), CONST64(0x597f299cfc657e2a),
    CONST64(0x5fcb6fab3ad6faec), CONST64(0x6c44198c4a475817)
};


/* compress 1024-bits */
int s_sha512_compress(sha512_state *state, const unsigned char *buf)
{
    uint64_t S[8], W[80], t0, t1;
    int i;

    /* copy state into S */
    for (i = 0; i < 8; i++) {
        S[i] = state->state[i];
    }

    /* copy the state into 1024-bits into W[0..15] */
    for (i = 0; i < 16; i++) {
        LOAD64H(W[i], buf + (8*i));
    }

    /* fill W[16..79] */
    for (i = 16; i < 80; i++) {
        W[i] = Gamma1_512(W[i - 2]) + W[i - 7] + Gamma0_512(W[i - 15]) + W[i - 16];
    }

    /* Compress */
#define RND(a,b,c,d,e,f,g,h,i)                    \
     t0 = h + Sigma1_512(e) + Ch(e, f, g) + K_512[i] + W[i];   \
     t1 = Sigma0_512(a) + Maj(a, b, c);                  \
     d += t0;                                        \
     h  = t0 + t1;

    for (i = 0; i < 80; i += 8) {
        RND(S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],i+0);
        RND(S[7],S[0],S[1],S[2],S[3],S[4],S[5],S[6],i+1);
        RND(S[6],S[7],S[0],S[1],S[2],S[3],S[4],S[5],i+2);
        RND(S[5],S[6],S[7],S[0],S[1],S[2],S[3],S[4],i+3);
        RND(S[4],S[5],S[6],S[7],S[0],S[1],S[2],S[3],i+4);
        RND(S[3],S[4],S[5],S[6],S[7],S[0],S[1],S[2],i+5);
        RND(S[2],S[3],S[4],S[5],S[6],S[7],S[0],S[1],i+6);
        RND(S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[0],i+7);
    }
#undef RND

    /* feedback */
    for (i = 0; i < 8; i++) {
        state->state[i] = state->state[i] + S[i];
    }
    return 0;
}


/**
   Process a block of memory though the hash
   @param state  The hash state
   @param in     The data to hash
   @param inlen  The length of the data (bytes)
   @return 0 if successful
*/
int sha512_process(sha512_state *state, const unsigned char *in, unsigned long inlen)
{
    unsigned long n;
    if (state->curlen > sizeof(state->buf)) {
       return 1;
    }
    if (((state->length + inlen * 8) < state->length) || ((inlen * 8) < inlen)) {
      return 1;
    }
    while (inlen > 0) {
        if (state->curlen == 0 && inlen >= 128) {
           if (s_sha512_compress(state, in) != 0) {
              return 1;
           }
           state->length += 128 * 8;
           in             += 128;
           inlen          -= 128;
        } else {
           n = MIN(inlen, (128 - state->curlen));
           memcpy(state->buf + state->curlen, in, (size_t)n);
           state->curlen += n;
           in             += n;
           inlen          -= n;
           if (state->curlen == 128) {
              if (s_sha512_compress(state, state->buf) != 0) {
                 return 1;
              }
              state->length += 8*128;
              state->curlen = 0;
           }
       }
    }
    return 0;
}


/**
   Terminate the hash to get the digest
   @param state     The hash state
   @param out [out] The destination of the hash (64 bytes)
   @return 0 if successful
*/
int sha512_done(sha512_state *state, unsigned char *out)
{
    int i;

    if (state->curlen >= sizeof(state->buf)) {
       return 1;
    }

    /* increase the length of the message */
    state->length += state->curlen * CONST64(8);

    /* append the '1' bit */
    state->buf[state->curlen++] = (unsigned char)0x80;

    /* if the length is currently above 112 bytes we append zeros
     * then compress.  Then we can fall back to padding zeros and length
     * encoding like normal.
     */
    if (state->curlen > 112) {
        while (state->curlen < 128) {
            state->buf[state->curlen++] = (unsigned char)0;
        }
        s_sha512_compress(state, state->buf);
        state->curlen = 0;
    }

    /* pad upto 120 bytes of zeroes
     * note: that from 112 to 120 is the 64 MSB of the length.  We assume that you won't hash
     * > 2^64 bits of data... :-)
     */
    while (state->curlen < 120) {
        state->buf[state->curlen++] = (unsigned char)0;
    }

    /* store length */
    STORE64H(state->length, state->buf+120);
    s_sha512_compress(state, state->buf);

    /* copy output */
    for (i = 0; i < 8; i++) {
        STORE64H(state->state[i], out+(8*i));
    }

    return 0;
}


/********************************* SHA-384 ************************************************/
/**
*  Initialize the hash state
*  @param state   The hash state you wish to initialize
*  @return 0 if successful
*/
int sha384_init(sha512_state *state)
{

    state->curlen = 0;
    state->length = 0;
    state->state[0] = CONST64(0xcbbb9d5dc1059ed8);
    state->state[1] = CONST64(0x629a292a367cd507);
    state->state[2] = CONST64(0x9159015a3070dd17);
    state->state[3] = CONST64(0x152fecd8f70e5939);
    state->state[4] = CONST64(0x67332667ffc00b31);
    state->state[5] = CONST64(0x8eb44a8768581511);
    state->state[6] = CONST64(0xdb0c2e0d64f98fa7);
    state->state[7] = CONST64(0x47b5481dbefa4fa4);
    return 0;
}


/*
   Terminate the hash to get the digest
   @param state     The hash state
   @param out [out] The destination of the hash (48 bytes)
   @return 0 if successful
*/
int sha384_done(sha512_state *state, unsigned char *out)
{
   unsigned char buf[64];

    if (state->curlen >= sizeof(state->buf)) {
       return 1;
    }

   sha512_done(state, buf);
   memcpy(out, buf, 48);
   return 0;
}


/********************** SHA2 wrapper functions *******************************************/
void sha2_256(unsigned char *out, const unsigned char *in, const uint64_t inlen) 
{
    sha256_state state;

    sha256_init(&state);
    sha256_process(&state, in, inlen);
    sha256_done(&state, out);
}

void sha2_384(unsigned char *out, const unsigned char *in, const uint64_t inlen) 
{
    sha512_state state;

    sha384_init(&state);
    sha512_process(&state, in, inlen);
    sha384_done(&state, out);
}

void sha2_512(unsigned char *out, const unsigned char *in, const uint64_t inlen) 
{
    sha512_state state;

    sha512_init(&state);
    sha512_process(&state, in, inlen);
    sha512_done(&state, out);
}
