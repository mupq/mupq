/*
 * Reference ISO-C11 Implementation of CROSS.
 *
 * @version 1.1 (May 2023)
 *
 * @author: Patrick Karl <patrick.karl@tum.de>
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

#pragma once

#include <stdint.h>
#include <string.h>


void sha2_256(unsigned char *out, const unsigned char *in, const uint64_t inlen);

void sha2_384(unsigned char *out, const unsigned char *in, const uint64_t inlen);

void sha2_512(unsigned char *out, const unsigned char *in, const uint64_t inlen);


typedef struct sha256_state {
    uint64_t length;
    uint32_t curlen;
    uint32_t state[8];
    unsigned char buf[64];
} sha256_state;


typedef struct sha512_state {
    uint64_t length;
    uint64_t curlen;
    uint64_t state[8];
    unsigned char buf[128];
} sha512_state;
