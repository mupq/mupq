#include "aes_ctr.h"
#include "aes-publicinputs.h"


int AES_128_CTR(unsigned char *output, size_t outputByteLen,
                const unsigned char *input, size_t inputByteLen){
    (void) inputByteLen;
    aes128ctx_publicinputs ctx;
    unsigned char iv[12] = {0};
    aes128_ctr_keyexp_publicinputs(&ctx, input);
    aes128_ctr_publicinputs(output, outputByteLen, iv, &ctx);
    return outputByteLen;
}