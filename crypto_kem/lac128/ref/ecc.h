#include "lac_param.h"

#if defined(LAC128)
//bch(511,264,33)
#define DATA_LEN 32//MESSAGE_LEN+one byte for message length
#define DATABUF_LEN 48//CODE_LEN-ECC_LEN
#define ECCBUF_LEN 64
#define ECC_LEN 18
#define MAX_ERROR 16
#define CODE_LEN 64
#define LOG_CODE_LEN 9
#endif

#if defined(LAC192)
//bch(511,264,33)
#define DATA_LEN 32
#define DATABUF_LEN 48
#define ECCBUF_LEN 64
#define ECC_LEN 9
#define MAX_ERROR 8
#define CODE_LEN 64
#define LOG_CODE_LEN 9
#endif

#if defined(LAC256)
//D2+bch(511,264,33)
#define DATA_LEN 32
#define DATABUF_LEN 48
#define ECCBUF_LEN 64
#define ECC_LEN 18
#define MAX_ERROR 16
#define CODE_LEN 64 
#define LOG_CODE_LEN 9
#endif

//error correction encode
int ecc_enc(const unsigned char *d, unsigned char *c);

//error corrction decode
int ecc_dec(unsigned char *d, const unsigned char *c);

