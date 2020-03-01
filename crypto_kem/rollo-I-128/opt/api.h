#ifndef API_H
#define API_H

/* ROLLO-I-128*/

#define  PARAM_D 						6
#define  PARAM_R 						5
#define  PARAM_M 						79
#define  PARAM_M_LEN_WORD				3
#define	 PARAM_M_LEN_BYTES 				10
#define  PARAM_M_LEN_BYTES_PADDING 		12
#define  PARAM_N 						47
#define NB_COEFF_MODULO_GF2_N 			2
#define MODULO_GF2M_LEN       			3
#define SHARED_SECRET_BYTES 			64
#define NM_WORDS                      PARAM_N*PARAM_M_LEN_WORD
#define NM_WORDS_PAD                      (PARAM_N + 1)*PARAM_M_LEN_WORD   //with Karatsuba padding

#define CRYPTO_SECRETKEYBYTES 40
#define CRYPTO_PUBLICKEYBYTES 465
#define CRYPTO_BYTES 64
#define CRYPTO_CIPHERTEXTBYTES 465

int crypto_kem_keypair(unsigned char* pk, unsigned char* sk);
int crypto_kem_enc(unsigned char* ct, unsigned char* ss, unsigned char* pk);
int crypto_kem_dec(unsigned char* ss, unsigned char* ct, unsigned char* sk);

#endif
