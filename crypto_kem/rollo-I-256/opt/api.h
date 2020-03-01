#ifndef API_H
#define API_H

/* ROLLO-I-256*/

#define  PARAM_D 						8
#define  PARAM_R 						7
#define  PARAM_M 						113
#define  PARAM_M_LEN_WORD				4
#define	 PARAM_M_LEN_BYTES 				15
#define  PARAM_M_LEN_BYTES_PADDING 		16
#define  PARAM_N 						67
#define NB_COEFF_MODULO_GF2_N 			4
#define MODULO_GF2M_LEN       			4
#define SHARED_SECRET_BYTES 			64
#define NM_WORDS                       PARAM_N *PARAM_M_LEN_WORD
#define NM_WORDS_PAD                  (PARAM_N + 6)*PARAM_M_LEN_WORD  //With Karatsuba padding

#define CRYPTO_SECRETKEYBYTES 40
#define CRYPTO_PUBLICKEYBYTES 947
#define CRYPTO_BYTES 64
#define CRYPTO_CIPHERTEXTBYTES 947

int crypto_kem_keypair(unsigned char* pk, unsigned char* sk);
int crypto_kem_enc(unsigned char* ct, unsigned char* ss, unsigned char* pk);
int crypto_kem_dec(unsigned char* ss, unsigned char* ct, unsigned char* sk);

#endif
