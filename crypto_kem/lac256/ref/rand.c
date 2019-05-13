#include <string.h>
#include "aes.h"
#include "sha2.h"
#include "lac_param.h"
#include "rand.h"
#include "randombytes.h"

//random bytes
int random_bytes(unsigned char *r, unsigned int len)
{
	//check parameter
	if(r==NULL)
	{
		return 1;
	}
	// call the random function
//	RAND_bytes(r,len);
	randombytes(r,len);
	return 0;
}

//pseudo-random bytes
int pseudo_random_bytes(unsigned char *r, unsigned int len, const unsigned char *seed)
{
    /*
	int c_len;
	unsigned char data[16],c[16];
	unsigned int *p=(unsigned int *)data;
    int i,loop=len/16;
    */
    unsigned char data[12] = {0};
	//check  parameter
	if(r==NULL || seed==NULL)
	{
		return 1;
	}
    /*
	memset(r,0,len);
	EVP_CIPHER_CTX *ctx;
	ctx = EVP_CIPHER_CTX_new();
	EVP_EncryptInit_ex(ctx, EVP_aes_256_ctr(), NULL, seed, NULL);
	memset(data,0,AES_BLOCK_SIZE);
	for(i=0;i<loop;i++)
	{
	//	*p=i;//set counter
		EVP_EncryptUpdate(ctx, r+i*AES_BLOCK_SIZE, &c_len, data, AES_BLOCK_SIZE);
	}
	//check tail
	if(len%AES_BLOCK_SIZE>0)
	{
	//	*p=loop;
		EVP_EncryptUpdate(ctx, c, &c_len, data, AES_BLOCK_SIZE);
	}
	memcpy(r+loop-1,c,len%AES_BLOCK_SIZE);
	EVP_CIPHER_CTX_free(ctx);
    */
    aes256ctx ctx;
    aes256_keyexp(&ctx, seed);
    aes256_ctr(r, len, data, &ctx);
	return 0;
}

//hash
int hash(const unsigned char *in, unsigned int len_in, unsigned char * out)
{
	//check  parameter
	if(in==NULL || out==NULL)
	{
		return 1;
	}

	#if defined LAC128
	sha256(out,in,len_in);
	#endif

	#if defined LAC192
	sha256(out,in,len_in);
	#endif

	#if defined LAC256
	sha256(out,in,len_in);
	#endif


	return 0;
}

//generate seed
int gen_seed(unsigned char *in, unsigned int len_in, unsigned char * out)
{
	//check  parameter
	if(in==NULL || out==NULL)
	{
		return 1;
	}

	sha256(out,in,len_in);
	return 0;
}

