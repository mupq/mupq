#include <string.h>
#include "api.h"
#include "rand.h"
#include "bch.h"
#include "ecc.h"
#include "bin-lwe.h"

#define RATIO 125

//key generation
int crypto_encrypt_keypair( unsigned char *pk, unsigned char *sk)
{
	//check parameter
	if(pk==NULL || sk==NULL)
	{
		return -1;
	}
	kg(pk,sk);

	return 0;
}

//encryption
int crypto_encrypt( unsigned char *c, unsigned long long *clen, const unsigned char *m, unsigned long long mlen, const unsigned char *pk)
{
	//check parameter
	if(c==NULL || m==NULL || pk==NULL)
	{
		return -1;
	}
	if(mlen>MESSAGE_LEN)
	{
		return -1;
	}

	//call pke encryption function
	pke_enc(pk,m, mlen,c,clen);

	return 0;
}
//decryption
int crypto_encrypt_open(unsigned char *m, unsigned long long *mlen,const unsigned char *c, unsigned long long clen,const unsigned char *sk)
{
	//check parameter
	if(sk==NULL || m==NULL || c==NULL || mlen==NULL)
	{
		return -1;
	}

	//call pke decryption function
	pke_dec(sk,c,clen,m,mlen);

	return 0;
}

//key generation with seed
int kg_seed(unsigned char *pk, unsigned char *sk, unsigned char *seed)
{
	unsigned char seeds[3*SEED_LEN];
	unsigned char a[DIM_N];
	signed char e[DIM_N];
	//check pointer
	if(pk==NULL || sk==NULL)
	{
		return -1;
	}
	//generate three seeds for a,sk,e
	pseudo_random_bytes(seeds,3*SEED_LEN,seed);
	//copy the seed of a to pk
	memcpy(pk,seeds,SEED_LEN);
	//generate a
	gen_a(a,pk);
	//generate  sk,e
	gen_psi_fix_ham((signed char*)sk,DIM_N,seeds+SEED_LEN);
	gen_psi_fix_ham((signed char*)e,DIM_N,seeds+2*SEED_LEN);
	//compute pk=a*sk+e
	poly_aff(a,(signed char *)sk,e,pk+SEED_LEN,DIM_N);
	//copy pk=as+e to the second part of sk, now sk=s|pk
	memcpy(sk+DIM_N,pk,PK_LEN);
	return 0;
}

//key generation
int kg(unsigned char *pk, unsigned char *sk)
{
	unsigned char seed[SEED_LEN];

	//generate seed
	random_bytes(seed,SEED_LEN);
	//key generation with seed
	kg_seed(pk,sk,seed);

	return 0;
}
// encryption
int pke_enc(const unsigned char *pk, const unsigned char *m, unsigned long long mlen, unsigned char *c, unsigned long long *clen)
{
	unsigned char seed[SEED_LEN];

	//generate seed
	random_bytes(seed,SEED_LEN);
	//encrypt with seed
	pke_enc_seed(pk,m,mlen,c,clen,seed);

	return 0;
}
// decrypt
int pke_dec(const unsigned char *sk, const unsigned char *c,unsigned long long clen, unsigned char *m, unsigned long long *mlen)
{
	unsigned char out[DIM_N];
	unsigned char code[CODE_LEN],*p_code;
	unsigned char c2[C2_VEC_NUM];
	int i,c2_len=(clen-DIM_N)*2;
	unsigned char m_buf[MESSAGE_LEN];

	//check parameter
	if(sk==NULL || m==NULL || c==NULL)
	{
		return -1;
	}

	#ifdef LAC256 //D2 decoding

	int vec_bound=c2_len/2;
	int temp1,temp2;
	int center_point=Q/2;
	int d2_bound=Q/2;
	//compute mlen
	*mlen=c2_len/16-ECC_LEN;
	//shif the pointer of ecc data
	p_code=code+(DATA_LEN-(*mlen));
	//init code
	memset(code,0,CODE_LEN);
	//c2 decompress
	poly_decompress(c+DIM_N,c2,c2_len);
	//c1*sk
	poly_mul(c,(signed char *)sk,out,DIM_N);
	//compute c2-c1*sk and recover data from m*q/2+e
	for(i=0;i<vec_bound;i++)
	{
		//D2 decoding:compute m*q/2+e1 + m*q/2+e2 in [0,2*Q]
		temp1=(c2[i]-out[i]+Q)%Q;
		temp2=(c2[i+vec_bound]-out[i+vec_bound]+Q)%Q;

		//shift
		if(temp1<center_point)
		{
			temp1=center_point-temp1+center_point;//mirror around Q/2
		}
		if(temp2<center_point)
		{
			temp2=center_point-temp2+center_point;//mirror around Q/2
		}
		//merge erors
		temp1+=temp2-Q;

		//recover m from m*q/2+e1 + m*q/2+e2, RATIO=q/2
		if(temp1<d2_bound)
		{
			p_code[i/8]=p_code[i/8]^(1<<(i%8));
		}
	}
	//bch decode to recover m
	ecc_dec(m_buf,code);
	//get plaintext
	memcpy(m,m_buf+(DATA_LEN-(*mlen)),*mlen);

	#else


	int temp;
	int low=Q/4;
	int high=Q*3/4;
	//c2 decompress
	poly_decompress(c+DIM_N,c2,c2_len);
	//c1*sk
	poly_mul(c,(signed char*)sk,out,DIM_N);
	//compute mlen
	*mlen=c2_len/8-ECC_LEN;
	//init code
	memset(code,0,CODE_LEN);
	//shift the pointer of code
	p_code=code+(DATA_LEN-(*mlen));

	for(i=0;i<c2_len;i++)
	{
		//compute m*q/2+e in [0,Q]
		temp=(c2[i]-out[i]+Q)%Q;

		//recover m from m*q/2+e, RATIO=q/2
		if(temp>=low && temp<high)
		{
			p_code[i/8]=p_code[i/8]^(1<<(i%8));
		}
	}

	//bch decode to recover m
	ecc_dec(m_buf,code);
	//get plaintext
	memcpy(m,m_buf+(DATA_LEN-(*mlen)),*mlen);
	#endif

	return 0;
}

// encryption with seed
int pke_enc_seed(const unsigned char *pk, const unsigned char *m, unsigned long long mlen, unsigned char *c, unsigned long long *clen, unsigned char *seed)
{
	unsigned char code[CODE_LEN],seeds[3*SEED_LEN],*p_code;
	signed char r[DIM_N];
	signed char e1[DIM_N],e2[C2_VEC_NUM];
	unsigned char c2[C2_VEC_NUM];
	unsigned char a[DIM_N];
	unsigned char m_buf[MESSAGE_LEN];
	int i,c2_len;

	//check parameter
	if(pk==NULL || m==NULL || c==NULL )
	{
		return -1;
	}
	if(mlen>MESSAGE_LEN)
	{
		return -1;
	}

	//generate  a from seed in the first part of pk
	gen_a(a,pk);
	//package m_buf
	memset(m_buf,0,MESSAGE_LEN);
	//set data
	memcpy(m_buf+(DATA_LEN-mlen),m,mlen);
	//encode m with ecc code
	ecc_enc(m_buf,code);
	//set p_code
	p_code=code+(DATA_LEN-mlen);
	//generate three seeds for r,e1,e2
	pseudo_random_bytes(seeds,3*SEED_LEN,seed);
	//generate random vector r
	gen_psi_fix_ham(r,DIM_N,seeds);
	//generate error vector e1
	gen_psi_fix_ham(e1,DIM_N,seeds+SEED_LEN);
	//compute c1=a*r+e1
	poly_aff(a,r,e1,c,DIM_N);

	//D2 encoding
	#ifdef LAC256

	//compute the length of c2
	c2_len=(mlen+ECC_LEN)*8*2;
	//generate error vector e2
	gen_psi_std(e2,c2_len,seeds+2*SEED_LEN);

	int vec_bound=c2_len/2;
	int8_t message;
	//compute  code*q/2+e2,
	for(i=0;i<vec_bound;i++)
	{
		//RATIO=q/2. add code*q/2 to e2
		message=RATIO*((p_code[i/8]>>(i%8))&1);
		e2[i]=e2[i]+message;
		//D2 encode, repeat at i+vec_bound
		e2[i+vec_bound]=e2[i+vec_bound]+message;
	}

	#else

	//compute the length of c2
	c2_len=(mlen+ECC_LEN)*8;
	//generate error vector e2
	gen_psi_std(e2,c2_len,seeds+2*SEED_LEN);
	//compute  code*q/2+e2,
	for(i=0;i<c2_len;i++)
	{
		//RATIO=q/2. add code*q/2 to e2
		e2[i]=e2[i]+RATIO*((p_code[i/8]>>(i%8))&1);
	}
	#endif
	//c2=b*r+e2+m*[q/2]
	poly_aff(pk+SEED_LEN,r,e2,c2,c2_len);
	//compress c2
	poly_compress(c2,c+DIM_N,c2_len);
	*clen=DIM_N+c2_len/2;

	return 0;

}

