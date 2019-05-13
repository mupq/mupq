#include <string.h>
#include "bin-lwe.h"
#include "rand.h"
#include "lac_param.h"

//generate the public parameter a from seed
int gen_a(unsigned char *a,  const unsigned char *seed)
{
	int i,j;
	unsigned char buf[MESSAGE_LEN];
	//check the pointers
	if(a==NULL || seed==NULL)
	{
		return 1;
	}

	pseudo_random_bytes(a,DIM_N,seed);

	hash(seed,SEED_LEN,buf);
	j=0;
	for(i=0;i<DIM_N;i++)
	{
		while(a[i]>=Q)
		{
			memcpy(a+i,buf+(j++),1);//replace a[i] with buf[j]
			if(j>=MESSAGE_LEN)
			{
				hash(buf,MESSAGE_LEN,buf);//use hash chain to refresh buf
				j=0;
			}
		}
	}

	return 0;
}

//generate the small random vector for secret and error, with fixed hamming weight
int gen_psi_fix_ham(signed char *e, unsigned int vec_num, unsigned char *seed)
{
	if(e==NULL)
	{
		return 1;
	}

	unsigned char buf[MESSAGE_LEN];
	int i,bound1,bound2,j;
	uint16_t *p_index,mask;


	#if defined PSI_SQUARE

	//Pr[e[i]=-1]=1/8,Pr[s[i]=0]=3/4,Pr[e[i]=1]=1/8,
	unsigned char r[vec_num/2];
	//generate vec_num/2 bytes, use 4 bits to generate one error item
	pseudo_random_bytes(r,vec_num/2,seed);
	//set the number  of 1 or -1
	bound1=vec_num/8;
	bound2=vec_num/4;

	#else

	//Pr[e[i]=-1]=1/4,Pr[s[i]=0]=1/2,Pr[e[i]=1]=1/4,
	unsigned char r[vec_num];
	//generate vec_num/2 bytes, use 4 bits to generate one error item
	pseudo_random_bytes(r,vec_num,seed);
	//set the number  of 1 or -1
	bound1=vec_num/4;
	bound2=vec_num/2;

	#endif

	//init e to be 0
	memset(e,0,vec_num);
	//set mask
	mask=DIM_N-1;
	//if collision, refresh index
	hash(seed,SEED_LEN,buf);
	//set index pointer
	p_index=(uint16_t*)r;
	j=0;
	//set 1
	for(i=0;i<bound1;i++)
	{
		while(e[p_index[i]&mask])
		{
			//refresh index
			memcpy(p_index+i,buf+j,2);//replace p_index[i] with buf[j]
			j+=2;
			if(j>=MESSAGE_LEN)
			{
				hash(buf,MESSAGE_LEN,buf);//use hash chain to refresh buf
				j=0;
			}
		}
		e[p_index[i]&mask]=1;
	}
	//set -1
	for(i=bound1;i<bound2;i++)
	{
		while(e[p_index[i]&mask])
		{
			//refresh index
			memcpy(p_index+i,buf+j,2);//replace p_index[i] with buf[j]
			j+=2;
			if(j>=MESSAGE_LEN)
			{
				hash(buf,MESSAGE_LEN,buf);//use hash chain to refresh buf
				j=0;
			}
		}
		e[p_index[i]&mask]=-1;
	}

	return 0;
}

//generate the small random vector for secret and error
int gen_psi_std(signed char *e, unsigned int vec_num, unsigned char *seed)
{
	unsigned int i;

	if(e==NULL)
	{
		return 1;
	}

	#if defined PSI_SQUARE
	//Pr[e[i]=-1]=1/8,Pr[s[i]=0]=3/4,Pr[e[i]=1]=1/8,
	unsigned char r[vec_num/2],*p1,*p2,*p3;
	int e_0,e_1;
	//generate vec_num/2 bytes, use 4 bits to generate one error item
	pseudo_random_bytes(r,vec_num/2,seed);
	//COMPUTE e from r
	p1=r+vec_num/8;
	p2=p1+vec_num/8;
	p3=p2+vec_num/8;
	for(i=0;i<vec_num;i++)
	{
		e_0=((r[i/8]>>(i%8))&1)-((p1[i/8]>>(i%8))&1);
		e_1=((p2[i/8]>>(i%8))&1)-((p3[i/8]>>(i%8))&1);
		e[i]=e_0*e_1;
	}

	#else
	//Pr[e[i]=-1]=1/4,Pr[s[i]=0]=1/2,Pr[e[i]=1]=1/4,
	unsigned char r[vec_num/4],*p;

	//generate vec_num/4 bytes, use two bits to generate one error item
	pseudo_random_bytes(r,vec_num/4,seed);
	//COMPUTE e from r
	p=r+vec_num/8;
	for(i=0;i<vec_num;i++)
	{
		e[i]=((r[i/8]>>(i%8))&1)-((p[i/8]>>(i%8))&1);
	}

	#endif

	return 0;
}

// poly_mul  b=[as]
int poly_mul(const unsigned char *a, const signed char *s, unsigned char *b, unsigned int vec_num)
{
	unsigned int i,j,loop=DIM_N;
	unsigned char v[DIM_N+DIM_N];
	unsigned char *v_p;
	int32_t sum;

	//construct matrix of a
	for(i=0;i<DIM_N;i++)
	{
		v[i]=a[DIM_N-1-i];
		v[i+DIM_N]=Q-v[i];
	}

	for(i=0;i<vec_num;i++)
	{
		sum=0;
		v_p=(v+DIM_N-i-1);

		for(j=0;j<loop;j++)
		{
			sum+=v_p[j]*s[j];
		}
		b[i]=(sum+BIG_Q)%Q;
	}

	return 0;
}
//b=as+e
int poly_aff(const unsigned char *a, const signed char *s, signed char *e, unsigned char *b, unsigned int vec_num)
{
	unsigned int i,j,loop=DIM_N;
	unsigned char v[DIM_N+DIM_N];
	unsigned char *v_p;
	int32_t sum;

	//construct matrix of a
	for(i=0;i<DIM_N;i++)
	{
		v[i]=a[DIM_N-1-i];
		v[i+DIM_N]=Q-v[i];
	}
	for(i=0;i<vec_num;i++)
	{
		v_p=(v+DIM_N-i-1);
		sum=0;
		for(j=0;j<loop;j++)
		{
			sum+=v_p[j]*s[j];
		}
		b[i]=(sum+e[i]+BIG_Q)%Q;
	}

	return 0;
}

// compress: cut the low 4bit
int poly_compress(const unsigned char *in,  unsigned char *out, const unsigned int vec_num)
{
	int i,loop;
	loop=vec_num/2;
	for(i=0;i<loop;i++)
	{
		out[i]=(in[i*2]+0x08)>>4;
		out[i]=out[i]^((in[i*2+1]+0x08)&0xf0);
	}

	return 0;
}
// decompress: set the low 4bits to be 0
int poly_decompress(const unsigned char *in,  unsigned char *out, const unsigned int vec_num)
{
	int i,loop;
	loop=vec_num/2;
	for(i=0;i<loop;i++)
	{
		out[i*2]=in[i]<<4;
		out[i*2+1]=in[i]&0xf0;
	}

	return 0;
}

