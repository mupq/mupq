#include <string.h>
#include "bin-lwe.h"
#include "rand.h"
#include "lac_param.h"

//generate the public parameter a from seed
int gen_a(uint8_t *a,  const uint8_t *seed)
{
	pseudo_random_bytes(a,DIM_N,seed);
	
	return 0;
}
 
//generate the small random vector for secret and error, with fixed hamming weight
//use for e,e1,e2
int gen_e(uint8_t *e,  uint8_t *seed)
{
	int i;
	uint16_t buf[NUM_ONE*2];
	gen_r((uint8_t *)buf,seed);
	memset(e,0,DIM_N);
	for(i=0;i<NUM_ONE;i++)
	{
		e[buf[i]]=1;
		e[buf[i+NUM_ONE]]=Q-1;
	}
	
	return 0;
}

//for r,s
int gen_r(uint8_t *r,  uint8_t *seed)
{
    int i,p;
    uint16_t tmp;
    uint16_t  r_buf[DIM_N],index[SAMPLE_LEN],tmp_index,index_mk;
    uint16_t mk=DIM_N-1;
    unsigned int mask_p,loop=SAMPLE_LEN;

    //init r to be 1,2,3,4,5
    for(i=0;i<DIM_N;i++)
    {
        r_buf[i]=i;
    }
    
    
    p=0;
    while(p<NUM_ONE*2)
    {
        pseudo_random_bytes((uint8_t *)index,SAMPLE_LEN*2,seed);
        //shuffle
        for(i=0;i<loop;i++)
        {
            //check index
            index_mk=index[i]&mk;
            tmp_index=index_mk>=p ? index_mk: p;
            mask_p=index_mk>=p ? 1:0;
            //swap
            tmp=r_buf[tmp_index];
            r_buf[tmp_index]=r_buf[p];
            r_buf[p]=tmp;
            //move to the next
            p+=mask_p;
        }
        //update seed
        if(p<NUM_ONE*2)
        {
            memcpy(seed,(uint8_t *)index,SEED_LEN);
        }
    }
    //compy the first NUM_ONE positions to r
    memcpy(r,r_buf,NUM_ONE*2*sizeof(uint16_t));
    
    return 0;
}

int mul_core(const uint8_t *a, const uint8_t *s, int32_t *sum1, int32_t *sum2, unsigned int vec_num)
{
	int i,j;
	uint8_t v[DIM_N+DIM_N];
	uint8_t *v1_p,*v2_p;

	uint16_t *s_one,*s_minusone;
	
	
	//construct matrix of a
	for(i=0;i<DIM_N;i++)
	{
		v[i]=Q-a[i];
	}	
	memcpy(v+DIM_N,a,DIM_N);
	
	s_one=(uint16_t*)s;
	s_minusone=(uint16_t*)s+NUM_ONE;
	memset(sum1,0,vec_num*sizeof(int32_t));
	memset(sum2,0,vec_num*sizeof(int32_t));
	
	for(i=0;i<NUM_ONE;i++)
	{
		v1_p=v+DIM_N-s_one[i];
		for(j=0;j<vec_num;j++)
		{
			sum1[j]+=v1_p[j];
		}
		
		v2_p=v+DIM_N-s_minusone[i];
		for(j=0;j<vec_num;j++)
		{
			sum2[j]-=v2_p[j];
		}
	}

	return 0;
}


// poly_mul  b=[as]
int poly_mul(const uint8_t *a, const uint8_t *s, uint8_t *b, unsigned int vec_num)
{
	int i;
	int32_t sum1[DIM_N],sum2[DIM_N];
	
	mul_core(a,s,sum1,sum2,vec_num);
	
	for(i=0;i<vec_num;i++)
	{
		b[i]=(sum1[i]+sum2[i]);
	}
	return 0;
}
//b=as+e
int poly_aff(const uint8_t *a, const uint8_t *s, uint8_t *e, uint8_t *b, unsigned int vec_num)
{
	int i;
	int32_t sum1[DIM_N],sum2[DIM_N];
	
	mul_core(a,s,sum1,sum2,vec_num);
	
	for(i=0;i<vec_num;i++)
	{
		b[i]=(sum1[i]+sum2[i]+e[i]);
	}
	return 0;
}
 
// compress: cut the low 4bit
int poly_compress(const uint8_t *in,  uint8_t *out, const unsigned int vec_num)
{
	int i,loop;
	loop=vec_num/2;
	for(i=0;i<loop;i++)
	{
		out[i]=(in[i*2])>>4;
		out[i]=out[i]^(in[i*2+1]&0xf0);
	}
	
	return 0;
}
// decompress: set the low 4bits to be 0
int poly_decompress(const uint8_t *in,  uint8_t *out, const unsigned int vec_num)
{
	int i,loop;
	loop=vec_num/2;
	for(i=0;i<loop;i++)
	{
		out[i*2]=(in[i]<<4)^0x08;
		out[i*2+1]=(in[i]&0xf0)^0x08;
	}
	
	return 0;
}


