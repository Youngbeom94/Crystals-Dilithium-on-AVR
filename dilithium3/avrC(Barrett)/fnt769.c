#include <stdint.h>
#include <stdlib.h>
#include "fnt769.h"

__flash
int16_t FNT_Q2_zetas[128] = {0, -164, -81, 361, -186, 3, 250, 120, -129, -308, 223, -16, -143, 362, -337, -131, -75, -36, 76, 98, 203, 282, -339, -255, 178, 270, 199, 34, -369, 192, -149, -10, -80, -346, -124, 2, 114, 147, -54, -272, -169, 288, 161, -15, -86, 51, -364, -267, 170, -226, -121, 188, -50, -24, 307, -191, 263, 157, -246, 128, 375, 180, -380, 279, -341, -379, 202, 220, 236, 21, 212, 71, -134, 151, 23, -112, -232, 227, -52, -148, 244, -252, -237, -83, -117, -333, -66, -247, -292, 352, -145, 238, -276, -194, -274, -70, 209, -115, -99, 14, 29, 260, -378, -366, 355, -291, 358, -105, 167, 357, -241, -331, -348, -44, -78, -222, -350, -168, -158, 201, 303, 330, -184, 127, 318, -278, -353, -354};

static int ct_lt_s32(uint32_t x, uint32_t y) {return (x ^ ((x ^ (x - y)) & (y ^ (x - y)))) >> 31;}

static int ct_gt_s32(uint32_t x, uint32_t y) {return ct_lt_s32(y, x);}

int32_t fnt_freeze(int32_t src)
{
	int32_t r  = src;	
	uint16_t mask;

	mask = 0 -ct_gt_s32(src, FNT_Q2HF); // if (src > FNT_Q2HF) r -= 769;
	r -= mask & 769;
	mask = 0 - ct_lt_s32(src, (uint32_t) FNT_Q2HFNEG); // if (src < -FNT_Q2HF) r += 769;
	r += mask & 769;
	return r;
}

static int16_t fnt769_montgomery_reduce(int32_t a) 
{
	int16_t t;
	int16_t hi = a >> 16;

	t = (int16_t)a * QinvR;
	t = ((int32_t)t * FNT_Q2) >> 16;	
	return hi - t;
}

//static int16_t fnt769_montgomery_reduce_origin(int32_t a) 
//{
//    int16_t t;
//
//	t = (int16_t)a * QinvR;
//	t = (a - (int32_t)t * FNT_Q2) >> 16;
//	return t;
//}

void fnt_ntt(int32_t r[256])
{
	unsigned int len, start, j, k;
	int16_t t, zeta;

	k = 1;
	for (len = 128; len >= 2; len >>= 1) {
		for (start = 0; start < 256; start = j + len) {
			zeta = FNT_Q2_zetas[k++];
			for (j = start; j < start + len; j++)
			{			 
				t = fnt769_montgomery_reduce((int32_t)zeta * r[j + len]);
				r[j + len] = (r[j] - t);
				r[j] = (r[j] + t);
			}
		}
	}
}

void fnt_invntt(int32_t r[256])
{
	unsigned int start, len, j, k;
	int16_t t, zeta;
	const int16_t f = 655; // (1 << 16)^2 / 128

	k = 127;
	for (len = 2; len <= 128; len <<= 1) {
		for (start = 0; start < 256; start = j + len)
		{
			//if (len == 128) zeta = FNT_Q2_zetas[0];
			//else zeta = FNT_Q2_zetas[k--];
			zeta = FNT_Q2_zetas[k--];
			for (j = start; j < start + len; j++)
			{
				t = r[j];
				r[j] = (t + r[j + len]);
				r[j + len] = (r[j + len] - t);
				r[j + len] = fnt769_montgomery_reduce((int32_t)zeta * r[j + len]);
			}
		}
	}

	for (int cnt_i = 0; cnt_i < 256; cnt_i++)
	{
		r[cnt_i] = fnt769_montgomery_reduce((int32_t)r[cnt_i] * f);
	}
}


void fnt769_basemul(int32_t r[2], const int32_t a[2], const int32_t b[2], int16_t zeta)
{
    int32_t tmp[2] = {0};
  
    tmp[0] =  (int32_t)fnt769_montgomery_reduce(a[1]  * b[1]);
    tmp[0] =  (int32_t)fnt769_montgomery_reduce(tmp[0]* zeta);
    tmp[0] += (int32_t)fnt769_montgomery_reduce(a[0]  * b[0]);
    tmp[1] =  (int32_t)fnt769_montgomery_reduce(a[0]  * b[1]);
    tmp[1] += (int32_t)fnt769_montgomery_reduce(a[1]  * b[0]);
    
    r[0] = tmp[0]; r[1] = tmp[1];
}

void point_mul_small(int32_t des[N], int32_t src1[N], int32_t src2[N])
{
  for(int cnt_i = 0 ; cnt_i < N / 4 ; cnt_i++)
  {
    fnt769_basemul(&des[4 * cnt_i    ], &src1[4 * cnt_i    ], &src2[4 * cnt_i    ], FNT_Q2_zetas[64 + cnt_i]);
    fnt769_basemul(&des[4 * cnt_i + 2], &src1[4 * cnt_i + 2], &src2[4 * cnt_i + 2], -FNT_Q2_zetas[64 + cnt_i]);
  } 
}