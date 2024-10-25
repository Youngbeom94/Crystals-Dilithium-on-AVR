#include <stdint.h>
#include <stdlib.h>
#include "fnt257.h"

__flash
int16_t fnt_q1_zetas[128] = { -32, 16, 4,  64, 2,  32,  8,  128, -60,  68,  17,  15, -120, -121,  34,  30, -35, -46, 117,  73, -70, -92, -23, -111,  44, -67, -81, -11,  88, 123,  95, -22, -42, 99,  89, -118, -84, -59, -79,  21, -50, -29,  57, -116, -100, -58, 114,  25, -72, -124, -31,  18, 113,   9, -62,  36, -49, -13,  61, -52, -98, -26, 122, -104, -27,  82, -108,  71, -54, -93,  41, -115,  78, -37,  55, 109, -101, -74, 110, -39, -83, -43, -75,  85,  91, -86, 107, -87,  97,  10, -126,  40, -63,  20,   5,  80, 106, -103, -90, 102, -45,  51,  77, -53,  65,  12,   3,  48, -127,  24,   6,  96, -112,   7,  66,  28,  33,  14, -125,  56,  38,  94, -105, 119,  76, -69,  47, -19 };

static int ct_lt_s32(uint32_t x, uint32_t y) {return (x ^ ((x ^ (x - y)) & (y ^ (x - y)))) >> 31;}

static int ct_gt_s32(uint32_t x, uint32_t y) {return ct_lt_s32(y, x);}

int16_t fnt_reduce(int32_t src)
{
	int32_t r;
	uint16_t hi;
	hi = src & 0xFFFF;
	r = hi + (src >> 16);

	uint8_t lo;
	lo = r & 0xFF;
	r = lo - (r >> 8);
	return (int16_t)r;
}

int32_t fnt_freeze(int32_t src)
{
	int32_t r  = src;	
	uint16_t mask;

	mask = 0 -ct_gt_s32(src, FNT_Q1HF); // if (src > FNT_Q1HF) r -= 257;
	r -= mask & 257;
	mask = 0 - ct_lt_s32(src, (uint32_t) FNT_Q1HFNEG); // if (src < -FNT_Q1HF) r += 257;
	r += mask & 257;
	return r;
}

void fnt_ntt(int32_t r[256])
{
	unsigned int len, start, j, k;
	int16_t t, zeta;

	k = 1;
	for (len = 128; len >= 2; len >>= 1) {
		for (start = 0; start < 256; start = j + len) {
			zeta = fnt_q1_zetas[k++];
			for (j = start; j < start + len; j++)
			{
				if ((len == 128) || (len == 64)) t = ((int32_t)zeta * r[j + len]);
				else t = fnt_reduce((int32_t)zeta * r[j + len]);
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
	const int16_t f = 255; // 1/128

	k = 127;
	for (len = 2; len <= 128; len <<= 1) {
		for (start = 0; start < 256; start = j + len)
		{
			if (len == 128) zeta = fnt_q1_zetas[0];
			else zeta = fnt_q1_zetas[k--];
			for (j = start; j < start + len; j++)
			{
				t = r[j];
				r[j] = (t + r[j + len]);
				r[j + len] = (r[j + len] - t);
				r[j + len] = fnt_reduce((int32_t)zeta * r[j + len]);
			}
		}
	}

	for (int cnt_i = 0; cnt_i < 128; cnt_i++)
	{
		r[cnt_i] = fnt_reduce((int32_t)r[cnt_i] * f);
	}
}

void fnt_basemul(int32_t r[2], const int32_t a[2], const int32_t b[2], int16_t zeta)
{
    int32_t tmp[2] = {0};
  
    tmp[0] =  (int32_t)fnt_reduce(a[1]  * b[1]);
    tmp[0] =  (int32_t)fnt_reduce(tmp[0]* zeta);
    tmp[0] += (int32_t)fnt_reduce(a[0]  * b[0]);
    tmp[1] =  (int32_t)fnt_reduce(a[0]  * b[1]);
    tmp[1] += (int32_t)fnt_reduce(a[1]  * b[0]);
    
    r[0] = tmp[0]; r[1] = tmp[1];
}

void point_mul_small(int32_t des[N], int32_t src1[N], int32_t src2[N])
{
  for(int cnt_i = 0 ; cnt_i < N / 4 ; cnt_i++)
  {
    fnt_basemul(&des[4 * cnt_i    ], &src1[4 * cnt_i    ], &src2[4 * cnt_i    ], fnt_q1_zetas[64 + cnt_i]);
    fnt_basemul(&des[4 * cnt_i + 2], &src1[4 * cnt_i + 2], &src2[4 * cnt_i + 2], -fnt_q1_zetas[64 + cnt_i]);
  } 
}