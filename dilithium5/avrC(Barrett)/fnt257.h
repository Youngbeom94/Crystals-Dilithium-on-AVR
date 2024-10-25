#ifndef FNT_H
#define FNT_H

#include "ntt.h"

#define FNT_Q1 (257)
#define FNT_Q1HF (257 / 2)
#define FNT_Q1HFNEG -(257 / 2)


extern __flash int16_t fnt_q1_zetas[128];
int16_t fnt_reduce(int32_t src);
void fnt_ntt(int32_t r[N]);
void fnt_invntt(int32_t r[N]);
void point_mul_small(int32_t des[N], int32_t src1[N], int32_t src2[N]);
int32_t fnt_freeze(int32_t src);

#endif
