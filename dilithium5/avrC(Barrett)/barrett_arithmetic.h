#ifndef BARRETT_ARITHMETIC_H
#define BARRETT_ARITHMETIC_H


#include "NTT_params.h"
#include "ntt.h"
#include "params.h"
#include <stdint.h>

// we assume -Q / 2 < a < Q / 2
int32_t gethi(int32_t a);

int32_t freeze_32(int32_t a);

int32_t Barrett_mul_approx(int32_t a, int32_t b, int32_t bp, int32_t q);

void point_mul_pre(int32_t res[N], int32_t src1[N]);

void point_mul_on_the_fly(int32_t des[N], int32_t src1[N], int32_t src2[N]);

void barrett_ntt(int32_t a[N]);

void barrett_invntt(int32_t a[N]);
#endif 
