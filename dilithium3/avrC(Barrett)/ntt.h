#ifndef NTT_H
#define NTT_H
#include "params.h"
#include "barrett_arithmetic.h"
#include "fnt769.h"
#include <stdint.h>

void ntt(int32_t a[N]);

void invntt_tomont(int32_t a[N]);

void invntt(int32_t a[N]);

void ntt_small(int32_t a[N]);

void invntt_small(int32_t a[N]);

#endif
