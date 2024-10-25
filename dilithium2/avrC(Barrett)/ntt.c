#include "ntt.h"
#include "params.h"
#include "reduce.h"
#include <stdint.h>

void ntt(int32_t a[N]) 
{
    barrett_ntt(a);
    for (int cnt_i = 0; cnt_i < N; cnt_i++)  a[cnt_i] = freeze_32(a[cnt_i]);
}

void invntt_tomont(int32_t a[N]) 
{
    barrett_invntt(a);
}

void invntt(int32_t a[N])
{
    barrett_invntt(a);
}


void ntt_small(int32_t a[N])
{
    fnt_ntt(a);
}

void invntt_small(int32_t a[N])
{
    fnt_invntt(a);
}