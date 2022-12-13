#ifndef __INCOMPLETE_BESSEL_K_H_
#define __INCOMPLETE_BESSEL_K_H_

#include <math.h>
#include <immintrin.h>
#if (defined __AVX__ || defined __SSE4_2__)
#include "math_x86.h"
#endif

double incomplete_bessel_K(int n, double x, double y);

#endif
