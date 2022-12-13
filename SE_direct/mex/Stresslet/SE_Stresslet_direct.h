#ifndef SE_STRESSLET_DIRECT_H
#define SE_STRESSLET_DIRECT_H

#include "../SE_direct_common.h"

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define __MALLOC mxMalloc
#define __PRINTF mexPrintf
#define __FREE mxFree
#else
#define __MALLOC malloc
#define __PRINTF printf
#define __FREE free
#endif

#define SE_DECLARE_DIRECT(NAME) void NAME(double*, const size_t, \
    const double*, const double*, const double*, const size_t, const ewald_opts*)

SE_DECLARE_DIRECT(SE3P_Stresslet_direct_real);
SE_DECLARE_DIRECT(SE3P_Stresslet_direct_fourier);
SE_DECLARE_DIRECT(SE3P_Stresslet_direct_fourier_k0);

SE_DECLARE_DIRECT(SE2P_Stresslet_direct_real);
SE_DECLARE_DIRECT(SE2P_Stresslet_direct_fourier);
SE_DECLARE_DIRECT(SE2P_Stresslet_direct_fourier_k0);

SE_DECLARE_DIRECT(SE1P_Stresslet_direct_full);
SE_DECLARE_DIRECT(SE1P_Stresslet_direct_real);
SE_DECLARE_DIRECT(SE1P_Stresslet_direct_fourier);
SE_DECLARE_DIRECT(SE1P_Stresslet_direct_fourier_k0);

SE_DECLARE_DIRECT(SE0P_Stresslet_direct_full);
SE_DECLARE_DIRECT(SE0P_Stresslet_direct_real);

#endif
