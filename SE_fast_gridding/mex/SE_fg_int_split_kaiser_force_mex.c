#include "mex.h"
#include "../src/SE_fgg.h"
#include "../src/SE_fkg.h"

void SE_FGG_MEX_params(SE_FGG_params*, const mxArray*, int);

#define X   prhs[0]
#define HH  prhs[1]
#define OPT prhs[2]
#define ZX  prhs[3]
#define ZY  prhs[4]
#define ZZ  prhs[5]
#define ZFX prhs[6]
#define ZFY prhs[7]
#define ZFZ prhs[8]
#define IDX prhs[9]

#define FORCE_OUT plhs[0]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] )
{

    const int N = mxGetM(X);
    double* restrict x = mxGetPr(X);
    const double* H_per = mxGetPr(HH);

    SE_FGG_params params;
    SE_FGG_MEX_params(&params, OPT, N);

    // scratch arrays
    SE_FGG_work work;
    SE_FKG_allocate_workspace(&work, &params, false);

    // attach pre-computed quantities
    work.zx = mxGetPr(ZX);
    work.zy = mxGetPr(ZY);
    work.zz = mxGetPr(ZZ);
    work.zfx = mxGetPr(ZFX);
    work.zfy = mxGetPr(ZFY);
    work.zfz = mxGetPr(ZFZ);
    work.idx = (int*)mxGetData(IDX);

    // output vector
    FORCE_OUT = mxCreateDoubleMatrix(N,3,mxREAL);
    double* force = mxGetPr(FORCE_OUT);

    // coordinates and charges
    SE_state st = {.x = x, .q = NULL};

    if(VERBOSE)
        mexPrintf("[SE%s FG(i)] N=%d, P=%d\n",PER_STR,N,params.P);

    // now do the work
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
    {
#ifdef THREE_PERIODIC
        SE_FGG_extend_fcn(&work, H_per, &params);
#endif
#ifdef TWO_PERIODIC
        SE2P_FGG_extend_fcn(&work, H_per, &params);
#endif
#ifdef ONE_PERIODIC
        SE1P_FGG_extend_fcn(&work, H_per, &params);
#endif
#ifdef __AVX__
        //SE_FKG_int_split_AVX_dispatch_force(force, &st, &work, &params); // TODO
#else
        //SE_FKG_int_split_SSE_dispatch_force(force, &st, &work, &params); // TODO
#endif
        SE_FKG_int_split_force(force, &st, &work, &params); // slower code
    }

    // done
    SE_FGG_free_workspace(&work);
}
