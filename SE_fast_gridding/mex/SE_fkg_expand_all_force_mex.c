#include "mex.h"
#include "../src/SE_fgg.h"
#include "../src/SE_fkg.h"

void SE_FGG_MEX_params(SE_FGG_params*, const mxArray*, int);

#define X   prhs[0] 
#define OPT prhs[1] 

#define ZX  plhs[0]  // Output
#define ZY  plhs[1]  // Output
#define ZZ  plhs[2]  // Output
#define ZFX plhs[3]  // Output
#define ZFY plhs[4]  // Output
#define ZFZ plhs[5]  // Output
#define IDX plhs[6]  // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

// Round up towards nearest multiple of 4
int ceil_mult_4(int n)
{
  const int mod = n % 4;
  if (mod == 0)
    return n;
  else
    return n + (4-mod);
}

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
    const int N = mxGetM(X);
    double* restrict x = mxGetPr(X);

    // pack parameters
    SE_FGG_params params;
    SE_FGG_MEX_params(&params, OPT, N);

    // Evaluating the polynomial window is much faster if P is a
    // multiple of 4, so round up for the evaluation.
    if (params.use_polynomial_window) {
      params.P_eval_window = ceil_mult_4(params.P);
    }

    // allocate output array
    ZX = mxCreateDoubleMatrix(params.P_eval_window,N,mxREAL);
    ZY = mxCreateDoubleMatrix(params.P_eval_window,N,mxREAL);
    ZZ = mxCreateDoubleMatrix(params.P_eval_window,N,mxREAL);
    ZFX = mxCreateDoubleMatrix(params.P_eval_window,N,mxREAL);
    ZFY = mxCreateDoubleMatrix(params.P_eval_window,N,mxREAL);
    ZFZ = mxCreateDoubleMatrix(params.P_eval_window,N,mxREAL);

    // output
    const size_t dims[2] = {N,1};
    IDX = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);

    // wrap in SE_work struct
    SE_FGG_work work;
    work.zx = mxGetPr(ZX);
    work.zy = mxGetPr(ZY);
    work.zz = mxGetPr(ZZ);
    work.zfx = mxGetPr(ZFX);
    work.zfy = mxGetPr(ZFY);
    work.zfz = mxGetPr(ZFZ);
    work.idx = (int*)mxGetData(IDX);

    // coordinates and charges
    const SE_state st = {.x = x,  .q = NULL};

    if(VERBOSE)
	mexPrintf("[SE%s FG(E)] N=%d, P=%d\n",PER_STR,N,params.P);

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
    {
	// now do the work (COMPILED for different periodicities)
	SE_FKG_expand_all_force(&work, &st, &params);
    }
    // done
}
