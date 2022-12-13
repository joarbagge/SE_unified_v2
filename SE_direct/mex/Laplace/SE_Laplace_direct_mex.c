#include "mex.h"
#include "SE_Laplace_direct.h"
#include "../SE_direct_common.c"

#define X   prhs[0] // Source locations
#define F   prhs[1] // Source strengths
#define OPT prhs[2] // Ewald options

#define U   plhs[0] // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

// Select what to compute from command-line defined variables.
// Typical compiler defines:
// (...) -DEWALD_REAL -DTHREE_PERIODIC
#ifdef EWALD_REAL
#ifdef THREE_PERIODIC
#define EWALD_KERNEL SE3P_Laplace_direct_real
#define EWALD_TAG "3P_REAL"
#elif TWO_PERIODIC
#define EWALD_KERNEL SE2P_Laplace_direct_real
#define EWALD_TAG "2P_REAL"
#elif ONE_PERIODIC
#define EWALD_KERNEL SE1P_Laplace_direct_real
#define EWALD_TAG "1P_REAL"
#elif ZERO_PERIODIC
#define EWALD_KERNEL SE0P_Laplace_direct_real
#define EWALD_TAG "0P_REAL"
#else
#error "Must give -D<ZERO/ONE/TWO/THREE>_PERIODIC to compiler"
#endif
#endif

#ifdef EWALD_FULL
#ifdef ONE_PERIODIC
#define EWALD_KERNEL SE1P_Laplace_direct_full
#define EWALD_TAG "1P_FULL"
#elif ZERO_PERIODIC
#define EWALD_KERNEL SE0P_Laplace_direct_full
#define EWALD_TAG "0P_FULL"
#else
#error "Must give -D<ZERO/ONE>_PERIODIC to compiler"
#endif
#endif

#ifdef EWALD_FOURIER
#ifdef THREE_PERIODIC
#define EWALD_KERNEL SE3P_Laplace_direct_fourier
#define EWALD_TAG "3P_FOURIER"
#elif TWO_PERIODIC
#define EWALD_KERNEL SE2P_Laplace_direct_fourier
#define EWALD_TAG "2P_FOURIER"
#elif ONE_PERIODIC
#define EWALD_KERNEL SE1P_Laplace_direct_fourier
#define EWALD_TAG "1P_FOURIER"
#else
#error "Must give -D<ONE/TWO/THREE>_PERIODIC to compiler"
#endif
#endif

#ifdef EWALD_FOURIER_K0
#ifdef TWO_PERIODIC
#define EWALD_KERNEL SE2P_Laplace_direct_fourier_k0
#define EWALD_TAG "2P_FOURIER_K0"
#elif ONE_PERIODIC
#define EWALD_KERNEL SE1P_Laplace_direct_fourier_k0
#define EWALD_TAG "1P_FOURIER_K0"
#else
#error "Must give -D<ONE/TWO>_PERIODIC to compiler"
#endif
#endif

#ifdef EWALD_SELF
#ifdef THREE_PERIODIC
#define EWALD_KERNEL SE3P_Laplace_direct_self
#define EWALD_TAG "3P_SELF"
#elif TWO_PERIODIC
#define EWALD_KERNEL SE2P_Laplace_direct_self
#define EWALD_TAG "2P_SELF"
#elif ONE_PERIODIC
#define EWALD_KERNEL SE1P_Laplace_direct_self
#define EWALD_TAG "1P_SELF"
#else
#error "Must give -D<ONE/TWO/THREE>_PERIODIC to compiler"
#endif
#endif

#ifdef GRADIENT
#define GRADIENT_TAG "GRAD"
#else
#define GRADIENT_TAG "POT"
#endif

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] )
{
    // Check input
    if (nrhs != 3)
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "Wrong number of input arguments (x, f, opt)");

    // Input dimensions, X should be N-by-3, F should be N-by-1
    if (mxGetNumberOfElements(X) > 0 && mxGetN(X) != 3)
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "x must be N-by-3");
    const size_t N = mxGetM(X);
    if (mxGetM(F) != N || mxGetN(F) > 1)
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "f must be N-by-1");

    // Parameter values
    ewald_opts opt;
    unpack_mex_opt(&opt, OPT);
    // Check options
#ifdef THREE_PERIODIC
    if (opt.box[0] <= 0 || opt.box[1] <= 0 || opt.box[2] <= 0)
#elif TWO_PERIODIC
    if (opt.box[0] <= 0 || opt.box[1] <= 0)
#elif ONE_PERIODIC
    if (opt.box[0] <= 0)
#else
    if (false)
#endif
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "opt.box missing or invalid");
#ifndef EWALD_FULL
    if (opt.xi <= 0)
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "opt.xi missing or invalid");
#endif
#ifdef CUTOFF
    if (opt.rc <= 0)
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "opt.rc missing or invalid");
#endif
#if !defined(ZERO_PERIODIC) && !defined(EWALD_FOURIER_K0) && !defined(EWALD_SELF)
    if (opt.layers == -1)
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "opt.layers missing or invalid");
#endif

    // Pointers
#ifdef EXTERNAL
    size_t num_eval = opt.N_ext;
#else
    size_t num_eval = opt.N_eval;
    if (num_eval == -1) num_eval = N;
#endif

#ifndef EWALD_SELF
    const double* x = mxGetPr(X);
#endif
    const double* f = mxGetPr(F);
#ifdef GRADIENT
    // Allocate output matrix for the potential gradient
    U = mxCreateDoubleMatrix(num_eval, 3, mxREAL);
#else
    // Allocate output matrix for the potential
    U = mxCreateDoubleMatrix(num_eval, 1, mxREAL);
#endif
    double* restrict u = mxGetPr(U);

    if (VERBOSE)
    {
#ifdef EXTERNAL
        __PRINTF("\n[DIRECT LAPLACE EWALD (%s EXT %s)] MEX N=(%d,%d) ", EWALD_TAG, GRADIENT_TAG, N, num_eval);
#elif CUTOFF
        __PRINTF("\n[DIRECT LAPLACE EWALD (%s CUT %s)] MEX N=(%d,%d) ", EWALD_TAG, GRADIENT_TAG, N, num_eval);
#else
        __PRINTF("\n[DIRECT LAPLACE EWALD (%s %s)] MEX N=(%d,%d) ", EWALD_TAG, GRADIENT_TAG, N, num_eval);
#endif
        __PRINTF("xi=%.2f layers=%d\n", opt.xi, opt.layers);
    }

    // Call kernel
#ifdef EWALD_SELF
    EWALD_KERNEL(u, num_eval, f, N, &opt);
#else
    EWALD_KERNEL(u, num_eval, x, f, N, &opt);
#endif

    free_ewald_opts(&opt);
}
