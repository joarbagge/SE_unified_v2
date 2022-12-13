#include <stdio.h>
#include <string.h>
#include <mex.h>
#include "SE_fgg.h"
#include "SE_fg_windows.c"

// Get field s from matlab struct p. Abort if field is missing.
void* get_arg(const mxArray* p, const char* s)
{
    mxArray* x = mxGetField(p,0,s);
    if (x) return mxGetData(x);
    else
    {
        mexErrMsgTxt("Missing mandatory parameter in struct");
        return (void*) NULL;
    }
}
void* get_arg_str(const mxArray* p, const char* s)
{
    mxArray* x = mxGetField(p,0,s);
    if (x) return mxArrayToString(x);
    else
    {
        mexErrMsgTxt("Missing mandatory parameter in struct");
        return (void*) NULL;
    }
}

KaiserlikeWindow get_kaiserlike_window(const char* window_id)
{
    if (strcmp(window_id, "expsemicirc") == 0) {
        return &SE_fg_window_expsemicirc;
    } else if (strcmp(window_id, "kaiser_exact") == 0) {
        return &SE_fg_window_kaiser_exact;
    } else if (strcmp(window_id, "kaiser_poly") == 0) {
        return &SE_fg_window_kaiser_poly;
    } else {
        mexErrMsgTxt("Unsupported window function");
        return (KaiserlikeWindow) NULL;
    }
}

KaiserlikeWindow get_kaiserlike_window_deriv(const char* window_id)
{
    if (strcmp(window_id, "kaiser_exact") == 0) {
        return &SE_fg_window_kaiser_exact_deriv;
    } else if (strcmp(window_id, "kaiser_poly") == 0) {
        return &SE_fg_window_kaiser_poly_deriv;
    } else {
        mexErrMsgTxt("Unsupported window function");
        return (KaiserlikeWindow) NULL;
    }
}

// Get relevant parameters from mxArray and populate FGG params struct.
void SE_FGG_MEX_params(SE_FGG_params* params, const mxArray* OPT, int N)
{
    const double* m    = (double*) get_arg(OPT, "grid");
    const double* p    = (double*) get_arg(OPT, "window_P");
    const double* h    = (double*) get_arg(OPT, "h");
#ifdef KAISER
    const double* beta = (double*) get_arg(OPT, "window_shape");
    params->beta = *beta;
    char* window = (char*) get_arg_str(OPT, "window");
    params->window = get_kaiserlike_window(window);
#ifdef FORCE
    params->window_deriv = get_kaiserlike_window_deriv(window);
#endif
    if (strcmp(window, "kaiser_poly") == 0) {
        params->use_polynomial_window = 1;
        params->polynomial_degree = (int) * (double*) get_arg(OPT, "polynomial_degree");
    } else {
        params->use_polynomial_window = 0;
    }
    mxFree(window);
#else
    const double w = *(double*) get_arg(OPT, "window_halfwidth");
    const double alpha = *(double*) get_arg(OPT, "window_shape");
    params->c = (alpha) / (w*w);
    params->d = pow(params->c/PI, 1.5);
#endif

    params->N = N;
    params->P = (int) *p;
    params->P_half=half(params->P);
    params->P_eval_window = params->P;
    params->h = *h;

#ifdef THREE_PERIODIC
    params->a = -FGG_INF;
    params->dims[0] = (int) m[0];
    params->dims[1] = (int) m[1];
    params->dims[2] = (int) m[2];
    params->npdims[0] = params->dims[0]+params->P;
    params->npdims[1] = params->dims[1]+params->P;
    params->npdims[2] = params->dims[2]+params->P;
#elif TWO_PERIODIC
    const double* mz  = (double*) get_arg(OPT, "grid3");
    const double* a   = (double*) get_arg(OPT, "free_offset"); /* z offset */
    params->a = a[0];
    params->dims[0] = (int)  m[0];
    params->dims[1] = (int)  m[1];
    params->dims[2] = (int) mz[0];
    params->npdims[0] = params->dims[0]+params->P;
    params->npdims[1] = params->dims[1]+params->P;
    params->npdims[2] = params->dims[2];
#elif ONE_PERIODIC
    const double* my  = (double*) get_arg(OPT, "grid2");
    const double* mz  = (double*) get_arg(OPT, "grid3");
    /* y and z offsets */
    const double* a   = (double*) get_arg(OPT, "free_offset");
    params->a = a[0];
    params->b = a[1];
    params->dims[0] = (int)  m[0];
    params->dims[1] = (int) my[0];
    params->dims[2] = (int) mz[0];
    params->npdims[0] = params->dims[0]+params->P;
    params->npdims[1] = params->dims[1];
    params->npdims[2] = params->dims[2];
#elif ZERO_PERIODIC
// NB: It seems this is never used; THREE_PERIODIC is used in the 0P case
    /* x, y and z offsets */
    const double* a   = (double*) get_arg(OPT, "free_offset");
    params->a = a[0]; // FIXME: we assume the same offset in each direction!
    params->dims[0] = (int) m[0];
    params->dims[1] = (int) m[1];
    params->dims[2] = (int) m[2];
    params->npdims[0] = params->dims[0];
    params->npdims[1] = params->dims[1];
    params->npdims[2] = params->dims[2];
#endif
}
