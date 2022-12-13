#include "SE_direct_common.h"

// MATLAB (one-based, doubles) to C (zero-based, integers) index translation
void index_translation(size_t* idx, const double* idx_d, size_t N)
{
    for (size_t n=0; n<N; n++)
    {
        idx[n] = (size_t)idx_d[n] - 1;
    }
}

// MEX option unpacking
void unpack_mex_opt(ewald_opts* opt, const mxArray* mx_opt)
{
    // box: mandatory except for ZERO_PERIODIC
    opt->box[0] = -1; opt->box[1] = -1; opt->box[2] = -1;
    const mxArray * mx_box = mxGetField(mx_opt, 0, "box");
    if (mx_box)
    {
#ifndef ZERO_PERIODIC
        double * box = mxGetPr(mx_box);
#endif
#ifdef THREE_PERIODIC
        opt->box[0] = box[0];
        opt->box[1] = box[1];
        opt->box[2] = box[2];
#elif TWO_PERIODIC
        opt->box[0] = box[0];
        opt->box[1] = box[1];
#elif ONE_PERIODIC
        opt->box[0] = box[0];
#elif ZERO_PERIODIC
#else
#error "Must give -D<ZERO/ONE/TWO/THREE>_PERIODIC to compiler"
#endif
    }

    // xi: mandatory except for EWALD_FULL
    opt->xi = -1;
    const mxArray * mx_xi = mxGetField(mx_opt, 0, "xi");
    if (mx_xi) opt->xi = mxGetScalar(mx_xi);

    // rc: mandatory for REAL_CUT Ewald sum
    opt->rc = -1;
    const mxArray * mx_rc = mxGetField(mx_opt, 0, "rc");
    if (mx_rc) opt->rc = mxGetScalar(mx_rc);

    // layers: mandatory for Ewald sums that are truncated
    opt->layers = -1;
    const mxArray * mx_layers = mxGetField(mx_opt, 0, "layers");
    if (mx_layers) opt->layers = (size_t)mxGetScalar(mx_layers);

    // eval_idx: optional (default is to evaluate at all source points)
    opt->N_eval = -1;
    opt->eval_idx = NULL;
    const mxArray * mx_eval_idx = mxGetField(mx_opt, 0, "eval_idx");
    if (mx_eval_idx)
    {
        opt->N_eval = mxGetNumberOfElements(mx_eval_idx);
        const double * eval_idx = mxGetPr(mx_eval_idx);
        opt->eval_idx = (size_t*)__MALLOC(sizeof(size_t)*opt->N_eval);
        index_translation(opt->eval_idx, eval_idx, opt->N_eval);
    }

    // eval_ext_x: optional (default is empty list)
    opt->N_ext = 0;
    opt->eval_ext_x = NULL;
    const mxArray * mx_eval_ext_x = mxGetField(mx_opt, 0, "eval_ext_x");
    if (mx_eval_ext_x)
    {
        opt->N_ext = mxGetM(mx_eval_ext_x);
        if (opt->N_ext > 0 && mxGetN(mx_eval_ext_x) != 3)
            mexErrMsgIdAndTxt("SE:DirectMEXOpt", "opt.eval_ext_x must be N_ext-by-3");
        opt->eval_ext_x = mxGetPr(mx_eval_ext_x);
    }

    // stokeslet_k0_constant: optional (defaults to zero)
    opt->stokeslet_k0_constant = 0;
    const mxArray * mx_s0c = mxGetField(mx_opt, 0, "stokeslet_k0_constant");
    if (mx_s0c) opt->stokeslet_k0_constant = mxGetScalar(mx_s0c);
}

// Free Ewald options
void free_ewald_opts(ewald_opts* opt)
{
    __FREE(opt->eval_idx);
    // opt->eval_ext_x was allocated by the MEX library itself, no need to free it
}
