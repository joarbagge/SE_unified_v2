#include <mex.h>
#include <math.h>
#include <string.h>
#include "SE_real_space_mex.h"
#define INV_SQRT_PI 0.564189583547756286948079 // = 1/sqrt(pi)
#define INV_8_PI 0.0397887357729738339422209 // = 1/(8*pi)

#define POT_OUT plhs[0] // potential

/* For documentation, see matlab/src/SpectralEwald/SE_real_space.m. */

// Entry point
// Call signature #1: pot = SE_real_space_mex("stokes_sl", opt, targets, sources, density)
// Call signature #2: pot = SE_real_space_mex("stokes_sl", ..., delta)
// Call signature #3: pot = SE_real_space_mex("stokes_dl", opt, targets, sources, normals, density)
// Call signature #4: pot = SE_real_space_mex("stokes_dl", ..., delta)
// Call signature #5: pot = SE_real_space_mex("stokes_dl", ..., delta, smoothing)
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Argument handling
    bool double_layer = false;
    if (nrhs < 1) {
        mexErrMsgIdAndTxt("CapsMex:NotEnoughInputs",
                          "Not enough input arguments.");
    }
    const char* restrict kernel = mxArrayToString(prhs[0]);
    kernel_fun_t kernel_fun = NULL;
    if (strcmp(kernel, "stokes_sl") == 0) {
        kernel_fun = &stokes_sl_kernel;
    } else if (strcmp(kernel, "stokes_dl") == 0) {
        double_layer = true;
        kernel_fun = &stokes_dl_kernel;
    } else {
        mexErrMsgIdAndTxt("CapsMex:UnknownOption",
                          "Unsupported kernel \"%s\"", kernel);
    }

    if (nrhs < 5+double_layer) {
        mexErrMsgIdAndTxt("CapsMex:NotEnoughInputs",
                          "Not enough input arguments, at least %d expected.",
                          5+double_layer);
    }

    Options opt;
    load_options(prhs[1], &opt);

    const double* restrict targets = mxGetPr(prhs[2]);
    const double* restrict sources = mxGetPr(prhs[3]);
    const double* restrict normals = NULL;
    if (double_layer) {
        normals = mxGetPr(prhs[4]);
    }
    const double* restrict density = mxGetPr(prhs[4+double_layer]);

    // Default values for optional arguments
    const double default_delta = 0;
    const double* restrict delta;
    size_t numel_delta = 1;
    const char default_smoothing[] = "auto";
    const char* restrict smoothing;
    if (nrhs >= 6+double_layer) {
        delta = mxGetPr(prhs[5+double_layer]);
        numel_delta = mxGetNumberOfElements(prhs[5+double_layer]);
    } else {
        delta = &default_delta;
    }
    if (nrhs >= 7+double_layer) {
        smoothing = mxArrayToString(prhs[6+double_layer]);
    } else {
        smoothing = default_smoothing;
    }
    if (!double_layer && strcmp(smoothing, "auto") != 0) {
        mexErrMsgIdAndTxt("CapsMex:UnknownOption",
                          "%s does not support smoothing \"%s\"",
                          kernel, smoothing);
    }
    // End argument handling

    // Check dimensions
    size_t Nt = mxGetM(prhs[2]);
    size_t Ns = mxGetM(prhs[3]);
    if (mxGetN(prhs[2]) != 3 || mxGetN(prhs[3]) != 3 || mxGetN(prhs[4+double_layer]) != 3) {
        mexErrMsgIdAndTxt("CapsMex:IncorrectDimensions",
                          "Incorrect number of columns in input (should be 3).");
    }
    if (double_layer && mxGetN(prhs[4]) != 3) {
        mexErrMsgIdAndTxt("CapsMex:IncorrectDimensions",
                          "Incorrect number of columns in input (should be 3).");
    }
    if (double_layer && mxGetM(prhs[4]) != Ns) {
        mexErrMsgIdAndTxt("CapsMex:IncorrectDimensions",
                          "Incorrect number of rows in normals (should match sources).");
    }
    if (mxGetM(prhs[4+double_layer]) != Ns) {
        mexErrMsgIdAndTxt("CapsMex:IncorrectDimensions",
                          "Incorrect number of rows in density (should match sources).");
    }
    if (numel_delta > 1 && numel_delta != Ns) {
        mexErrMsgIdAndTxt("CapsMex:IncorrectDimensions",
                          "Incorrect number of delta values given (should match sources).");
    }

    // Select smoothing function
    int num_terms = 4;
    Smoothing smooth[2];
    if (!double_layer) {
        smooth[0].fun = &smoothfun1;
        smooth[0].taylor = &smoothfun1_taylor;
        // For smoothfun1, 1 Taylor term is enough to get machine precision
        smooth[0].threshold = get_smoothfun1_threshold(1);
        smooth[0].terms = num_terms;
        smooth[1].fun = &smoothfun2;
        smooth[1].taylor = &smoothfun2_taylor;
        smooth[1].threshold = get_smoothfun2_threshold(num_terms);
        smooth[1].terms = num_terms;
    } else {
        bool full_smoothing = true;
        if (strcmp(smoothing, "simple") == 0) {
            full_smoothing = false;
        } else if (strcmp(smoothing, "auto") != 0 && strcmp(smoothing, "full") != 0) {
            mexErrMsgIdAndTxt("CapsMex:UnknownOption",
                              "Unknown smoothing function \"%s\"",
                              smoothing);
        }
        if (full_smoothing) {
            smooth[0].fun = &smoothfun3_full;
            smooth[0].taylor = &smoothfun3_full_taylor;
            smooth[0].threshold = get_smoothfun3_full_threshold(num_terms);
            smooth[0].terms = num_terms;
        } else {
            smooth[0].fun = &smoothfun3_simple;
            smooth[0].taylor = &smoothfun3_simple_taylor;
            smooth[0].threshold = get_smoothfun3_simple_threshold(num_terms);
            smooth[0].terms = num_terms;
        }
    }

    // Allocate output
    POT_OUT = mxCreateDoubleMatrix(Nt, 3, mxREAL);
    double* restrict pot = mxGetPr(POT_OUT);
    int num_shifts = get_num_shifts(opt.periodicity);

    // Compute potential in helper function
    core_computation(pot, kernel_fun, double_layer, &opt, targets, Nt,
                     sources, Ns, normals, density, delta, numel_delta,
                     smooth, num_shifts);
}

void load_options(const mxArray* s, Options* restrict opt)
{
    // Periodicity
    mxArray* x = mxGetField(s, 0, "periodicity");
    if (!x) {
        mexErrMsgIdAndTxt("CapsMex:NotEnoughInputs",
                          "Mandatory field \"periodicity\" missing from option struct.");
    }
    const double* per = (double*) mxGetData(x);
    opt->periodicity = (int) *per;
    // Box
    x = mxGetField(s, 0, "box");
    if (!x) {
        mexErrMsgIdAndTxt("CapsMex:NotEnoughInputs",
                          "Mandatory field \"box\" missing from option struct.");
    }
    if (mxGetNumberOfElements(x) != 3) {
        mexErrMsgIdAndTxt("CapsMex:NotEnoughInputs",
                          "Field \"box\" should be a vector of exactly size three.");
    }
    const double* box = (double*) mxGetData(x);
    opt->box[0] = box[0];
    opt->box[1] = box[1];
    opt->box[2] = box[2];
    // Ewald split parameter xi
    x = mxGetField(s, 0, "xi");
    if (!x) {
        mexErrMsgIdAndTxt("CapsMex:NotEnoughInputs",
                          "Mandatory field \"xi\" missing from option struct.");
    }
    const double* xi = (double*) mxGetData(x);
    opt->xi = *xi;
    // Cutoff distance rc
    x = mxGetField(s, 0, "rc");
    if (!x) {
        mexErrMsgIdAndTxt("CapsMex:NotEnoughInputs",
                          "Mandatory field \"rc\" missing from option struct.");
    }
    const double* rc = (double*) mxGetData(x);
    opt->rc = *rc;
}

int get_num_shifts(const int periodicity)
{
    int shifts = 1;
    for (int i=0; i<periodicity; i++) {
        shifts *= 3;
    }
    return shifts;
}

void core_computation(double* restrict pot, kernel_fun_t kernel_fun, bool double_layer,
                      const Options* restrict opt, const double* restrict targets,
                      size_t Nt, const double* restrict sources, size_t Ns,
                      const double* restrict normals, const double* restrict density,
                      const double* restrict delta, size_t numel_delta,
                      const Smoothing* restrict smooth, int num_shifts)
{
    double xi = opt->xi;
    double xi2 = xi*xi;
    double rc = opt->rc;
    double rc2 = rc*rc;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t j=0; j<Nt; j++) { // loop over targets
        double target[] = {targets[j], targets[j+Nt], targets[j+Nt*2]};
        double rvec[3], rshvec[3], r2, r, delt;
        double xir, xir2, xiexp, rdotf, A, B, C, term1, term2, term3, term4;
        int* shift;
        double sum[] = {0, 0, 0};
        double normal[] = {0, 0, 0};
        for (size_t i=0; i<Ns; i++) { // loop over sources
            double source[] = {sources[i], sources[i+Ns], sources[i+Ns*2]};
            if (double_layer) {
                for (size_t d=0; d<3; d++) {
                    normal[d] = normals[i+Ns*d];
                }
            }
            double dens[] = {density[i], density[i+Ns], density[i+Ns*2]};
            if (numel_delta == 1) {
                delt = delta[0];
            } else {
                delt = delta[i];
            }
            for (size_t d=0; d<3; d++) {
                rvec[d] = target[d] - source[d];
            }
            for (size_t k=0; k<num_shifts; k++) { // loop over shifts
                shift = ALL_SHIFTS[k];
                for (size_t d=0; d<3; d++) {
                    rshvec[d] = rvec[d] + shift[d]*opt->box[d];
                }
                // Check if source is within rc of target
                r2 = dot3(rshvec, rshvec);
                if (r2 > rc2) continue;
                r = sqrt(r2);
                // Common computations (regardless of kernel)
                xir = xi*r;
                xir2 = xi2*r2;
                xiexp = xi * exp(-xir2);
                rdotf = dot3(rshvec, dens);
                // Replace A by its limit close to zero (relative error at 1e-14
                // should be around 4e-29, so this is very accurate)
                if (xir < 1e-14) {
                    A = 2*xi*INV_SQRT_PI;
                } else {
                    A = erf(xir)/r;
                }
                B = 2*xiexp*INV_SQRT_PI;
                // Replace C by its Taylor expansion close to zero
                // (relative error at most around 3e-13 for four terms)
                // NOTE: This might be a bit unnecessary, as C goes into terms
                // that will anyway go to zero as r -> 0. But we do it
                // anyway for good measure.
                if (xir < 4.75e-2) {
                    term1 = xi2*xi*INV_SQRT_PI;
                    term2 = xir2*term1;
                    term3 = xir2*term2;
                    term4 = xir2*term3;
                    C = (4.0/3.0)*term1 - (4.0/5.0)*term2 + (2.0/7.0)*term3 - (2.0/27.0)*term4;
                } else {
                    C = (A-B)/r2;
                }
                // Kernel-dependent computations
                (*kernel_fun)(sum, xi, xi2, rshvec, r, r2, normal, dens, rdotf,
                              delt, A, B, C, smooth);
            }
        }
        pot[j     ] = sum[0] * INV_8_PI;
        pot[j+Nt  ] = sum[1] * INV_8_PI;
        pot[j+Nt*2] = sum[2] * INV_8_PI;
    }
}

void stokes_sl_kernel(double* restrict pot, const double xi, const double xi2,
                      const double* restrict rvec, const double r, const double r2,
                      const double* restrict n, const double* restrict f,
                      const double rdotf, const double delt,
                      const double A, const double B, const double C,
                      const Smoothing* smooth)
{
    // Singular part, i.e., terms coming from the full (free-space) stokeslet
    double s1=0, s2=0;
    if (delt == 0) {
        // Punctured trapezoidal rule
        // Remove point where r==0 (1e-14 is ad hoc)
        if (r >= 1e-14) {
            s1 = 1/r;
            s2 = 1/(r2*r);
        }
    } else {
        // Beale regularization
        double rd = r/delt;
        // Replace by known Taylor expansion for r close to zero
        if (rd < smooth[0].threshold) {
            s1 = (*smooth[0].taylor)(delt, rd, smooth[0].terms);
        } else {
            s1 = (*smooth[0].fun)(rd) / r;
        }
        if (rd < smooth[1].threshold) {
            s2 = (*smooth[1].taylor)(delt, rd, smooth[1].terms);
        } else {
            s2 = (*smooth[1].fun)(rd) / (r2*r);
        }
    }
    // Sum up all terms
    double tmp = s1 - B - A;
    double t1[] = {tmp*f[0], tmp*f[1], tmp*f[2]};
    tmp = rdotf*(s2-C);
    double t2[] = {tmp*rvec[0], tmp*rvec[1], tmp*rvec[2]};
    for (size_t d=0; d<3; d++) {
        pot[d] += t1[d] + t2[d];
    }
}

void stokes_dl_kernel(double* restrict pot, const double xi, const double xi2,
                      const double* restrict rvec, const double r, const double r2,
                      const double* restrict n, const double* restrict f,
                      const double rdotf, const double delt,
                      const double A, const double B, const double C,
                      const Smoothing* smooth)
{
    double r4 = r2*r2;
    double rdotn = dot3(rvec, n);
    double fdotn = dot3(f, n);
    // Singular part, i.e., terms coming from the full (free-space) stresslet
    double s=0;
    if (delt == 0) {
        // Punctured trapezoidal rule
        // Remove point where r==0 (1e-14 is ad hoc)
        if (r >= 1e-14) {
            s = 1/(r4*r);
        }
    } else {
        // Beale regularization
        double rd = r/delt;
        // Replace by known Taylor expansion for r close to zero
        if (rd < smooth[0].threshold) {
            s = (*smooth[0].taylor)(delt, rd, smooth[0].terms);
        } else {
            s = (*smooth[0].fun)(rd) / (r4*r);
        }
    }
    // Sum up all terms
    double D = 2*xi2*B;
    // At r==0, this part will become zero (1e-14 is ad hoc)
    double tmp = 0;
    if (r >= 1e-14) {
        tmp = (6*C-2*D)*rdotf*rdotn/r2;
    }
    tmp += -6*s*rdotf*rdotn + D*fdotn;
    double t1[] = {tmp*rvec[0], tmp*rvec[1], tmp*rvec[2]};
    tmp = D*rdotn;
    double t2[] = {tmp*f[0], tmp*f[1], tmp*f[2]};
    tmp = D*rdotf;
    double t3[] = {tmp*n[0], tmp*n[1], tmp*n[2]};
    for (size_t d=0; d<3; d++) {
        pot[d] += t1[d] + t2[d] + t3[d];
    }
}

double smoothfun1(const double x)
{
    double x2 = x*x;
    double b = 2.0*x2 - 5.0;
    double c = 2.0/3.0;
    return erf(x) - c*INV_SQRT_PI*x*b*exp(-x2);
}

double smoothfun1_taylor(const double d, const double rd, const int terms)
{
    double term = INV_SQRT_PI/d;
    double s1 = (16.0/3.0)*term;
    if (terms == 1) {
        return s1;
    }
    double rd2 = rd*rd;
    term = term * rd2;
    s1 += -(16.0/3.0)*term;
    if (terms == 2) {
        return s1;
    }
    term = term * rd2;
    s1 += (16.0/5.0)*term;
    if (terms == 3) {
        return s1;
    }
    term = term * rd2;
    s1 += -(80.0/63.0)*term;
    if (terms == 4) {
        return s1;
    }
    return 0; // not supported
}

double get_smoothfun1_threshold(const int terms)
{
    if (terms == 1) {
        return 2.11e-8;
    } else if (terms == 2) {
        return 1.65e-4;
    } else if (terms == 3) {
        return 3.51e-3;
    } else if (terms == 4) {
        return 1.68e-2;
    } else {
        return 0; // not supported
    }
}

double smoothfun2(const double x)
{
    double x2 = x*x;
    double b = 4.0*x2*x2 - 14.0*x2 + 3.0;
    double c = 2.0/3.0;
    return erf(x) - c*INV_SQRT_PI*x*b*exp(-x2);
}

double smoothfun2_taylor(const double d, const double rd, const int terms)
{
    double term = INV_SQRT_PI/(d*d*d);
    double s2 = (32.0/3.0)*term;
    if (terms == 1) {
        return s2;
    }
    double rd2 = rd*rd;
    term = term * rd2;
    s2 += -(64.0/5.0)*term;
    if (terms == 2) {
        return s2;
    }
    term = term * rd2;
    s2 += (160.0/21.0)*term;
    if (terms == 3) {
        return s2;
    }
    term = term * rd2;
    s2 += -(80.0/27.0)*term;
    if (terms == 4) {
        return s2;
    }
    return 0; // not supported
}

double get_smoothfun2_threshold(const int terms)
{
    if (terms == 1) {
        return 9.13e-5;
    } else if (terms == 2) {
        return 2.21e-3;
    } else if (terms == 3) {
        return 1.15e-2;
    } else if (terms == 4) {
        return 3.18e-2;
    } else {
        return 0; // not supported
    }
}

double smoothfun3_full(const double x)
{
    double x2 = x*x;
    double x4 = x2*x2;
    double b = 8.0*x4*x2 - 36.0*x4 + 6.0*x2 + 9.0;
    double c = 2.0/9.0;
    return erf(x) - c*INV_SQRT_PI*x*b*exp(-x2);
}

double smoothfun3_full_taylor(const double d, const double rd, const int terms)
{
    double d2 = d*d;
    double d3 = d2*d;
    double d5 = d2*d3;
    double term = INV_SQRT_PI/d5;
    double s3 = (128.0/15.0)*term;
    if (terms == 1) {
        return s3;
    }
    double rd2 = rd*rd;
    term = term * rd2;
    s3 += -(640.0/63.0)*term;
    if (terms == 2) {
        return s3;
    }
    term = term * rd2;
    s3 += (160.0/27.0)*term;
    if (terms == 3) {
        return s3;
    }
    term = term * rd2;
    s3 += -(224.0/99.0)*term;
    if (terms == 4) {
        return s3;
    }
    return 0; // not supported
}

double get_smoothfun3_full_threshold(const int terms)
{
    if (terms == 1) {
        return 2.11e-3;
    } else if (terms == 2) {
        return 1.05e-2;
    } else if (terms == 3) {
        return 2.88e-2;
    } else if (terms == 4) {
        return 5.78e-2;
    } else {
        return 0; // not supported
    }
}

double smoothfun3_simple(const double x)
{
    double x2 = x*x;
    double x4 = x2*x2;
    double b = -4.0*x4 + 6.0*x2 + 9.0;
    double c = 2.0/9.0;
    return erf(x) - c*INV_SQRT_PI*x*b*exp(-x2);
}

double smoothfun3_simple_taylor(const double d, const double rd, const int terms)
{
    double d2 = d*d;
    double d3 = d2*d;
    double d5 = d2*d3;
    double term = INV_SQRT_PI/d5;
    double s3 = (64.0/45.0)*term;
    if (terms == 1) {
        return s3;
    }
    double rd2 = rd*rd;
    term = term * rd2;
    s3 += -(80.0/63.0)*term;
    if (terms == 2) {
        return s3;
    }
    term = term * rd2;
    s3 += (16.0/27.0)*term;
    if (terms == 3) {
        return s3;
    }
    term = term * rd2;
    s3 += -(56.0/297.0)*term;
    if (terms == 4) {
        return s3;
    }
    return 0; // not supported
}

double get_smoothfun3_simple_threshold(const int terms)
{
    if (terms == 1) {
        return 2.98e-3;
    } else if (terms == 2) {
        return 1.40e-2;
    } else if (terms == 3) {
        return 3.69e-2;
    } else if (terms == 4) {
        return 7.20e-2;
    } else {
        return 0; // not supported
    }
}

// vim:set shiftwidth=4 softtabstop=4:
