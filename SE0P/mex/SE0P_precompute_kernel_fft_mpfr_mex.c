// MPFR version -- this code uses extended precision to compute
// the kernel more accurately. However, it is very slow compared
// to the standard code `SE0P_precompute_kernel_fft_mex.c`.

/* Compute the modified Green's function FFT (on upsampled grid):
 *
 * HARMONIC_KERNEL:
 *   GR = 4*pi*( 1 - cos(R*k) )/k^2
 *   GR(k==0) = 2*pi*R^2
 *   (which regularizes 4*pi/k^2)
 *
 * BIHARMONIC_KERNEL:
 *   GR = 4*pi*( (2 - (aB+R)*R*k^2 - bB*R*(-6+R^2*k^2))*cos(R*k) + (aB+2*R)*k*sin(R*k) + 3*bB*R*(-2+R^2*k^2)*sin(R*k)/(R*k) - 2 )/k^4
 *   GR(k==0) = pi*R^3*((4/3)*aB + R + (4/5)*bB*R^2)
 *   (which regularizes -8*pi/k^4)
 */

/* Equivalent MATLAB code (requires much more memory):
R = opt.greens_truncation_R;
[k1, k2, k3] = SE_k_vectors(opt.upsampled_grid, opt.upsampled_box, 'shifted');
[K1, K2, K3] = ndgrid(k1, k2, k3);
k = sqrt(K1.^2 + K2.^2 + K3.^2);
GR = <formula for GR>
GR(k==0) = <formula for GR(k==0)>
*/

#include "mex.h"
#include <math.h>
#include <mpfr.h>

#define GRID    prhs[0] // Grid size (of upsampled grid) -- vector of size 3
#define BOX     prhs[1] // Box side lengths (of upsampled box) -- vector of size 3
#define R_TRUNC prhs[2] // Green's truncation length -- scalar
#define PREC    prhs[3] // Variable precision (bits) -- scalar

#define OUT     plhs[0] // Output -- modified Green's function

#define SHIFTED_K_VECTORS 1
#define BIHARMONIC_VARIANT 3 // 3 is the best

void mexFunction(const int, mxArray *[], const int, const mxArray *[]);
void k_vectors(const mwSize, const mpfr_t, mpfr_t *, mpfr_t *, const mpfr_t, const mpfr_rnd_t);


void mexFunction(const int nlhs,       mxArray *plhs[],
                 const int nrhs, const mxArray *prhs[])
{
    // Verify input dimensions
    if (nrhs != 4)
        mexErrMsgIdAndTxt("SE0PMex:WrongNumberOfInputs",
                          "SE0P_precompute_kernel_fft_mpfr_mex takes exactly 4 input arguments");
    if (mxGetNumberOfElements(GRID) != 3)
        mexErrMsgIdAndTxt("SE0PMex:WrongDimensions", "Input GRID must have exactly 3 elements");
    if (mxGetNumberOfElements(BOX) != 3)
        mexErrMsgIdAndTxt("SE0PMex:WrongDimensions", "Input BOX must have exactly 3 elements");
    if (mxGetNumberOfElements(R_TRUNC) != 1)
        mexErrMsgIdAndTxt("SE0PMex:WrongDimensions", "Input R_TRUNC must have exactly 1 element");
    if (mxGetNumberOfElements(PREC) != 1)
        mexErrMsgIdAndTxt("SE0PMex:WrongDimensions", "Input PREC must have exactly 1 element");

    int precision = (int)*mxGetPr(PREC);
    mpfr_rnd_t round_type = MPFR_RNDN;

    mpfr_t pi, R, RR, tmp;
    mpfr_inits2(precision, pi, R, RR, tmp, (mpfr_ptr) 0);
    mpfr_const_pi(pi, round_type);

    const double* Grid = mxGetPr(GRID);
    const double* box = mxGetPr(BOX);
    mpfr_t boxL[3];
    for (int j=0; j<3; j++)
    {
        mpfr_init2(boxL[j], precision);
        mpfr_set_d(boxL[j], box[j], round_type);
    }
    mpfr_set_d(R, mxGetPr(R_TRUNC)[0], round_type);
    mpfr_mul(RR, R, R, round_type);
    const mwSize grid[3] = {Grid[0], Grid[1], Grid[2]};

    // Allocate output
    OUT = mxCreateNumericArray(3, grid, mxDOUBLE_CLASS, mxREAL);
    double* restrict out = mxGetPr(OUT);

    // Compute k vectors
    mpfr_t* restrict k0 = mxMalloc(grid[0]*sizeof(mpfr_t));
    for (mwSize j=0; j<grid[0]; j++)
    {
        mpfr_init2(k0[j], precision);
    }
    mpfr_t* restrict k1 = mxMalloc(grid[1]*sizeof(mpfr_t));
    for (mwSize j=0; j<grid[1]; j++)
    {
        mpfr_init2(k1[j], precision);
    }
    mpfr_t* restrict k2 = mxMalloc(grid[2]*sizeof(mpfr_t));
    for (mwSize j=0; j<grid[2]; j++)
    {
        mpfr_init2(k2[j], precision);
    }
    k_vectors(grid[0], boxL[0], k0, &tmp, pi, round_type);
    k_vectors(grid[1], boxL[1], k1, &tmp, pi, round_type);
    k_vectors(grid[2], boxL[2], k2, &tmp, pi, round_type);

    // Compute modified Green's function
    mpfr_t kk, k, Rk, grval;
    mpfr_inits2(precision, kk, k, Rk, grval, (mpfr_ptr) 0);
    mwSize a = grid[0], b = grid[0]*grid[1];
    for (mwSize n2=0; n2<grid[2]; n2++)
        for (mwSize n1=0; n1<grid[1]; n1++)
            for (mwSize n0=0; n0<grid[0]; n0++)
            {
                // Compute kk = k0[n0]*k0[n0] + k1[n1]*k1[n1] + k2[n2]*k2[n2]
                mpfr_mul(kk, k0[n0], k0[n0], round_type);
                mpfr_mul(tmp, k1[n1], k1[n1], round_type);
                mpfr_add(kk, kk, tmp, round_type);
                mpfr_mul(tmp, k2[n2], k2[n2], round_type);
                mpfr_add(kk, kk, tmp, round_type);
                // Check kk < __DBL_EPSILON__
                if (mpfr_cmp_d(kk, __DBL_EPSILON__) < 0)
                {
#ifdef HARMONIC_KERNEL
                    // Compute grval = 2*pi*RR
                    mpfr_mul(grval, pi, RR, round_type);
                    mpfr_mul_si(grval, grval, 2, round_type);
#elif BIHARMONIC_KERNEL
#if BIHARMONIC_VARIANT == 1 /* aB==0, bB==0 */
                    // Compute grval = pi*RR*RR
                    mpfr_mul(grval, RR, RR, round_type);
                    mpfr_mul(grval, grval, pi, round_type);
#elif BIHARMONIC_VARIANT == 2 /* aB==-R, bB==0 */
                    // Compute grval = -pi*RR*RR/3
                    mpfr_mul(grval, RR, RR, round_type);
                    mpfr_mul(grval, grval, pi, round_type);
                    mpfr_div_si(grval, grval, -3, round_type);
#elif BIHARMONIC_VARIANT == 3 /* aB==-0.5*R, bB==-0.5/R */
                    // Compute grval = -pi*RR*RR/15
                    mpfr_mul(grval, RR, RR, round_type);
                    mpfr_mul(grval, grval, pi, round_type);
                    mpfr_div_si(grval, grval, -15, round_type);
#else
#error "BIHARMONIC_VARIANT must be 1, 2 or 3"
#endif /* BIHARMONIC_VARIANT */
#else
#error "Must give -D<HARMONIC/BIHARMONIC>_KERNEL to compiler"
#endif
                }
                else
                {
                    // Compute k = sqrt(kk)
                    mpfr_sqrt(k, kk, round_type);
                    // Compute Rk = R*k
                    mpfr_mul(Rk, R, k, round_type);
                    /* Computing trigonometric functions in GR is
                     * the most time-consuming part of the computation.
                     */
#ifdef HARMONIC_KERNEL
                    // Compute grval = 4*pi*( 1 - cos(Rk) )/kk
                    mpfr_cos(tmp, Rk, round_type);
                    mpfr_si_sub(tmp, 1, tmp, round_type);
                    mpfr_div(tmp, tmp, kk, round_type);
                    mpfr_mul(tmp, tmp, pi, round_type);
                    mpfr_mul_si(grval, tmp, 4, round_type);
#elif BIHARMONIC_KERNEL
#if BIHARMONIC_VARIANT == 1 /* aB==0, bB==0 */
                    // Compute grval = 4*pi*( (2-Rk*Rk)*cos(Rk) + 2*Rk*sin(Rk) - 2 )/(kk*kk)
                    mpfr_cos(grval, Rk, round_type);
                    mpfr_mul(tmp, Rk, Rk, round_type);
                    mpfr_si_sub(tmp, 2, tmp, round_type);
                    mpfr_mul(grval, grval, tmp, round_type);
                    mpfr_sin(tmp, Rk, round_type);
                    mpfr_mul(tmp, tmp, Rk, round_type);
                    mpfr_mul_si(tmp, tmp, 2, round_type);
                    mpfr_add(tmp, tmp, grval, round_type);
                    mpfr_sub_si(tmp, tmp, 2, round_type);
                    mpfr_mul(grval, kk, kk, round_type);
                    mpfr_div(tmp, tmp, grval, round_type);
                    mpfr_mul(tmp, tmp, pi, round_type);
                    mpfr_mul_si(grval, tmp, 4, round_type);
#elif BIHARMONIC_VARIANT == 2 /* aB==-R, bB==0 */
                    // Compute grval = 4*pi*( 2*cos(Rk) + Rk*sin(Rk) - 2 )/(kk*kk)
                    mpfr_cos(grval, Rk, round_type);
                    mpfr_mul_si(grval, grval, 2, round_type);
                    mpfr_sin(tmp, Rk, round_type);
                    mpfr_mul(tmp, tmp, Rk, round_type);
                    mpfr_add(tmp, tmp, grval, round_type);
                    mpfr_sub_si(tmp, tmp, 2, round_type);
                    mpfr_mul(grval, kk, kk, round_type);
                    mpfr_div(tmp, tmp, grval, round_type);
                    mpfr_mul(tmp, tmp, pi, round_type);
                    mpfr_mul_si(grval, tmp, 4, round_type);
#elif BIHARMONIC_VARIANT == 3 /* aB==-0.5*R, bB==-0.5/R */
                    // Compute grval = 4*pi*( -cos(Rk) + 3*sin(Rk)/Rk - 2 )/(kk*kk)
                    mpfr_cos(grval, Rk, round_type);
                    mpfr_mul_si(grval, grval, -1, round_type);
                    mpfr_sin(tmp, Rk, round_type);
                    mpfr_div(tmp, tmp, Rk, round_type);
                    mpfr_mul_si(tmp, tmp, 3, round_type);
                    mpfr_add(tmp, tmp, grval, round_type);
                    mpfr_sub_si(tmp, tmp, 2, round_type);
                    mpfr_mul(grval, kk, kk, round_type);
                    mpfr_div(tmp, tmp, grval, round_type);
                    mpfr_mul(tmp, tmp, pi, round_type);
                    mpfr_mul_si(grval, tmp, 4, round_type);
#else
#error "BIHARMONIC_VARIANT must be 1, 2 or 3"
#endif /* BIHARMONIC_VARIANT */
#else
#error "Must give -D<HARMONIC/BIHARMONIC>_KERNEL to compiler"
#endif
                }
                out[n0 + n1*a + n2*b] = mpfr_get_d(grval, round_type);
            }

    for (int j=0; j<3; j++)
    {
        mpfr_clear(boxL[j]);
    }
    for (mwSize j=0; j<grid[0]; j++)
    {
        mpfr_clear(k0[j]);
    }
    for (mwSize j=0; j<grid[1]; j++)
    {
        mpfr_clear(k1[j]);
    }
    for (mwSize j=0; j<grid[2]; j++)
    {
        mpfr_clear(k2[j]);
    }

    mxFree(k0);
    mxFree(k1);
    mxFree(k2);

    mpfr_clears(pi, R, RR, tmp, kk, k, Rk, grval, (mpfr_ptr) 0);
    mpfr_free_cache();
}

void k_vectors(const mwSize M, const mpfr_t L, mpfr_t *k, mpfr_t *tmp,
               const mpfr_t pi, const mpfr_rnd_t round_type)
{
    mpfr_mul_si(*tmp, pi, 2, round_type);
    mpfr_div(*tmp, *tmp, L, round_type);
    mwSize Mh = M/2;
#if SHIFTED_K_VECTORS
    if (M % 2 != 0)
        Mh = (M-1)/2;
    for (mwSize i=0; i<M; i++)
    {
        mpfr_mul_si(k[i], *tmp, (long int)i-Mh, round_type);
    }
#else
    if (M % 2 != 0)
        Mh = (M+1)/2;
    for (mwSize i=0; i<Mh; i++)
    {
        mpfr_mul_si(k[i], *tmp, i, round_type);
    }
    for (mwSize i=Mh; i<M; i++)
    {
        mpfr_mul_si(k[i], *tmp, (long int)i-M, round_type);
    }
#endif
}
