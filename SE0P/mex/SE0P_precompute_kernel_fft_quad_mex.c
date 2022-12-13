// libquadmath version

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
#include <quadmath.h>

#define GRID    prhs[0] // Grid size (of upsampled grid) -- vector of size 3
#define BOX     prhs[1] // Box side lengths (of upsampled box) -- vector of size 3
#define R_TRUNC prhs[2] // Green's truncation length -- scalar

#define OUT     plhs[0] // Output -- modified Green's function

#define SHIFTED_K_VECTORS 1
#define PI M_PIq

#ifdef HARMONIC_KERNEL
#define GR 4*PI*( 1 - cosq(Rk) )/kk
#define GR_k0 2*PI*RR
#elif BIHARMONIC_KERNEL
// aB==0, bB==0
//#define GR 4*PI*( (2-Rk*Rk)*cosq(Rk) + 2*Rk*sinq(Rk) - 2 )/(kk*kk)
//#define GR_k0 PI*RR*RR
// aB==-R, bB==0
//#define GR 4*PI*( 2*cosq(Rk) + Rk*sinq(Rk) - 2 )/(kk*kk)
//#define GR_k0 -PI*RR*RR/3
// aB==-0.5*R, bB==-0.5/R -- the best
#define GR 4*PI*( -cosq(Rk) + 3*sinq(Rk)/Rk - 2 )/(kk*kk)
#define GR_k0 -PI*RR*RR/15
#else
#error "Must give -D<HARMONIC/BIHARMONIC>_KERNEL to compiler"
#endif

void mexFunction(const int, mxArray *[], const int, const mxArray *[]);
void k_vectors(const mwSize, const __float128, __float128 *);


void mexFunction(const int nlhs,       mxArray *plhs[],
                 const int nrhs, const mxArray *prhs[])
{
    // Verify input dimensions
    if (nrhs != 3)
        mexErrMsgIdAndTxt("SE0PMex:WrongNumberOfInputs",
                          "SE0P_precompute_kernel_fft_mex takes exactly 3 input arguments");
    if (mxGetNumberOfElements(GRID) != 3)
        mexErrMsgIdAndTxt("SE0PMex:WrongDimensions", "Input GRID must have exactly 3 elements");
    if (mxGetNumberOfElements(BOX) != 3)
        mexErrMsgIdAndTxt("SE0PMex:WrongDimensions", "Input BOX must have exactly 3 elements");
    if (mxGetNumberOfElements(R_TRUNC) != 1)
        mexErrMsgIdAndTxt("SE0PMex:WrongDimensions", "Input R_TRUNC must have exactly 1 element");

    const double* Grid = mxGetPr(GRID);
    const double* box = mxGetPr(BOX);
    __float128 boxL[3];
    for (int j=0; j<3; j++)
    {
        boxL[j] = box[j];
    }
    const __float128 R = mxGetPr(R_TRUNC)[0];
    const __float128 RR = R*R;
    const mwSize grid[3] = {Grid[0], Grid[1], Grid[2]};

    // Allocate output
    OUT = mxCreateNumericArray(3, grid, mxDOUBLE_CLASS, mxREAL);
    double* restrict out = mxGetPr(OUT);

    // Compute k vectors
    __float128* restrict k0 = mxMalloc(grid[0]*sizeof(__float128));
    __float128* restrict k1 = mxMalloc(grid[1]*sizeof(__float128));
    __float128* restrict k2 = mxMalloc(grid[2]*sizeof(__float128));
    k_vectors(grid[0], boxL[0], k0);
    k_vectors(grid[1], boxL[1], k1);
    k_vectors(grid[2], boxL[2], k2);

    // Compute modified Green's function
    __float128 kk, k, Rk, grval;
    mwSize a = grid[0], b = grid[0]*grid[1];
    for (mwSize n2=0; n2<grid[2]; n2++)
        for (mwSize n1=0; n1<grid[1]; n1++)
            for (mwSize n0=0; n0<grid[0]; n0++)
            {
                kk = k0[n0]*k0[n0] + k1[n1]*k1[n1] + k2[n2]*k2[n2];
                if (kk < __DBL_EPSILON__)
                {
                    grval = GR_k0;
                }
                else
                {
                    k = sqrtq(kk);
                    Rk = R*k;
                    /* Computing trigonometric functions in GR is
                     * the most time-consuming part of the computation.
                     */
                    grval = GR;
                }
                out[n0 + n1*a + n2*b] = (double)grval;
            }

    mxFree(k0);
    mxFree(k1);
    mxFree(k2);
}

void k_vectors(const mwSize M, const __float128 L, __float128 *k)
{
    __float128 c = 2*PI/L;
    mwSize Mh = M/2;
#if SHIFTED_K_VECTORS
    if (M % 2 != 0)
        Mh = (M-1)/2;
    for (mwSize i=0; i<M; i++)
        k[i] = c*((__float128)i-Mh);
#else
    if (M % 2 != 0)
        Mh = (M+1)/2;
    for (mwSize i=0; i<Mh; i++)
        k[i] = c*i;
    for (mwSize i=Mh; i<M; i++)
        k[i] = c*((__float128)i-M);
#endif
}
