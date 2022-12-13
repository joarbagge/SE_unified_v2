#include <math.h>
#include "SE_Stokeslet_direct.h"
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include "../incomplete_bessel_K.h"

#ifdef HASIMOTO
#include "hasimoto_decomp.h"
#elif BEENAKKER
#error "Beenakker k-space sum not implemented for 1P!"
#else
#error "Must provide -D<decomposition> to compiler"
#endif

#ifdef EWALD_FULL
inline double dot(double * a, double * b)
{
    return ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

void SE1P_Stokeslet_direct_full(double* restrict u, size_t nout,
                                const double* restrict x, const double* restrict f, size_t N,
                                const ewald_opts* restrict opt)
{
    const ptrdiff_t nbox = opt->layers;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
#ifdef EXTERNAL
        double xm[] = {xt[m], xt[m+nout], xt[m+2*nout]};
#else
        size_t m_idx;
        if (opt->N_eval == -1)
            m_idx = m;            // default is to evaluate at all source points
        else
            m_idx = idx[m];                // indirect indexing OK in outer loop
        double xm[] = {x[m_idx], x[m_idx+N], x[m_idx+2*N]};
#endif

        double rvec[3];
        double r[3];
        double p[] = {0,0,0};
        for (size_t n=0; n<N; n++)                          // for all particles
        {
            double fn[] = {f[n], f[n+N], f[n+2*N]};
            rvec[0] = xm[0] - x[n    ];
            rvec[1] = xm[1] - x[n+N  ];
            rvec[2] = xm[2] - x[n+2*N];

            for (ptrdiff_t p0 = -nbox; p0<=nbox; p0++)            // image boxes
            {
#ifndef EXTERNAL
                if (m_idx == n && p0 == 0)                          // skip self
                    continue;
#endif
                // EXTERNAL: Assuming that r != 0 in home box

                r[0] = rvec[0] + p0*opt->box[0];
                r[1] = rvec[1];
                r[2] = rvec[2];

                double rdotf = dot(r,fn);
                double ri = 1.0/sqrt(dot(r,r));
                double ri3 = ri*ri*ri;

                p[0] += ri*fn[0] + rdotf*ri3*r[0];
                p[1] += ri*fn[1] + rdotf*ri3*r[1];
                p[2] += ri*fn[2] + rdotf*ri3*r[2];
            }

            u[m       ] = p[0];
            u[m+nout  ] = p[1];
            u[m+2*nout] = p[2];
        }
    }
}
#endif // EWALD_FULL

#ifdef EWALD_REAL
void SE1P_Stokeslet_direct_real(double* restrict u, size_t nout,
                                const double* restrict x, const double* restrict f, size_t N,
                                const ewald_opts* restrict opt)
{
#ifdef CUTOFF
    const double rc2 = opt->rc * opt->rc;
#endif
    const ptrdiff_t nbox = opt->layers;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
        double rvec[3];
        double r[3];
        double A[3][3];
        u[m       ] = 0;
        u[m+nout  ] = 0;
        u[m+2*nout] = 0;

#ifdef EXTERNAL
        double xm[] = {xt[m], xt[m+nout], xt[m+2*nout]};
#else
        size_t m_idx;
        if (opt->N_eval == -1)
            m_idx = m;            // default is to evaluate at all source points
        else
            m_idx = idx[m];                // indirect indexing OK in outer loop
        double xm[] = {x[m_idx], x[m_idx+N], x[m_idx+2*N]};
#endif

        for (size_t n=0; n<N; n++)                          // for all particles
        {
            rvec[0] = xm[0] - x[n    ];
            rvec[1] = xm[1] - x[n+N  ];
            rvec[2] = xm[2] - x[n+2*N];

            for (ptrdiff_t p0 = -nbox; p0<=nbox; p0++)            // image boxes
            {
#ifndef EXTERNAL
                if (m_idx == n && p0 == 0)                          // skip self
                    continue;
#endif
                // EXTERNAL: Assuming that r != 0 in home box

                r[0] = rvec[0] + p0*opt->box[0];
                r[1] = rvec[1];
                r[2] = rvec[2];

#ifdef CUTOFF
                if (r[0]*r[0] + r[1]*r[1] + r[2]*r[2] >= rc2)
                    continue;                                 // skip outside rc
#endif

                op_A(A, r, opt->xi);                                 // u += A*f
                u[m       ] +=
                    A[0][0]*f[n    ]+
                    A[0][1]*f[n+N  ]+
                    A[0][2]*f[n+2*N];
                u[m+nout  ] +=
                    A[1][0]*f[n    ]+
                    A[1][1]*f[n+N  ]+
                    A[1][2]*f[n+2*N];
                u[m+2*nout] +=
                    A[2][0]*f[n    ]+
                    A[2][1]*f[n+N  ]+
                    A[2][2]*f[n+2*N];
            }
        }
    }
}
#endif // EWALD_REAL

#ifdef EWALD_FOURIER
void SE1P_Stokeslet_direct_fourier(double* restrict u, size_t nout,
                                   const double* restrict x, const double* restrict f, size_t N,
                                   const ewald_opts* restrict opt)
{
    const double xi = opt->xi;
    const double xi2 = xi*xi;
    const double D = 2.0/(PI*opt->box[0]);
    double p[3], xm[3], r[3];
    double pif0, pif1, pif2; // PI times force (for convenience)
    double k0, k0r0;
    double Qf_re[3], Qf_im[3];

    const ptrdiff_t kmax = opt->layers;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif

#ifdef _OPENMP
#pragma omp parallel for private(p,xm,r,pif0,pif1,pif2,k0,k0r0,Qf_re,Qf_im)
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
        double bessel_K0, bessel_K1, bessel_Km1, U, V;
        double A, B, C;
#ifdef EXTERNAL
        xm[0] = xt[m       ];
        xm[1] = xt[m+nout  ];
        xm[2] = xt[m+2*nout];
#else
        size_t m_idx;
        if (opt->N_eval == -1)
            m_idx = m;            // default is to evaluate at all source points
        else
            m_idx = idx[m];                // indirect indexing OK in outer loop
        xm[0] = x[m_idx    ];
        xm[1] = x[m_idx+N  ];
        xm[2] = x[m_idx+2*N];
#endif

        p[0] = 0; p[1] = 0; p[2] = 0;

        for (size_t n = 0; n<N; n++)                        // for all particles
        {
            r[0] = xm[0] - x[n    ];
            r[1] = xm[1] - x[n+N  ];
            r[2] = xm[2] - x[n+2*N];
            pif0 = PI*f[n    ];
            pif1 = PI*f[n+N  ];
            pif2 = PI*f[n+2*N];
            for (ptrdiff_t j=-kmax; j<=kmax; j++)             // for wavenumbers
            {
                if (j == 0) continue;                       // exclude zero mode

                k0 = 2*PI*j/opt->box[0];
                k0r0 = k0*r[0];

                /* Compute incomplete modified Bessel functions */
                U = 0.25*k0*k0/xi2;
                V = (r[1]*r[1]+r[2]*r[2])*xi2;
                bessel_K0 = incomplete_bessel_K(0, U, V);
                bessel_K1 = incomplete_bessel_K(1, U, V);
                bessel_Km1 = (exp(-U-V) + V*bessel_K1)/U;

                A = -0.5*k0*bessel_K0;
                B = r[1]*r[2]*xi2*bessel_K1;
                C = 0.5*bessel_K0 + U*bessel_Km1;

                /* Real part of Q*f */
                Qf_re[0] = pif0*(bessel_K0 - V*bessel_K1);
                Qf_re[1] = pif1*(C - r[2]*r[2]*xi2*bessel_K1) + pif2*B;
                Qf_re[2] = pif2*(C - r[1]*r[1]*xi2*bessel_K1) + pif1*B;

                /* Imag part of Q*f */
                Qf_im[0] = A*(r[1]*pif1 + r[2]*pif2);
                Qf_im[1] = A*r[1]*pif0;
                Qf_im[2] = A*r[2]*pif0;

                /* Sum contributions --
                 * imag part enters with negative sign (from i*i)
                 */
                p[0] += cos(k0r0)*Qf_re[0] - sin(k0r0)*Qf_im[0];
                p[1] += cos(k0r0)*Qf_re[1] - sin(k0r0)*Qf_im[1];
                p[2] += cos(k0r0)*Qf_re[2] - sin(k0r0)*Qf_im[2];
            }
        }
        u[m       ] = D*p[0];
        u[m+nout  ] = D*p[1];
        u[m+2*nout] = D*p[2];
    }
}
#endif // EWALD_FOURIER

#ifdef EWALD_FOURIER_K0
void SE1P_Stokeslet_direct_fourier_k0(double* restrict u, size_t nout,
                                      const double* restrict x, const double* restrict f, size_t N,
                                      const ewald_opts* restrict opt)
{
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif
    const double bi_c = opt->stokeslet_k0_constant;
    const double bi_c1 = 2*bi_c - 2;
    const double bi_c2 = 2*bi_c - 3;
    const double euler_gamma = 0.57721566490153286061;
    const double xi = opt->xi;
    const double D = 1/opt->box[0];
    double p[3], y, ym, z, zm, r2;
    double v, A, B, C;
    gsl_set_error_handler_off();
    int gsl_error_code = GSL_SUCCESS;

#ifdef _OPENMP
#pragma omp parallel for private(p,y,ym,z,zm,r2,v,A,B,C)
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
#ifdef EXTERNAL
        ym = xt[m+nout  ];
        zm = xt[m+2*nout];
#else
        size_t m_idx;
        if (opt->N_eval == -1)
            m_idx = m;            // default is to evaluate at all source points
        else
            m_idx = idx[m];                // indirect indexing OK in outer loop
        ym = x[m_idx+N  ];
        zm = x[m_idx+2*N];
#endif

        int status;
        gsl_sf_result result;
        p[0] = 0; p[1] = 0; p[2] = 0;
        for (size_t n=0; n<N; n++)                          // for all particles
        {
            y = ym - x[n+N  ];
            z = zm - x[n+2*N];
            r2 = y*y + z*z;
            v = r2*xi*xi;
            A = exp(-v);
            if (r2 > __DBL_EPSILON__)
            {
                status = gsl_sf_expint_E1_e(v, &result);
                if (status == GSL_EUNDRFLW) // underflow
                {
                    printf("SE:DirectGSL Warning: GSL underflow occurred.\n");
                }
                else if (status != GSL_SUCCESS)
                {
#ifdef _OPENMP
#pragma omp critical
#endif
                    if (gsl_error_code == GSL_SUCCESS) {
                        gsl_error_code = status;
                    }
                    break;
                }
                B = A - log(r2) - result.val;
                C = (2.0/r2)*(y*f[n+N] + z*f[n+2*N])*(1-A);
            }
            else
            {
                /* Limit for r2=0 */
                B = 1 + euler_gamma + 2*log(xi);
                C = 0;
            }
            p[0] += (B+bi_c1)*f[n];
            // bi_c1 is to account for the difference from bi_c=1
            p[1] += (A+B+bi_c2)*f[n+N  ] + C*y;
            p[2] += (A+B+bi_c2)*f[n+2*N] + C*z;
            // bi_c2 is to account for the difference from bi_c=3/2
        }
        u[m       ] = 2*D*p[0];
        u[m+nout  ] = D*p[1];
        u[m+2*nout] = D*p[2];
    }

    if (gsl_error_code != GSL_SUCCESS)
    {
        char errMsg [64];
        sprintf(errMsg, "GSL error code: %d", gsl_error_code);
        mexErrMsgIdAndTxt("SE:DirectGSL", errMsg);
    }
}
#endif // EWALD_FOURIER_K0

#ifdef EWALD_SELF
void SE1P_Stokeslet_direct_self(double* restrict u, size_t nout,
                                const double* restrict f, size_t N,
                                const ewald_opts* restrict opt)
{
    const size_t* restrict idx = opt->eval_idx;
    const double c = -4*opt->xi/sqrt(PI);
    size_t m_idx;
    for (size_t m=0; m<nout; m++)
    {
        if (opt->N_eval == -1)
            m_idx = m;            // default is to evaluate at all source points
        else
            m_idx = idx[m];

        u[m       ] = c*f[m_idx    ];
        u[m+nout  ] = c*f[m_idx+N  ];
        u[m+2*nout] = c*f[m_idx+2*N];
    }
}
#endif // EWALD_SELF
