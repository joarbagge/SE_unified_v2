#include <math.h>
#include "SE_Stokeslet_direct.h"

#ifdef HASIMOTO
#include "hasimoto_decomp.h"
#elif BEENAKKER
#error "Beenakker k-space sum not implemented for 2P!"
#else
#error "Must provide -D<decomposition> to compiler"
#endif

#define __EXP_ARG_MAX 600

void SE2P_Stokeslet_direct_real(double* restrict u, size_t nout,
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
    double rvec[3];
    double r[3];
    double A[3][3];

#ifdef _OPENMP
#pragma omp parallel for private(rvec,r,A)
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
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
                for (ptrdiff_t p1 = -nbox; p1<=nbox; p1++)
                {
#ifndef EXTERNAL
                    if (m_idx == n && p1 == 0 && p0 == 0)           // skip self
                        continue;
#endif
                    // EXTERNAL: Assuming that r != 0 in home box

                    r[0] = rvec[0] + p0*opt->box[0];
                    r[1] = rvec[1] + p1*opt->box[1];
                    r[2] = rvec[2];

#ifdef CUTOFF
                    if (r[0]*r[0] + r[1]*r[1] + r[2]*r[2] >= rc2)
                        continue;                             // skip outside rc
#endif

                    op_A(A, r, opt->xi);                             // u += A*f
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

static inline double theta_plus(double z, double k, double xi)
{
    /* idea for a more stable form [LK] */
    /* exp( k*z + log( erfc(k/(2.0*xi) + xi*z) ) ); */

    if(k*z <  __EXP_ARG_MAX)
        return exp( k*z)*erfc(k/(2.0*xi) + xi*z);
    else
        return 0.0;
}

static inline double theta_minus(double z, double k, double xi)
{
    /* exp(-k*z + log( erfc(k/(2.0*xi) - xi*z) ) ); */

    if(-k*z <  __EXP_ARG_MAX)
        return exp(-k*z)*erfc(k/(2.0*xi) - xi*z);
    else
        return 0.0;
}

static inline double Lambda_plus(double z, double k, double cp, double cm)
{
    return (cm+cp)/k + (cm-cp)*z;
}

static inline double Lambda_minus(double z, double k, double cp, double cm)
{
    return -(cm+cp)/k + (cm-cp)*z;
}

void SE2P_Stokeslet_direct_fourier(double* restrict u, size_t nout,
                                   const double* restrict x, const double* restrict f, size_t N,
                                   const ewald_opts* restrict opt)
{
    double xm[3];
    double k_dot_r, z;
    double f0, f1, f2;
    double k0, k1, nk;
    double Qf_re[3];
    double Qf_im[3];
    double p[3];
    double ce, cp, cm, Lp, Lm, cL;
    const double C = 4/(opt->box[0]*opt->box[1]);
    const double xi = opt->xi;

    const ptrdiff_t kmax = opt->layers;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif

#ifdef _OPENMP
#pragma omp parallel for private(xm,k_dot_r,z,f0,f1,f2,k0,k1,nk,Qf_re,Qf_im,p,ce,cp,cm,Lp,Lm,cL)
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
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

        p[0] = 0;
        p[1] = 0;
        p[2] = 0;

        for (size_t n = 0; n<N; n++)                        // for all particles
        {
            f0 = f[n    ];
            f1 = f[n+N  ];
            f2 = f[n+2*N];
            for (ptrdiff_t j0 = -kmax; j0<=kmax; j0++)        // for wavenumbers
                for (ptrdiff_t j1 = -kmax; j1<=kmax; j1++)
                {
                    if (j0 == 0 && j1 == 0)                 // exclude zero mode
                        continue;

                    k0 = 2*PI*j0/opt->box[0];
                    k1 = 2*PI*j1/opt->box[1];
                    nk = sqrt(k0*k0 + k1*k1);
                    k_dot_r = k0*(xm[0]-x[n]) + k1*(xm[1]-x[n+N]);

                    z = xm[2]-x[n+2*N];

                    ce = exp( -(nk*nk)/(4.0*(xi*xi)) - (xi*xi)*(z*z) );
                    cp = theta_plus(z,nk,xi);
                    cm = theta_minus(z,nk,xi);
                    Lp = Lambda_plus(z,nk,cp,cm);
                    Lm = Lambda_minus(z,nk,cp,cm);
                    cL = 2*sqrt(PI)*ce + PI*xi*Lp;

                    /* Real part of Q*f */
                    Qf_re[0] = cL*k1*(k1*f0 - k0*f1)/(4.*xi*nk*nk) - 0.25*PI*Lm*f0;
                    Qf_re[1] = cL*k0*(k0*f1 - k1*f0)/(4.*xi*nk*nk) - 0.25*PI*Lm*f1;
                    Qf_re[2] = cL*f2/(4.*xi);

                    /* Imag part of Q*f */
                    Qf_im[0] = -(cm + cp)*f2*k0*PI*z/(4.*nk);
                    Qf_im[1] = -(cm + cp)*f2*k1*PI*z/(4.*nk);
                    Qf_im[2] = -(cm + cp)*(f0*k0 + f1*k1)*PI*z/(4.*nk);

                    /*
                       Sum contributions --
                       imag part enters with negative sign (from i*i)
                    */
                    p[0] += cos(k_dot_r)*Qf_re[0] - sin(k_dot_r)*Qf_im[0];
                    p[1] += cos(k_dot_r)*Qf_re[1] - sin(k_dot_r)*Qf_im[1];
                    p[2] += cos(k_dot_r)*Qf_re[2] - sin(k_dot_r)*Qf_im[2];
                }
        }
        u[m       ] = C*p[0];
        u[m+nout  ] = C*p[1];
        u[m+2*nout] = C*p[2];
    }
}

void SE2P_Stokeslet_direct_fourier_k0(double* restrict u, size_t nout,
                                      const double* restrict x, const double* restrict f, size_t N,
                                      const ewald_opts* restrict opt)
{
    double z, zm, p[2], q;
    const double C = -4*sqrt(PI)/(opt->box[0]*opt->box[1]);
    const double xi = opt->xi;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif

#ifdef _OPENMP
#pragma omp parallel for private(z,zm,p,q)
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
#ifdef EXTERNAL
        zm = xt[m+2*nout];
#else
        size_t m_idx;
        if (opt->N_eval == -1)
            m_idx = m;            // default is to evaluate at all source points
        else
            m_idx = idx[m];                // indirect indexing OK in outer loop
        zm = x[m_idx+2*N];
#endif

        p[0] = 0; p[1] = 0;
        for (size_t n=0; n<N; n++)                          // for all particles
        {
            z = zm - x[n+2*N];
            q = exp(-xi*xi*z*z)/(2*xi) + sqrt(PI)*z*erf(xi*z);
            p[0] += q*f[n  ];
            p[1] += q*f[n+N];
        }
        u[m     ] = C*p[0];
        u[m+nout] = C*p[1];
    }
}

void SE2P_Stokeslet_direct_self(double* restrict u, size_t nout,
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
