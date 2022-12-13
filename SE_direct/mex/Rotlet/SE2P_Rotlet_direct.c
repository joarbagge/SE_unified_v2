#include <math.h>
#include "SE_Rotlet_direct.h"
#include "rotlet_decomp.h"

#define __EXP_ARG_MAX 600

void SE2P_Rotlet_direct_real(double* restrict u, const size_t nout,
                             const double* restrict x, const double* restrict t,
                             const size_t N, const ewald_opts* restrict opt)
{
    const ptrdiff_t nbox = opt->layers;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif
#ifdef CUTOFF
    const double rc2 = opt->rc * opt->rc;
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
        double um[] = {0, 0, 0};
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

        double r[3], rsh[3], tn[3];
        for (size_t n=0; n<N; n++)                          // for all particles
        {
            r[0] = xm[0] - x[n    ];
            r[1] = xm[1] - x[n+N  ];
            r[2] = xm[2] - x[n+2*N];
            tn[0] = t[n    ];
            tn[1] = t[n+N  ];
            tn[2] = t[n+2*N];

            for (ptrdiff_t j0=-nbox; j0<=nbox; j0++)              // image boxes
                for (ptrdiff_t j1=-nbox; j1<=nbox; j1++)
                {
#ifndef EXTERNAL
                    if (n == m_idx && j0 == 0 && j1 == 0)
                        continue;                                   // skip self
#endif
                    // EXTERNAL: Assuming that r != 0 in home box

                    rsh[0] = r[0] + j0*opt->box[0];
                    rsh[1] = r[1] + j1*opt->box[1];
                    rsh[2] = r[2];
#ifdef CUTOFF
                    if (dot(rsh,rsh) >= rc2) continue;        // skip outside rc
#endif
                    real_part_rotlet_contribution(um, rsh, tn, opt->xi);
                }
        }
        u[m       ] = um[0];
        u[m+nout  ] = um[1];
        u[m+2*nout] = um[2];
    }
}

static inline double theta_plus(double z, double k, double xi)
{
    /* idea for a more stable form [LK] */
    /* exp( k*z + log( erfc(k/(2.0*xi) + xi*z) ) ); */
    double exp_arg = k*z;
    if (exp_arg < __EXP_ARG_MAX)
        return exp(exp_arg)*erfc(k/(2.0*xi) + xi*z);
    else
        return 0.0;
}

static inline double theta_minus(double z, double k, double xi)
{
    /* exp(-k*z + log( erfc(k/(2.0*xi) - xi*z) ) ); */
    double exp_arg = -k*z;
    if (exp_arg < __EXP_ARG_MAX)
        return exp(exp_arg)*erfc(k/(2.0*xi) - xi*z);
    else
        return 0.0;
}

void SE2P_Rotlet_direct_fourier(double* restrict u, const size_t nout,
                                const double* restrict x, const double* restrict t,
                                const size_t N, const ewald_opts* restrict opt)
{
    const ptrdiff_t kmax = opt->layers;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif
    const double kfac[] = {2*PI/opt->box[0], 2*PI/opt->box[1]};
    const double invvol = 1.0/(2*PI*opt->box[0]*opt->box[1]);
    const double xi = opt->xi;
    const double PI2 = PI*PI;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
        double um[] = {0, 0, 0};
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

        double k0, k1, nk, k_dot_r, z, thp, thm, tmp, ckr, skr;
        double Q01, Q02, Q12;
        double t0, t1, t2;
        double Qt_re[3], Qt_im[3];
        for (ptrdiff_t j0=-kmax; j0<=kmax; j0++)             // for k-space cube
            for (ptrdiff_t j1=-kmax; j1<=kmax; j1++)
            {
                if (j0 == 0 && j1 == 0)
                    continue;                                     // exclude k=0

                k0 = kfac[0]*j0;
                k1 = kfac[1]*j1;
                nk = sqrt(k0*k0 + k1*k1);

                for (size_t n=0; n<N; n++)                  // for all particles
                {
                    k_dot_r = k0*(xm[0]-x[n]) + k1*(xm[1]-x[n+N]);
                    z = xm[2] - x[n+2*N];
                    thp = theta_plus(z, nk, xi);
                    thm = theta_minus(z, nk, xi);
                    tmp = 2*PI2 * (thm+thp) / nk;

                    /* Fill Q matrix entries (antisymmetric matrix) */
                    Q01 = 2*PI2 * (thm-thp); // REAL
                    Q02 = k1*tmp; // IMAG
                    Q12 = -k0*tmp; // IMAG

                    t0 = t[n    ];
                    t1 = t[n+N  ];
                    t2 = t[n+2*N];

                    /* Real part of Q*t */
                    Qt_re[0] = Q01*t1;
                    Qt_re[1] = -Q01*t0;
                    Qt_re[2] = 0;

                    /* Imag part of Q*t */
                    Qt_im[0] = Q02*t2;
                    Qt_im[1] = Q12*t2;
                    Qt_im[2] = -Q02*t0 - Q12*t1;

                    /*
                       Sum contributions --
                       imag part enters with negative sign (from i*i)
                    */
                    ckr = cos(k_dot_r);
                    skr = sin(k_dot_r);
                    um[0] += ckr*Qt_re[0] - skr*Qt_im[0];
                    um[1] += ckr*Qt_re[1] - skr*Qt_im[1];
                    um[2] += ckr*Qt_re[2] - skr*Qt_im[2];
                }
            }
        u[m       ] = invvol*um[0];
        u[m+nout  ] = invvol*um[1];
        u[m+2*nout] = invvol*um[2];
    }
}

void SE2P_Rotlet_direct_fourier_k0(double* restrict u, const size_t nout,
                                   const double* restrict x, const double* restrict t,
                                   const size_t N, const ewald_opts* restrict opt)
{
    /* Zero-mode Fourier-part rotlet */
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif
    const double prefac = 2*PI/(opt->box[0]*opt->box[1]);
    const double xi = opt->xi;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
        double um[] = {0, 0, 0};
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

        double z, erfval, t0, t1;
        for (size_t n=0; n<N; n++)                          // for all particles
        {
            z = xm[2] - x[n+2*N];
            erfval = erf(xi*z);

            t0 = t[n  ];
            t1 = t[n+N];

            um[0] += t1 * erfval;
            um[1] += -t0 * erfval;
            //um[2] += 0;
        }
        u[m       ] = prefac*um[0];
        u[m+nout  ] = prefac*um[1];
        u[m+2*nout] = 0;
    }
}
