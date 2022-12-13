#include <math.h>
#include "SE_Stresslet_direct.h"
#include "hasimoto_decomp.h"

#define __EXP_ARG_MAX 600

void SE2P_Stresslet_direct_real(double* restrict u, const size_t nout,
                                const double* restrict x,
                                const double* restrict q, const double* restrict no,
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

        double r[3], rsh[3], qn[3], non[3];
        for (size_t n=0; n<N; n++)                          // for all particles
        {
            r[0] = xm[0] - x[n    ];
            r[1] = xm[1] - x[n+N  ];
            r[2] = xm[2] - x[n+2*N];
            qn[0]  =  q[n    ];
            qn[1]  =  q[n+N  ];
            qn[2]  =  q[n+2*N];
            non[0] = no[n    ];
            non[1] = no[n+N  ];
            non[2] = no[n+2*N];

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
                    real_part_stresslet_contribution(um, rsh, qn, non, opt->xi);
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

static inline double lambda(double z, double k, double xi)
{
    return exp( -k*k/(4.0*xi*xi) - xi*xi*z*z );
}

void SE2P_Stresslet_direct_fourier(double* restrict u, const size_t nout,
                                   const double* restrict x,
                                   const double* restrict q, const double* restrict no,
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
    const double PI32 = PI*sqrt(PI);
    const double PI32oxi = PI32/xi;

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

        double k0, k1, nk2, nk, onk, K0, K1, k_dot_r, z, thp, thm, bp, bm, lam, ckr, skr;
        double Q000, Q111, Q222, Q001, Q002, Q011, Q012, Q022, Q112, Q122;
        double q0, q1, q2, n0, n1, n2;
        double Qf_re[3], Qf_im[3];
        for (ptrdiff_t j0=-kmax; j0<=kmax; j0++)             // for k-space cube
            for (ptrdiff_t j1=-kmax; j1<=kmax; j1++)
            {
                if (j0 == 0 && j1 == 0)
                    continue;                                 // exclude k=0

                k0 = kfac[0]*j0;
                k1 = kfac[1]*j1;
                nk2 = k0*k0 + k1*k1;
                nk = sqrt(nk2);
                onk = 1.0/nk;
                K0 = k0*onk;
                K1 = k1*onk;

                for (size_t n=0; n<N; n++)              // for all particles
                {
                    k_dot_r = k0*(xm[0]-x[n]) + k1*(xm[1]-x[n+N]);
                    z = xm[2] - x[n+2*N];
                    thp = theta_plus(z, nk, xi);
                    thm = theta_minus(z, nk, xi);
                    bp = thp + thm;
                    bm = thp - thm;
                    lam = lambda(z, nk, xi);

                    /* Fill tensor entries (symmetric tensor) */
                    Q000 = 4*k0*(PI32oxi*lam*(1+2*K1*K1) + PI2*(onk*(2+K1*K1)*bp + K0*K0*bm*z)); // IMAG
                    Q111 = 4*k1*(PI32oxi*lam*(1+2*K0*K0) + PI2*(onk*(2+K0*K0)*bp + K1*K1*bm*z)); // IMAG
                    Q222 = -8*PI32*lam*xi*z - 4*PI2*(nk*bp*z - bm); // REAL
                    Q001 = -4*k1*(PI32oxi*lam*(1-2*K1*K1) - PI2*(onk*K1*K1*bp + K0*K0*bm*z)); // IMAG
                    Q002 = -8*PI32*lam*xi*z + 4*PI2*(nk*K0*K0*bp*z + bm); // REAL
                    Q011 = -4*k0*(PI32oxi*lam*(1-2*K0*K0) - PI2*(onk*K0*K0*bp + K1*K1*bm*z)); // IMAG
                    Q012 = 4*PI2*nk*K0*K1*bp*z; // REAL
                    Q022 = 4*k0*(PI32oxi*lam - PI2*bm*z); // IMAG
                    Q112 = -8*PI32*lam*xi*z + 4*PI2*(nk*K1*K1*bp*z + bm); // REAL
                    Q122 = 4*k1*(PI32oxi*lam - PI2*bm*z); // IMAG

                    q0 = q[n    ];
                    q1 = q[n+N  ];
                    q2 = q[n+2*N];
                    n0 = no[n    ];
                    n1 = no[n+N  ];
                    n2 = no[n+2*N];

                    /* Real part of Q*f */
                    Qf_re[0] =                           Q002*q0*n2
                                                       + Q012*q1*n2
                             + Q002*q2*n0 + Q012*q2*n1             ;
                    Qf_re[1] =                           Q012*q0*n2
                                                       + Q112*q1*n2
                             + Q012*q2*n0 + Q112*q2*n1             ;
                    Qf_re[2] = Q002*q0*n0 + Q012*q0*n1
                             + Q012*q1*n0 + Q112*q1*n1
                                                       + Q222*q2*n2;

                    /* Imag part of Q*f */
                    Qf_im[0] = Q000*q0*n0 + Q001*q0*n1
                             + Q001*q1*n0 + Q011*q1*n1
                                                       + Q022*q2*n2;
                    Qf_im[1] = Q001*q0*n0 + Q011*q0*n1
                             + Q011*q1*n0 + Q111*q1*n1
                                                       + Q122*q2*n2;
                    Qf_im[2] =                           Q022*q0*n2
                                                       + Q122*q1*n2
                             + Q022*q2*n0 + Q122*q2*n1             ;

                    /*
                       Sum contributions --
                       imag part enters with negative sign (from i*i)
                    */
                    ckr = cos(k_dot_r);
                    skr = sin(k_dot_r);
                    um[0] += ckr*Qf_re[0] - skr*Qf_im[0];
                    um[1] += ckr*Qf_re[1] - skr*Qf_im[1];
                    um[2] += ckr*Qf_re[2] - skr*Qf_im[2];
                }
            }
        u[m       ] = invvol*um[0];
        u[m+nout  ] = invvol*um[1];
        u[m+2*nout] = invvol*um[2];
    }
}

void SE2P_Stresslet_direct_fourier_k0(double* restrict u, const size_t nout,
                                      const double* restrict x,
                                      const double* restrict q, const double* restrict no,
                                      const size_t N, const ewald_opts* restrict opt)
{
    /* Zero-mode Fourier-part stresslet */
    /* Computes (1/V) sum_n TF(0)_jlm q_l n_m */
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif
    const double invvol = 1.0/(2*PI*opt->box[0]*opt->box[1]);
    const double xi = opt->xi;
    const double PI2 = PI*PI;
    const double PI32 = PI*sqrt(PI);

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

        double z, Q;
        double q0, q1, q2, n0, n1, n2;
        for (size_t n=0; n<N; n++)              // for all particles
        {
            z = xm[2] - x[n+2*N];

            /* Limit of non-zero elements of Q tensor */
            Q = -8*PI32*exp(-xi*xi*z*z)*xi*z - 8*PI2*erf(xi*z);

            q0 = q[n    ];
            q1 = q[n+N  ];
            q2 = q[n+2*N];
            n0 = no[n    ];
            n1 = no[n+N  ];
            n2 = no[n+2*N];

            um[0] += Q*(q0*n2 + q2*n0);
            um[1] += Q*(q1*n2 + q2*n1);
            um[2] += Q*(q0*n0 + q1*n1 + q2*n2);

        }
        u[m       ] = invvol*um[0];
        u[m+nout  ] = invvol*um[1];
        u[m+2*nout] = invvol*um[2];
    }
}
