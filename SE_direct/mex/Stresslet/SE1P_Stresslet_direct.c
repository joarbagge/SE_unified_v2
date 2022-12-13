#include <math.h>
#include "SE_Stresslet_direct.h"
#include "hasimoto_decomp.h"
#include "../incomplete_bessel_K.h"

#ifdef EWALD_FULL
void SE1P_Stresslet_direct_full(double* restrict u, const size_t nout,
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
            {
#ifndef EXTERNAL
                if (n == m_idx && j0 == 0)
                    continue;                                       // skip self
#endif
                // EXTERNAL: Assuming that r != 0 in home box

                rsh[0] = r[0] + j0*opt->box[0];
                rsh[1] = r[1];
                rsh[2] = r[2];

                double rdotq = dot(rsh, qn);
                double rdotn = dot(rsh, non);
                double ri2 = 1.0/dot(rsh,rsh);
                double ri = sqrt(ri2);
                double ri5 = ri2*ri2*ri;
                double C = -6*rdotq*rdotn*ri5;

                um[0] += C*rsh[0];
                um[1] += C*rsh[1];
                um[2] += C*rsh[2];
            }
        }
        u[m       ] = um[0];
        u[m+nout  ] = um[1];
        u[m+2*nout] = um[2];
    }
}
#endif // EWALD_FULL

#ifdef EWALD_REAL
void SE1P_Stresslet_direct_real(double* restrict u, const size_t nout,
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
            {
#ifndef EXTERNAL
                if (n == m_idx && j0 == 0)
                    continue;                                       // skip self
#endif
                // EXTERNAL: Assuming that r != 0 in home box

                rsh[0] = r[0] + j0*opt->box[0];
                rsh[1] = r[1];
                rsh[2] = r[2];
#ifdef CUTOFF
                if (dot(rsh,rsh) >= rc2) continue;            // skip outside rc
#endif
                real_part_stresslet_contribution(um, rsh, qn, non, opt->xi);
            }
        }
        u[m       ] = um[0];
        u[m+nout  ] = um[1];
        u[m+2*nout] = um[2];
    }
}
#endif // EWALD_REAL

#ifdef EWALD_FOURIER
void SE1P_Stresslet_direct_fourier(double* restrict u, const size_t nout,
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
    const double PI2 = PI*PI;
    const double kfac = 2*PI/opt->box[0];
    const double invvol = 1.0/(4*PI2*opt->box[0]);
    const double xi = opt->xi;
    const double xi2 = xi*xi;
    const double xi4 = xi2*xi2;
    const double oxi2 = 1.0/xi2;

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

        double k, kk, kx, y, yy, z, zz, rho2, ckx, skx;
        double U, V, Euv, Km1, K0, K1, K2; // incomplete modified Bessel functions
        double Q000, Q111, Q222, Q001, Q002, Q011, Q012, Q022, Q112, Q122;
        double q0, q1, q2, n0, n1, n2;
        double Qf_re[3], Qf_im[3];
        for (ptrdiff_t j0=-kmax; j0<=kmax; j0++)             // for k-space cube
        {
            if (j0 == 0) continue;                                // exclude k=0

            k = kfac*j0;
            kk = k*k;
            U = 0.25*kk*oxi2;

            for (size_t n=0; n<N; n++)                      // for all particles
            {
                kx = k*(xm[0] - x[n    ]);
                y  =    xm[1] - x[n+N  ];
                z  =    xm[2] - x[n+2*N];
                yy = y*y;
                zz = z*z;
                rho2 = yy + zz;
                V = rho2*xi2;
                Euv = exp(-U-V);

                /* Compute incomplete modified Bessel functions */
                K0 = incomplete_bessel_K(0, U, V);
                K1 = incomplete_bessel_K(1, U, V);
                K2 = incomplete_bessel_K(2, U, V);
                Km1 = (Euv + V*K1)/U;

                /* Fill tensor entries (symmetric tensor) */
                Q000 = 2*PI2*k * (kk*oxi2*Km1 + 12*K0 - 12*xi2*rho2*K1); // IMAG
                Q111 = 4*PI2*y * (-3*kk*K0 - 12*xi2*K1 + 4*xi4*(rho2+2*zz)*K2); // REAL
                Q222 = 4*PI2*z * (-3*kk*K0 - 12*xi2*K1 + 4*xi4*(rho2+2*yy)*K2); // REAL
                Q001 = 4*PI2*y * (kk*K0 - 8*xi2*K1 + 4*xi4*rho2*K2); // REAL
                Q002 = 4*PI2*z * (kk*K0 - 8*xi2*K1 + 4*xi4*rho2*K2); // REAL
                Q011 = 2*PI2*k * (kk*oxi2*Km1 + 4*xi2*(yy-zz)*K1); // IMAG
                Q012 = 16*PI2*xi2*k*y*z*K1; // IMAG
                Q022 = 2*PI2*k * (kk*oxi2*Km1 + 4*xi2*(zz-yy)*K1); // IMAG
                Q112 = 4*PI2*z * (-kk*K0 - 4*xi2*K1 + 4*xi4*(zz-yy)*K2); // REAL
                Q122 = 4*PI2*y * (-kk*K0 - 4*xi2*K1 + 4*xi4*(yy-zz)*K2); // REAL

                q0 = q[n    ];
                q1 = q[n+N  ];
                q2 = q[n+2*N];
                n0 = no[n    ];
                n1 = no[n+N  ];
                n2 = no[n+2*N];

                /* Real part of Q*f */
                Qf_re[0] =              Q001*q0*n1 + Q002*q0*n2
                         + Q001*q1*n0
                         + Q002*q2*n0                          ;
                Qf_re[1] = Q001*q0*n0
                                      + Q111*q1*n1 + Q112*q1*n2
                                      + Q112*q2*n1 + Q122*q2*n2;
                Qf_re[2] = Q002*q0*n0
                                      + Q112*q1*n1 + Q122*q1*n2
                                      + Q122*q2*n1 + Q222*q2*n2;

                /* Imag part of Q*f */
                Qf_im[0] = Q000*q0*n0
                                      + Q011*q1*n1 + Q012*q1*n2
                                      + Q012*q2*n1 + Q022*q2*n2;
                Qf_im[1] =              Q011*q0*n1 + Q012*q0*n2
                         + Q011*q1*n0
                         + Q012*q2*n0                          ;
                Qf_im[2] =              Q012*q0*n1 + Q022*q0*n2
                         + Q012*q1*n0
                         + Q022*q2*n0                          ;

                /*
                   Sum contributions --
                   imag part enters with negative sign (from i*i)
                */
                ckx = cos(kx);
                skx = sin(kx);
                um[0] += ckx*Qf_re[0] - skx*Qf_im[0];
                um[1] += ckx*Qf_re[1] - skx*Qf_im[1];
                um[2] += ckx*Qf_re[2] - skx*Qf_im[2];
            }
        }
        u[m       ] = invvol*um[0];
        u[m+nout  ] = invvol*um[1];
        u[m+2*nout] = invvol*um[2];
    }
}
#endif // EWALD_FOURIER

#ifdef EWALD_FOURIER_K0
void SE1P_Stresslet_direct_fourier_k0(double* restrict u, const size_t nout,
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
    const double invvol = 1.0/(4*PI*PI*opt->box[0]);
    const double xi = opt->xi;
    const double xi2 = xi*xi;
    const double xi4_4 = 4*xi2*xi2;
    const double PI2_4 = 4*PI*PI;

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

        double y, yy, z, zz, rho2;
        double V, K1, K2; // incomplete modified Bessel functions
        double Q111, Q222, Q001, Q002, Q112, Q122;
        double q0, q1, q2, n0, n1, n2;
        double Qf[3];
        for (size_t n=0; n<N; n++)              // for all particles
        {
            y  = xm[1] - x[n+N  ];
            z  = xm[2] - x[n+2*N];
            yy = y*y;
            zz = z*z;
            rho2 = yy + zz;
            V = rho2*xi2;

            /* Compute incomplete modified Bessel functions */
            K1 = incomplete_bessel_K(1, 0, V);
            K2 = incomplete_bessel_K(2, 0, V);

            /* Fill tensor entries (symmetric tensor) */
            Q111 = PI2_4*y * (-12*xi2*K1 + xi4_4*(rho2+2*zz)*K2);
            Q222 = PI2_4*z * (-12*xi2*K1 + xi4_4*(rho2+2*yy)*K2);
            Q001 = PI2_4*y * (-8*xi2*K1 + xi4_4*rho2*K2);
            Q002 = PI2_4*z * (-8*xi2*K1 + xi4_4*rho2*K2);
            Q112 = PI2_4*z * (-4*xi2*K1 + xi4_4*(zz-yy)*K2);
            Q122 = PI2_4*y * (-4*xi2*K1 + xi4_4*(yy-zz)*K2);

            q0 = q[n    ];
            q1 = q[n+N  ];
            q2 = q[n+2*N];
            n0 = no[n    ];
            n1 = no[n+N  ];
            n2 = no[n+2*N];

            Qf[0] =              Q001*q0*n1 + Q002*q0*n2
                  + Q001*q1*n0
                  + Q002*q2*n0                          ;
            Qf[1] = Q001*q0*n0
                               + Q111*q1*n1 + Q112*q1*n2
                               + Q112*q2*n1 + Q122*q2*n2;
            Qf[2] = Q002*q0*n0
                               + Q112*q1*n1 + Q122*q1*n2
                               + Q122*q2*n1 + Q222*q2*n2;

            um[0] += Qf[0];
            um[1] += Qf[1];
            um[2] += Qf[2];
        }
        u[m       ] = invvol*um[0];
        u[m+nout  ] = invvol*um[1];
        u[m+2*nout] = invvol*um[2];
    }
}
#endif // EWALD_FOURIER_K0
