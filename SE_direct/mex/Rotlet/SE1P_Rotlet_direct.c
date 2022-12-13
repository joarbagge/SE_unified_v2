#include <math.h>
#include "SE_Rotlet_direct.h"
#include "rotlet_decomp.h"
#include "../incomplete_bessel_K.h"

#ifdef EWALD_FULL
void SE1P_Rotlet_direct_full(double* restrict u, const size_t nout,
                             const double* restrict x, const double* restrict t,
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
            {
#ifndef EXTERNAL
                if (n == m_idx && j0 == 0)
                    continue;                                       // skip self
#endif
                // EXTERNAL: Assuming that r != 0 in home box

                rsh[0] = r[0] + j0*opt->box[0];
                rsh[1] = r[1];
                rsh[2] = r[2];

                double ri2 = 1.0/dot(rsh,rsh);
                double ri = sqrt(ri2);
                double ri3 = ri2*ri;

                // um += ri3 * cross(tn, rsh)
                um[0] += ri3 * (tn[1]*rsh[2] - tn[2]*rsh[1]);
                um[1] += ri3 * (tn[2]*rsh[0] - tn[0]*rsh[2]);
                um[2] += ri3 * (tn[0]*rsh[1] - tn[1]*rsh[0]);
            }
        }
        u[m       ] = um[0];
        u[m+nout  ] = um[1];
        u[m+2*nout] = um[2];
    }
}
#endif // EWALD_FULL

#ifdef EWALD_REAL
void SE1P_Rotlet_direct_real(double* restrict u, const size_t nout,
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
                real_part_rotlet_contribution(um, rsh, tn, opt->xi);
            }
        }
        u[m       ] = um[0];
        u[m+nout  ] = um[1];
        u[m+2*nout] = um[2];
    }
}
#endif // EWALD_REAL

#ifdef EWALD_FOURIER
void SE1P_Rotlet_direct_fourier(double* restrict u, const size_t nout,
                                const double* restrict x, const double* restrict t,
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

        double k, kx, y, z, tmp, ckx, skx;
        double U, V, K0, K1; // incomplete modified Bessel functions
        double Q01, Q02, Q12;
        double t0, t1, t2;
        double Qt_re[3], Qt_im[3];
        for (ptrdiff_t j0=-kmax; j0<=kmax; j0++)             // for k-space cube
        {
            if (j0 == 0) continue;                                // exclude k=0

            k = kfac*j0;
            U = k*k/(4*xi2);

            for (size_t n=0; n<N; n++)                      // for all particles
            {
                kx = k*(xm[0] - x[n    ]);
                y  =    xm[1] - x[n+N  ];
                z  =    xm[2] - x[n+2*N];
                V = xi2*(y*y + z*z);

                /* Compute incomplete modified Bessel functions */
                K0 = incomplete_bessel_K(0, U, V);
                K1 = incomplete_bessel_K(1, U, V);

                /* Fill Q matrix entries (antisymmetric matrix) */
                tmp = 8*PI2*xi2*K1;
                Q01 = z*tmp; // REAL
                Q02 = -y*tmp; // REAL
                Q12 = -4*PI2*k*K0; // IMAG

                t0 = t[n    ];
                t1 = t[n+N  ];
                t2 = t[n+2*N];

                /* Real part of Q*t */
                Qt_re[0] = Q01*t1 + Q02*t2;
                Qt_re[1] = -Q01*t0;
                Qt_re[2] = -Q02*t0;

                /* Imag part of Q*t */
                Qt_im[0] = 0;
                Qt_im[1] = Q12*t2;
                Qt_im[2] = -Q12*t1;

                /*
                   Sum contributions --
                   imag part enters with negative sign (from i*i)
                */
                ckx = cos(kx);
                skx = sin(kx);
                um[0] += ckx*Qt_re[0] - skx*Qt_im[0];
                um[1] += ckx*Qt_re[1] - skx*Qt_im[1];
                um[2] += ckx*Qt_re[2] - skx*Qt_im[2];
            }
        }
        u[m       ] = invvol*um[0];
        u[m+nout  ] = invvol*um[1];
        u[m+2*nout] = invvol*um[2];
    }
}
#endif // EWALD_FOURIER

#ifdef EWALD_FOURIER_K0
void SE1P_Rotlet_direct_fourier_k0(double* restrict u, const size_t nout,
                                   const double* restrict x, const double* restrict t,
                                   const size_t N, const ewald_opts* restrict opt)
{
    /* Zero-mode Fourier-part rotlet */
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif
    const double xi = opt->xi;
    const double xi2 = xi*xi;
    const double prefac = 2*xi2/opt->box[0];

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

        double y, z, t0, t1, t2;
        double V, K1; // incomplete modified Bessel function
        for (size_t n=0; n<N; n++)                          // for all particles
        {
            y = xm[1] - x[n+N  ];
            z = xm[2] - x[n+2*N];
            V = xi2*(y*y + z*z);
            K1 = incomplete_bessel_K(1, 0, V);

            t0 = t[n    ];
            t1 = t[n+N  ];
            t2 = t[n+2*N];

            um[0] += K1 * (z*t1 - y*t2);
            um[1] += -K1 * z*t0;
            um[2] += K1 * y*t0;
        }
        u[m       ] = prefac*um[0];
        u[m+nout  ] = prefac*um[1];
        u[m+2*nout] = prefac*um[2];
    }
}
#endif // EWALD_FOURIER_K0
