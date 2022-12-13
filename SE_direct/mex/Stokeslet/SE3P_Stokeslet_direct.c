#include <math.h>
#include "SE_Stokeslet_direct.h"

#ifdef HASIMOTO
#include "hasimoto_decomp.h"
#elif BEENAKKER
#include "beenakker_decomp.h"
#else
#error "Must provide -D<decomposition> to compiler"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

void SE3P_Stokeslet_direct_real(double* restrict u, size_t nout,
                                const double* restrict x, const double* restrict f, size_t N,
                                const ewald_opts* restrict opt)
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

    double r[3];
    double A[3][3];
    ptrdiff_t i1, i2, i3;
    size_t n;

#ifdef _OPENMP
#pragma omp parallel for private(r,A,i1,i2,i3,n)
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

        for (i1 = -nbox; i1<=nbox; i1++)                          // image boxes
            for (i2 = -nbox; i2<=nbox; i2++)
                for (i3 = -nbox; i3<=nbox; i3++)
                {
                    for (n=0; n<N; n++)                     // for all particles
                    {
#ifndef EXTERNAL
                        if (i1==0 && i2==0 && i3==0 && n==m_idx)    // skip self
                            continue;
#endif
                        // EXTERNAL: Assuming that r != 0 in home box

                        r[0] = xm[0]-x[n    ]+opt->box[0]*i1;
                        r[1] = xm[1]-x[n+  N]+opt->box[1]*i2;
                        r[2] = xm[2]-x[n+2*N]+opt->box[2]*i3;

#ifdef CUTOFF
                        if (r[0]*r[0] + r[1]*r[1] + r[2]*r[2] >= rc2)
                            continue;                         // skip outside rc
#endif

                        op_A(A, r, opt->xi);                         // u += A*f
                        um[0] +=
                            A[0][0]*f[n]+A[0][1]*f[n+N]+A[0][2]*f[n+2*N];
                        um[1] +=
                            A[1][0]*f[n]+A[1][1]*f[n+N]+A[1][2]*f[n+2*N];
                        um[2] +=
                            A[2][0]*f[n]+A[2][1]*f[n+N]+A[2][2]*f[n+2*N];
                    }
                }
        u[m       ] = um[0];
        u[m+nout  ] = um[1];
        u[m+2*nout] = um[2];
    }
}

void SE3P_Stokeslet_direct_fourier(double* restrict u, size_t nout,
                                   const double* restrict x, const double* restrict f, size_t N,
                                   const ewald_opts* restrict opt)
{
    double B[3][3];
    double z[3];
    double k[3];
    size_t n;
    ptrdiff_t i1, i2, i3;
    double q, k_dot_r;

    const double vol = opt->box[0]*opt->box[1]*opt->box[2];
    const ptrdiff_t kmax = opt->layers;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif
    const double kc0 = 2.0*PI/opt->box[0];
    const double kc1 = 2.0*PI/opt->box[1];
    const double kc2 = 2.0*PI/opt->box[2];

#ifdef _OPENMP
#pragma omp parallel for private(B,z,k,n,i1,i2,i3,q,k_dot_r)
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

        for (i1 = -kmax; i1<=kmax; i1++)                     // for k-space cube
            for (i2 = -kmax; i2<=kmax; i2++)
                for (i3 = -kmax; i3<=kmax; i3++)
                {
                    if (i3 != 0 || i2 != 0 || i1 != 0)            // exclude k=0
                    {
                        z[0] = 0; z[1] = 0; z[2] = 0;
                        k[0] = kc0*i1;
                        k[1] = kc1*i2;
                        k[2] = kc2*i3;

                        for (n=0; n<N; n++)                 // for all particles
                        {
                            k_dot_r =
                                k[0]*(xm[0]-x[n    ])+
                                k[1]*(xm[1]-x[n+N  ])+
                                k[2]*(xm[2]-x[n+2*N]);
                            q = cos(k_dot_r);
                            z[0] += q*f[n    ];
                            z[1] += q*f[n+N  ];
                            z[2] += q*f[n+2*N];
                        }
                        op_B(B, k, opt->xi);                   // multiplication
                        um[0] += B[0][0]*z[0]+B[0][1]*z[1]+B[0][2]*z[2];
                        um[1] += B[1][0]*z[0]+B[1][1]*z[1]+B[1][2]*z[2];
                        um[2] += B[2][0]*z[0]+B[2][1]*z[1]+B[2][2]*z[2];
                    }
                }
        u[m       ] = um[0]/vol;
        u[m+nout  ] = um[1]/vol;
        u[m+2*nout] = um[2]/vol;
    }
}

void SE3P_Stokeslet_direct_self(double* restrict u, size_t nout,
                                const double* restrict f, size_t N,
                                const ewald_opts* restrict opt)
{
    const size_t* restrict idx = opt->eval_idx;
    const double c = self_coeff(opt->xi);
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
