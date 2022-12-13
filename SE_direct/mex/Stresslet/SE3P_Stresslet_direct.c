#include <math.h>
#include "SE_Stresslet_direct.h"
#include "hasimoto_decomp.h"

void SE3P_Stresslet_direct_real(double* restrict u, const size_t nout,
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
                    for (ptrdiff_t j2=-nbox; j2<=nbox; j2++)
                    {
#ifndef EXTERNAL
                        if (n == m_idx && j0 == 0 && j1 == 0 && j2 == 0)
                            continue;                               // skip self
#endif
                        // EXTERNAL: Assuming that r != 0 in home box

                        rsh[0] = r[0] + j0*opt->box[0];
                        rsh[1] = r[1] + j1*opt->box[1];
                        rsh[2] = r[2] + j2*opt->box[2];
#ifdef CUTOFF
                        if (dot(rsh,rsh) >= rc2) continue;    // skip outside rc
#endif
                        real_part_stresslet_contribution(um, rsh, qn, non, opt->xi);
                    }
        }
        u[m       ] = um[0];
        u[m+nout  ] = um[1];
        u[m+2*nout] = um[2];
    }
}

void SE3P_Stresslet_direct_fourier(double* restrict u, const size_t nout,
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
    const double kfac[] = {2*PI/opt->box[0], 2*PI/opt->box[1], 2*PI/opt->box[2]};
    const double invvol = 1.0/(opt->box[0]*opt->box[1]*opt->box[2]);

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

        double k[3], Z[3][3], k_dot_r, sinkr;
        for (ptrdiff_t j0=-kmax; j0<=kmax; j0++)             // for k-space cube
            for (ptrdiff_t j1=-kmax; j1<=kmax; j1++)
                for (ptrdiff_t j2=-kmax; j2<=kmax; j2++)
                {
                    if (j0 == 0 && j1 == 0 && j2 == 0)
                        continue;                                 // exclude k=0

                    k[0] = kfac[0]*j0;
                    k[1] = kfac[1]*j1;
                    k[2] = kfac[2]*j2;

                    for (int i=0; i<3; i++) for (int j=0; j<3; j++) Z[i][j] = 0;
                    for (size_t n=0; n<N; n++)              // for all particles
                    {
                        k_dot_r = k[0]*(xm[0] - x[n    ])+
                                  k[1]*(xm[1] - x[n+N  ])+
                                  k[2]*(xm[2] - x[n+2*N]);
                        /* Only the real part remains here, the imaginary part cancels */
                        /* Real part will be multiplied by -sin */
                        sinkr = -sin(k_dot_r);
                        for (int i=0; i<3; i++)
                            for (int j=0; j<3; j++)
                            {
                                Z[i][j] += sinkr * q[n+i*N] * no[n+j*N];
                            }
                    }
                    fourier_part_stresslet_contribution(um, k, Z, opt->xi);
                }
        u[m       ] = invvol*um[0];
        u[m+nout  ] = invvol*um[1];
        u[m+2*nout] = invvol*um[2];
    }
}

void SE3P_Stresslet_direct_fourier_k0(double* restrict u, const size_t nout,
                                      const double* restrict x,
                                      const double* restrict q, const double* restrict no,
                                      const size_t N, const ewald_opts* restrict opt)
{
    /* Zero-mode Fourier-part stresslet TF(0)_jlm = -8 pi d_lm r_j */
    /* Computes (1/V) sum_n TF(0)_jlm q_l n_m */
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif
    const double volfac = -8*PI/(opt->box[0]*opt->box[1]*opt->box[2]);

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

        double r[3], q_dot_no;
        for (size_t n=0; n<N; n++)                          // for all particles
        {
            r[0] = xm[0] - x[n    ];
            r[1] = xm[1] - x[n+N  ];
            r[2] = xm[2] - x[n+2*N];
            q_dot_no = q[n    ]*no[n    ]+
                       q[n+N  ]*no[n+N  ]+
                       q[n+2*N]*no[n+2*N];
            um[0] += r[0] * q_dot_no;
            um[1] += r[1] * q_dot_no;
            um[2] += r[2] * q_dot_no;
        }
        u[m       ] = volfac*um[0];
        u[m+nout  ] = volfac*um[1];
        u[m+2*nout] = volfac*um[2];
    }
}
