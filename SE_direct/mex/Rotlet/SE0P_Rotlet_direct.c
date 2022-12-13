#include <math.h>
#include "SE_Rotlet_direct.h"
#include "rotlet_decomp.h"

#ifdef EWALD_FULL
void SE0P_Rotlet_direct_full(double* restrict u, const size_t nout,
                             const double* restrict x, const double* restrict t,
                             const size_t N, const ewald_opts* restrict opt)
{
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

        double r[3], tn[3];
        for (size_t n=0; n<N; n++)                          // for all particles
        {
            r[0] = xm[0] - x[n    ];
            r[1] = xm[1] - x[n+N  ];
            r[2] = xm[2] - x[n+2*N];
            tn[0] = t[n    ];
            tn[1] = t[n+N  ];
            tn[2] = t[n+2*N];

#ifndef EXTERNAL
            if (n == m_idx) continue;                               // skip self
#endif
            // EXTERNAL: Assuming that r != 0 in home box

            double ri2 = 1.0/dot(r,r);
            double ri = sqrt(ri2);
            double ri3 = ri2*ri;

            // um += ri3 * cross(tn, r)
            um[0] += ri3 * (tn[1]*r[2] - tn[2]*r[1]);
            um[1] += ri3 * (tn[2]*r[0] - tn[0]*r[2]);
            um[2] += ri3 * (tn[0]*r[1] - tn[1]*r[0]);
        }
        u[m       ] = um[0];
        u[m+nout  ] = um[1];
        u[m+2*nout] = um[2];
    }
}
#endif // EWALD_FULL

#ifdef EWALD_REAL
void SE0P_Rotlet_direct_real(double* restrict u, const size_t nout,
                             const double* restrict x, const double* restrict t,
                             const size_t N, const ewald_opts* restrict opt)
{
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

        double r[3], tn[3];
        for (size_t n=0; n<N; n++)                          // for all particles
        {
            r[0] = xm[0] - x[n    ];
            r[1] = xm[1] - x[n+N  ];
            r[2] = xm[2] - x[n+2*N];
            tn[0] = t[n    ];
            tn[1] = t[n+N  ];
            tn[2] = t[n+2*N];

#ifndef EXTERNAL
            if (n == m_idx) continue;                               // skip self
#endif
            // EXTERNAL: Assuming that r != 0 in home box

#ifdef CUTOFF
            if (dot(r,r) >= rc2) continue;                    // skip outside rc
#endif
            real_part_rotlet_contribution(um, r, tn, opt->xi);
        }
        u[m       ] = um[0];
        u[m+nout  ] = um[1];
        u[m+2*nout] = um[2];
    }
}
#endif // EWALD_REAL
