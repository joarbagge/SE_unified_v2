#include <math.h>
#include "SE_Stresslet_direct.h"
#include "hasimoto_decomp.h"

#ifdef EWALD_FULL
/* Equivalent MATLAB code (not external):
M = numel(opt.eval_idx);
N = size(x, 1);
u = zeros(M, 3);
for i=1:M
  target = opt.eval_idx(i);
  source = [1:target-1 target+1:N];
  rvec = bsxfun(@minus, x(target,:), x(source,:));
  ri5 = 1./sqrt(sum(rvec.^2, 2)).^5;
  rdotq = sum(rvec .* q(source,:), 2);
  rdotn = sum(rvec .* no(source,:), 2);
  u(i,:) = -6 * sum( bsxfun(@times, rvec, rdotn.*rdotq.*ri5), 1);
end
 */
void SE0P_Stresslet_direct_full(double* restrict u, const size_t nout,
                                const double* restrict x,
                                const double* restrict q, const double* restrict no,
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

        double r[3], qn[3], non[3];
        for (size_t n=0; n<N; n++)                          // for all particles
        {
#ifndef EXTERNAL
            if (n == m_idx) continue;                               // skip self
#endif
            // EXTERNAL: Assuming that r != 0 in home box

            r[0] = xm[0] - x[n    ];
            r[1] = xm[1] - x[n+N  ];
            r[2] = xm[2] - x[n+2*N];
            qn[0]  =  q[n    ];
            qn[1]  =  q[n+N  ];
            qn[2]  =  q[n+2*N];
            non[0] = no[n    ];
            non[1] = no[n+N  ];
            non[2] = no[n+2*N];

            double rdotq = dot(r, qn);
            double rdotn = dot(r, non);
            double ri2 = 1.0/dot(r, r);
            double ri = sqrt(ri2);
            double ri5 = ri2*ri2*ri;
            double C = -6*rdotq*rdotn*ri5;

            um[0] += C*r[0];
            um[1] += C*r[1];
            um[2] += C*r[2];
        }
        u[m       ] = um[0];
        u[m+nout  ] = um[1];
        u[m+2*nout] = um[2];
    }
}
#endif // EWALD_FULL

#ifdef EWALD_REAL
void SE0P_Stresslet_direct_real(double* restrict u, const size_t nout,
                                const double* restrict x,
                                const double* restrict q, const double* restrict no,
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

        double r[3], qn[3], non[3];
        for (size_t n=0; n<N; n++)                          // for all particles
        {
#ifndef EXTERNAL
            if (n == m_idx) continue;                               // skip self
#endif
            // EXTERNAL: Assuming that r != 0 in home box

            r[0] = xm[0] - x[n    ];
            r[1] = xm[1] - x[n+N  ];
            r[2] = xm[2] - x[n+2*N];

#ifdef CUTOFF
            if (dot(r,r) >= rc2) continue;                    // skip outside rc
#endif

            qn[0]  =  q[n    ];
            qn[1]  =  q[n+N  ];
            qn[2]  =  q[n+2*N];
            non[0] = no[n    ];
            non[1] = no[n+N  ];
            non[2] = no[n+2*N];

            real_part_stresslet_contribution(um, r, qn, non, opt->xi);
        }
        u[m       ] = um[0];
        u[m+nout  ] = um[1];
        u[m+2*nout] = um[2];
    }
}
#endif // EWALD_REAL
