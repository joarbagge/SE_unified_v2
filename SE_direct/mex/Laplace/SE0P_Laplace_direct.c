#include <math.h>
#include "SE_Laplace_direct.h"

inline double dot(double * a, double * b)
{
    return ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

#ifdef EWALD_FULL
/*
 * Equivalent MATLAB code for potential (not external):
M = numel(opt.eval_idx);
N = size(x, 1);
u = zeros(N, 1);
for i=1:M
  target = opt.eval_idx(i);
  source = [1:target-1 target+1:N];
  rvec = bsxfun(@minus, x(target,:), x(source,:));
  dist = sqrt(sum(rvec.^2, 2));
  u(i) = sum(f(source)./dist);
end

 * Equivalent MATLAB code for gradient (not external):
M = numel(opt.eval_idx);
N = size(x, 1);
du = zeros(N, 3);
for i=1:M
  target = opt.eval_idx(i);
  source = [1:target-1 target+1:N];
  rvec = bsxfun(@minus, x(target,:), x(source,:));
  dist = sqrt(sum(rvec.^2, 2));
  u0 = f(source)./dist.^3;
  du(i,:) = -sum(bsxfun(@times, u0, rvec), 1);
end
 */
void SE0P_Laplace_direct_full(double* restrict u, const size_t nout,
                              const double* restrict x, const double* restrict f,
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
    for (size_t mi=0; mi<nout; mi++)                // for all evaluation points
    {
#ifdef GRADIENT
        double p[] = {0, 0, 0};
#else
        double p = 0;
#endif

#ifndef EXTERNAL
        size_t m;
        if (opt->N_eval == -1)
            m = mi;               // default is to evaluate at all source points
        else
            m = idx[mi];                   // indirect indexing OK in outer loop
#endif

        for (size_t n=0; n<N; n++)                          // for all particles
        {
#ifndef EXTERNAL
            if (m == n) continue;                                   // skip self
#endif
            // EXTERNAL: Assuming that r != 0

            double fn = f[n];
#ifdef EXTERNAL
            double r[] = {xt[mi]-x[n], xt[mi+nout]-x[n+N], xt[mi+2*nout]-x[n+2*N]};
#else
            double r[] = {x[m]-x[n], x[m+N]-x[n+N], x[m+2*N]-x[n+2*N]};
#endif
            double ri = 1.0/sqrt(dot(r,r));

#ifdef GRADIENT
            double ri3 = -ri*ri*ri;
            p[0] += ri3*fn*r[0];
            p[1] += ri3*fn*r[1];
            p[2] += ri3*fn*r[2];
#else
            p += ri*fn;
#endif
        }
#ifdef GRADIENT
        u[mi       ] = p[0];
        u[mi+nout  ] = p[1];
        u[mi+2*nout] = p[2];
#else
        u[mi] = p;
#endif
    }
}
#endif // EWALD_FULL

#ifdef EWALD_REAL
void SE0P_Laplace_direct_real(double* restrict u, const size_t nout,
                              const double* restrict x, const double* restrict f,
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
#ifdef GRADIENT
    const double xi2 = opt->xi * opt->xi;
    const double a = 2*opt->xi/sqrt(PI);
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
#ifdef GRADIENT
        double um[] = {0, 0, 0};
#else
        double um = 0;
#endif

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

        double r[3], r2, r1, fn;
        for (size_t n=0; n<N; n++)                          // for all particles
        {
#ifndef EXTERNAL
            if (n == m_idx) continue;                               // skip self
#endif
            // EXTERNAL: Assuming that r != 0

            r[0] = xm[0] - x[n    ];
            r[1] = xm[1] - x[n+N  ];
            r[2] = xm[2] - x[n+2*N];
            r2 = dot(r,r);

#ifdef CUTOFF
            if (r2 >= rc2) continue;                          // skip outside rc
#endif
            r1 = sqrt(r2);
            fn = f[n];

#ifdef GRADIENT
            double c = -fn*(a*exp(-xi2*r2) + erfc(opt->xi*r1)/r1)/r2;
            um[0] += c*r[0];
            um[1] += c*r[1];
            um[2] += c*r[2];
#else
            um += fn*erfc(opt->xi*r1)/r1;
#endif
        }
#ifdef GRADIENT
        u[m       ] = um[0];
        u[m+nout  ] = um[1];
        u[m+2*nout] = um[2];
#else
        u[m] = um;
#endif
    }
}
#endif // EWALD_REAL
