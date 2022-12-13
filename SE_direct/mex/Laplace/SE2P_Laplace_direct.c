#include <math.h>
#include "SE_Laplace_direct.h"

#define __EXP_ARG_MAX 600

void SE2P_Laplace_direct_real(double* restrict u, size_t nout,
                              const double* restrict x, const double* restrict f, size_t N,
                              const ewald_opts* restrict opt)
{
    const ptrdiff_t nbox = opt->layers;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif

    double xm[3], rvec[3];
    double fn;
    double r, r2, rx, ry, rz;
#ifdef CUTOFF
    const double rc2 = opt->rc * opt->rc;
#endif
#ifdef GRADIENT
    const double xi2 = opt->xi * opt->xi;
    const double a = 2*opt->xi/sqrt(PI);
#endif

#ifdef _OPENMP
#pragma omp parallel for private(xm, rvec, fn, r, r2, rx, ry, rz)
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
#ifdef GRADIENT
        double p[] = {0, 0, 0};
#else
        double p = 0;
#endif

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

        for (size_t n=0; n<N; n++)                          // for all particles
        {
            rvec[0] = xm[0] - x[n    ];
            rvec[1] = xm[1] - x[n+N  ];
            rvec[2] = xm[2] - x[n+2*N];
            fn = f[n];

            for (ptrdiff_t j0 = -nbox; j0<=nbox; j0++)            // image boxes
                for (ptrdiff_t j1 = -nbox; j1<=nbox; j1++)
                {
#ifndef EXTERNAL
                    if (m_idx == n && j0 == 0 && j1 == 0)
                        continue;                                   // skip self
#endif
                    // EXTERNAL: Assuming that r != 0 in home box

                    rx = rvec[0] + j0*opt->box[0];
                    ry = rvec[1] + j1*opt->box[1];
                    rz = rvec[2];
                    r2 = rx*rx + ry*ry + rz*rz;
#ifdef CUTOFF
                    if (r2 >= rc2) continue;                  // skip outside rc
#endif
                    r = sqrt(r2);

#ifdef GRADIENT
                    double tmp = erfc(opt->xi*r)/r + a*exp(-xi2*r2);
                    tmp = -fn*tmp/r2;
                    p[0] += rx*tmp;
                    p[1] += ry*tmp;
                    p[2] += rz*tmp;
#else
                    p += fn*erfc(opt->xi*r)/r;
#endif
                }
        }
#ifdef GRADIENT
        u[m       ] = p[0];
        u[m+nout  ] = p[1];
        u[m+2*nout] = p[2];
#else
        u[m] = p;
#endif
    }
}

static inline double theta_plus(double z, double k, double xi)
{
    /* idea for a more stable form [LK] */
    /* exp( k*z + log( erfc(k/(2.0*xi) + xi*z) ) ); */
    if (k*z < __EXP_ARG_MAX)
        return exp(k*z)*erfc(k/(2.0*xi) + xi*z);
    else
        return 0.0;
}

static inline double theta_minus(double z, double k, double xi)
{
    /* exp(-k*z + log( erfc(k/(2.0*xi) - xi*z) ) ); */
    if (-k*z < __EXP_ARG_MAX)
        return exp(-k*z)*erfc(k/(2.0*xi) - xi*z);
    else
        return 0.0;
}

static inline double theta_mod_plus(double z, double k, double xi)
{
    const double v = k/(2.0*xi) + xi*z;
    const double arg = k*z - v*v;
    if (arg < __EXP_ARG_MAX)
        return (2*xi/sqrt(PI))*exp(arg);
    else
        return 0.0;
}

static inline double theta_mod_minus(double z, double k, double xi)
{
    const double v = k/(2.0*xi) - xi*z;
    const double arg = -k*z - v*v;
    if (arg < __EXP_ARG_MAX)
        return (2*xi/sqrt(PI))*exp(arg);
    else
        return 0.0;
}

void SE2P_Laplace_direct_fourier(double* restrict u, size_t nout,
                                 const double* restrict x, const double* restrict f, size_t N,
                                 const ewald_opts* restrict opt)
{
    double k[2], xm[3];
    double kn, k_dot_r, z;
    double cm, cp;
    const double xi = opt->xi;
    const ptrdiff_t kmax = opt->layers;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif
    const double C = PI/(opt->box[0]*opt->box[1]);

#ifdef _OPENMP
#pragma omp parallel for private(k, xm, kn, k_dot_r, z, cm, cp)
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
#ifdef GRADIENT
        double p[] = {0, 0, 0};
        double dp, dm, tmp;
#else
        double p = 0;
#endif

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

        for (size_t n=0; n<N; n++)                          // for all particles
        {
            for (ptrdiff_t j0 = -kmax; j0<=kmax; j0++)       // for k-space cube
                for (ptrdiff_t j1 = -kmax; j1<=kmax; j1++)
                {
                    if (j0 == 0 && j1 == 0)                       // exclude k=0
                        continue;

                    k[0] = 2*PI*j0/opt->box[0];
                    k[1] = 2*PI*j1/opt->box[1];
                    kn = sqrt(k[0]*k[0] + k[1]*k[1]);
                    k_dot_r = k[0]*(xm[0]-x[n]) + k[1]*(xm[1]-x[n+N]);
                    z = xm[2]-x[n+2*N];
                    cp = theta_plus(z,kn,xi);
                    cm = theta_minus(z,kn,xi);
#ifdef GRADIENT
                    dp = theta_mod_plus(z,kn,xi);
                    dm = theta_mod_minus(z,kn,xi);
                    tmp = -f[n]*sin(k_dot_r)*(cp+cm)/kn;
                    p[0] += k[0]*tmp;
                    p[1] += k[1]*tmp;
                    p[2] += f[n]*cos(k_dot_r)*(kn*cp-dp-kn*cm+dm)/kn;
#else
                    p += f[n]*cos(k_dot_r)*(cm+cp)/kn;
#endif
                }
        }
#ifdef GRADIENT
        u[m       ] = C*p[0];
        u[m+nout  ] = C*p[1];
        u[m+2*nout] = C*p[2];
#else
        u[m] = C*p;
#endif
    }
}

/*
Equivalent MATLAB code:

function u = SE2P_Laplace_direct_fourier_k0(x, f, opt)
N = numel(f);
xi = opt.xi;
u = zeros(numel(opt.eval_idx), 1);
for target=opt.eval_idx
    zt = x(target,3);
    u_t = 0;
    for source=1:N
        r = zt - x(source,3);
        u_t = u_t + f(source) * (exp(-xi^2*r^2)/xi+sqrt(pi)*r*erf(xi*r));
    end
    u(target) = -2*sqrt(pi)*u_t/opt.box(1)/opt.box(2);
end
*/
void SE2P_Laplace_direct_fourier_k0(double* restrict u, size_t nout,
                                    const double* restrict x, const double* restrict f, size_t N,
                                    const ewald_opts* restrict opt)
{
    double z, zm;
    const double xi = opt->xi;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif
    const double C = -2*sqrt(PI)/(opt->box[0]*opt->box[1]);

#ifdef _OPENMP
#pragma omp parallel for private(z, zm)
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
#ifdef GRADIENT
        double p[] = {0, 0, 0};
#else
        double p = 0;
#endif

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

        for (size_t n=0; n<N; n++)                          // for all particles
        {
            z = zm - x[n+2*N];
#ifdef GRADIENT
            p[2] += f[n]*sqrt(PI)*erf(xi*z);
#else
            p += f[n]*(exp(-xi*xi*z*z)/xi +
                       sqrt(PI)*z*erf(xi*z));
#endif
        }
#ifdef GRADIENT
        //u[m       ] = C*p[0]; // is zero
        //u[m+nout  ] = C*p[1]; // is zero
        u[m+2*nout] = C*p[2];
#else
        u[m] = C*p;
#endif
    }
}

void SE2P_Laplace_direct_self(double* restrict u, size_t nout,
                              const double* restrict f, size_t N,
                              const ewald_opts* restrict opt)
{
    const size_t* restrict idx = opt->eval_idx;
    const double c = -2*opt->xi/sqrt(PI);
    size_t m_idx;
#ifdef _OPENMP
#pragma omp parallel for private(m_idx)
#endif
    for (size_t m=0; m<nout; m++)
    {
        if (opt->N_eval == -1)
            m_idx = m;            // default is to evaluate at all source points
        else
            m_idx = idx[m];
        u[m] = c*f[m_idx];
    }
}
