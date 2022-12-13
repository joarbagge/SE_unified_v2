#include <math.h>
#include "SE_Laplace_direct.h"
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include "../incomplete_bessel_K.h"

inline double dot(double * a, double * b)
{
    return ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

#ifdef EWALD_FULL
void SE1P_Laplace_direct_full(double* restrict u, const size_t nout,
                              const double* restrict x, const double* restrict f,
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

        double rvec[3], r[3], ri, fn;
        for (size_t n=0; n<N; n++)                          // for all particles
        {
            fn = f[n];
            rvec[0] = xm[0] - x[n    ];
            rvec[1] = xm[1] - x[n+N  ];
            rvec[2] = xm[2] - x[n+2*N];

            for (ptrdiff_t j0 = -nbox; j0<=nbox; j0++)            // image boxes
            {
#ifndef EXTERNAL
                if (m_idx == n && j0 == 0)
                    continue;                                       // skip self
#endif
                // EXTERNAL: Assuming that r != 0 in home box

                r[0] = rvec[0] + j0*opt->box[0];
                r[1] = rvec[1];
                r[2] = rvec[2];

                ri = 1.0/sqrt(dot(r,r));

#ifdef GRADIENT
                double ri3 = -ri*ri*ri;
                um[0] += ri3*fn*r[0];
                um[1] += ri3*fn*r[1];
                um[2] += ri3*fn*r[2];
#else
                um += ri*fn;
#endif
            }
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
#endif // EWALD_FULL

#ifdef EWALD_REAL
void SE1P_Laplace_direct_real(double* restrict u, size_t nout,
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
            {
#ifndef EXTERNAL
                if (m_idx == n && j0 == 0)
                    continue;                                       // skip self
#endif
                // EXTERNAL: Assuming that r != 0 in home box

                rx = rvec[0] + j0*opt->box[0];
                ry = rvec[1];
                rz = rvec[2];
                r2 = rx*rx + ry*ry + rz*rz;
#ifdef CUTOFF
                if (r2 >= rc2) continue;                      // skip outside rc
#endif
                r = sqrt(r2);

#ifdef GRADIENT
                double c = -fn*(a*exp(-xi2*r2) + erfc(opt->xi*r)/r)/r2;
                p[0] += c*rx;
                p[1] += c*ry;
                p[2] += c*rz;
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
#endif // EWALD_REAL

#ifdef EWALD_FOURIER
/*
Equivalent MATLAB code (for potential):

function u = SE1P_Laplace_direct_fourier(x, f, opt)
N = numel(f);
xi = opt.xi;
xi2 = xi*xi;
TwoPiOverL = 2*pi/opt.box(1);
u = zeros(N,1);
for target=opt.eval_idx
    xt = x(target,:);
    u_t = 0;
    for source=1:N
        X = xt(1) - x(source,1);
        rho2 = norm(xt(2:3) - x(source,2:3))^2;
        b = rho2*xi2;
        fn = f(source);
        for p=1:opt.layers
            k = TwoPiOverL*p;
            a = k*k/(4*xi2);
            K0 = computeK0(a, b);
            kX = -k*X;
            u_t = u_t + 2*fn*cos(kX)*K0;
        end
    end
    u(target) = u_t/opt.box(1);
end

function v = computeK0(a, b)
func = @(t) exp(-a./t-b.*t)./t;
v = integral(func, 0, 1, 'reltol', 1e-15);
*/
void SE1P_Laplace_direct_fourier(double* restrict u, size_t nout,
                                 const double* restrict x, const double* restrict f, size_t N,
                                 const ewald_opts* restrict opt)
{
    double xm[3];
    const double xi = opt->xi;
    const double xi2 = xi*xi;
    const double TwoPiOverL = 2*PI/opt->box[0];
    const ptrdiff_t kmax = opt->layers;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif

#ifdef _OPENMP
#pragma omp parallel for private(xm)
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

        for (ptrdiff_t j0=-kmax; j0<=kmax; j0++)              // for wavenumbers
        {
            if (j0 == 0) continue;                             // skip zero mode
            double k = TwoPiOverL*j0;
            double a = k*k/(4*xi2);
            for (size_t n=0; n<N; n++)                      // for all particles
            {
                double r = xm[0] - x[n    ];
                double y = xm[1] - x[n+N  ];
                double z = xm[2] - x[n+2*N];
                double rho2 = y*y + z*z;
                double b = rho2*xi2;
                double kr = -k*r;
#ifdef GRADIENT
                double K0 = incomplete_bessel_K(0, a, b);
                double K1 = incomplete_bessel_K(1, a, b);
                p[0] += f[n]*k*sin(kr)*K0;
                p[1] += -2*f[n]*xi2*cos(kr)*y*K1;
                p[2] += -2*f[n]*xi2*cos(kr)*z*K1;
#else
                double K0 = incomplete_bessel_K(0, a, b);
                p += f[n]*cos(kr)*K0;
#endif
            }
        }
#ifdef GRADIENT
        u[m       ] = p[0]/opt->box[0];
        u[m+nout  ] = p[1]/opt->box[0];
        u[m+2*nout] = p[2]/opt->box[0];
#else
        u[m] = p/opt->box[0];
#endif
    }
}
#endif // EWALD_FOURIER

#ifdef EWALD_FOURIER_K0
/*
Equivalent MATLAB code (for potential):

function u = SE1P_Laplace_direct_fourier_k0(x, f, opt)
N = numel(f);
xi = opt.xi;
egamma = 0.57721566490153286061;
u = zeros(numel(opt.eval_idx), 1);
for j=1:numel(opt.eval_idx);
    target = opt.eval_idx(j);
    xt = x(target,2:3);
    u_t = 0;
    for source=1:N
        rho2 = norm(xt-x(source,2:3))^2;
        v = rho2*xi*xi;
        if v > eps
            u_t = u_t - f(source)*(expint(v) + log(v) + egamma);
        end
    end
    u(j) = u_t/opt.box(1);
end
*/
void SE1P_Laplace_direct_fourier_k0(double* restrict u, size_t nout,
                                    const double* restrict x, const double* restrict f, size_t N,
                                    const ewald_opts* restrict opt)
{
    double xm[3];
    const double xi = opt->xi;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif
#ifndef GRADIENT
    const double egamma = 0.57721566490153286061;
#endif
    gsl_set_error_handler_off();
    int gsl_error_code = GSL_SUCCESS;

#ifdef _OPENMP
#pragma omp parallel for private(xm)
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

        int status;
        gsl_sf_result result;
        for (size_t n=0; n<N; n++)                          // for all particles
        {
            double y = xm[1] - x[n+N  ];
            double z = xm[2] - x[n+2*N];
            double rho2 = y*y + z*z;
            double v = rho2*xi*xi;
#ifdef GRADIENT
            if (v > __DBL_EPSILON__)
            {
                p[1] += -f[n]*2*y/rho2*(1-exp(-v));
                p[2] += -f[n]*2*z/rho2*(1-exp(-v));
            }
#else
            /* If v > 34 then E1(v) < 5e-17, so save some time by not computing it then. */
            if (v > 34)
            {
                p += -f[n]*(log(v)+egamma);
            }
            else if (v > __DBL_EPSILON__)
            {
                status = gsl_sf_expint_E1_e(v, &result);
                if (status == GSL_EUNDRFLW) // underflow
                {
                    printf("SE:DirectGSL Warning: GSL underflow occurred.\n");
                }
                else if (status != GSL_SUCCESS)
                {
#ifdef _OPENMP
#pragma omp critical
#endif
                    if (gsl_error_code == GSL_SUCCESS) {
                        gsl_error_code = status;
                    }
                    break;
                }
                p += -f[n]*(result.val+log(v)+egamma);
            }
#endif
        }
#ifdef GRADIENT
        u[m       ] = 0;
        u[m+nout  ] = p[1]/opt->box[0];
        u[m+2*nout] = p[2]/opt->box[0];
#else
        u[m] = p/opt->box[0];
#endif
    }

    if (gsl_error_code != GSL_SUCCESS)
    {
        char errMsg [64];
        sprintf(errMsg, "GSL error code: %d", gsl_error_code);
        mexErrMsgIdAndTxt("SE:DirectGSL", errMsg);
    }
}
#endif // EWALD_FOURIER_K0

#ifdef EWALD_SELF
void SE1P_Laplace_direct_self(double* restrict u, size_t nout,
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
#endif // EWALD_SELF
