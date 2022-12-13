#include <math.h>
#include "SE_Laplace_direct.h"

void SE3P_Laplace_direct_real(double* restrict u, size_t nout,
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
                    for (ptrdiff_t j2 = -nbox; j2<=nbox; j2++)
                    {
#ifndef EXTERNAL
                        if (m_idx == n && j0 == 0 && j1 == 0 && j2 == 0)
                            continue;                               // skip self
#endif
                        // EXTERNAL: Assuming that r != 0 in home box

                        rx = rvec[0] + j0*opt->box[0];
                        ry = rvec[1] + j1*opt->box[1];
                        rz = rvec[2] + j2*opt->box[2];
                        r2 = rx*rx + ry*ry + rz*rz;
#ifdef CUTOFF
                        if (r2 >= rc2) continue;              // skip outside rc
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

/*
Equivalent MATLAB code (for potential, not external):

function u = SE3P_Laplace_direct_fourier(x, f, opt)
  N = numel(f);
  xi = opt.xi;
  xi2 = xi*xi;
  TwoPiOverL = 2*pi./opt.box;
  c = 4*pi/(opt.box(1)*opt.box(2)*opt.box(3));
  a = 1/(4*xi2);
  u = zeros(N,1);
  for target=opt.eval_idx
    xt = x(target,:);
    u_t = 0;
    for j1=-opt.layers:opt.layers
      for j2=-opt.layers:opt.layers
        for j3=-opt.layers:opt.layers
          if j1 == 0 && j2 == 0 && j3 == 0
            continue
          end
          k = TwoPiOverL .* [j1 j2 j3];
          k2 = dot(k,k);
          z = 0;
          for source=1:N
            X = xt - x(source,:);
            fn = f(source);
            z = z + fn*cos(-dot(k,X)); % real part of exp(-i*dot(k,X))
                                       % (imaginary part cancels)
          end
          u_t = u_t + c * z * exp(-a*k2)/k2;
        end
      end
    end
    u(target) = u_t;
  end
end
*/
void SE3P_Laplace_direct_fourier(double* restrict u, size_t nout,
                                 const double* restrict x, const double* restrict f, size_t N,
                                 const ewald_opts* restrict opt)
{
    double k[3];
    double k2, z;
    const double c = 4*PI/(opt->box[0]*opt->box[1]*opt->box[2]);
    const double fac[3] = {2*PI/opt->box[0], 2*PI/opt->box[1], 2*PI/opt->box[2]};
    const double a = 1.0/(4*opt->xi*opt->xi);
    const ptrdiff_t kmax = opt->layers;
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif

#ifdef _OPENMP
#pragma omp parallel for private(k, k2, z)
#endif
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
#ifdef GRADIENT
        double p[] = {0, 0, 0};
#else
        double p = 0;
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

        for (ptrdiff_t j0 = -kmax; j0<=kmax; j0++)           // for k-space cube
            for (ptrdiff_t j1 = -kmax; j1<=kmax; j1++)
                for (ptrdiff_t j2 = -kmax; j2<=kmax; j2++)
                {
                    if (j0 == 0 && j1 == 0 && j2 == 0)            // exclude k=0
                        continue;

                    k[0] = fac[0]*j0;
                    k[1] = fac[1]*j1;
                    k[2] = fac[2]*j2;
                    k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];

                    z = 0;
                    double k_dot_r;
                    for (size_t n=0; n<N; n++)              // for all particles
                    {
                        k_dot_r = k[0]*(xm[0]-x[n    ])+
                                  k[1]*(xm[1]-x[n+N  ])+
                                  k[2]*(xm[2]-x[n+2*N]);
                        /* Only the real part remains here, since the imaginary part cancels */
#ifdef GRADIENT
                        z += -f[n]*sin(k_dot_r);
#else
                        z += f[n]*cos(k_dot_r);
#endif
                    }
#ifdef GRADIENT
                    double tmp = z*exp(-a*k2)/k2;
                    p[0] += k[0]*tmp;
                    p[1] += k[1]*tmp;
                    p[2] += k[2]*tmp;
#else
                    p += z*exp(-a*k2)/k2;
#endif

                }
#ifdef GRADIENT
        u[m       ] = c*p[0];
        u[m+nout  ] = c*p[1];
        u[m+2*nout] = c*p[2];
#else
        u[m] = c*p;
#endif
    }
}

void SE3P_Laplace_direct_self(double* restrict u, size_t nout,
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
