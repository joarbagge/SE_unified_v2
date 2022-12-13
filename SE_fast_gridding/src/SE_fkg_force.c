#ifdef FORCE
#include "SE_fkg.h"

#ifdef THREE_PERIODIC
static
int kaiser_expansion_3p_force(const double x[3], const double q,
                              const SE_FGG_params* params,
                              double z2_0[P_MAX],
                              double z2_1[P_MAX],
                              double z2_2[P_MAX],
                              double zf_0[P_MAX],
                              double zf_1[P_MAX],
                              double zf_2[P_MAX])
{
    // unpack params
    const int p = params->P;
    const int p_half = params->P_half;
    const double h = params->h;
    const double oh = 1./h;
    const double w = params->P/2.;
    const double ow2 = 1./(w*w);
    const double beta = params->beta;
    double t0[3];
    int shift_t0;
    double scale_factor;
    KaiserlikeWindow window = params->window;
    KaiserlikeWindow window_deriv = params->window_deriv;
    SE_window_params window_opt;
    if (params->use_polynomial_window) {
      window_opt.polynomial_degree = params->polynomial_degree;
      shift_t0 = 0;
      scale_factor = oh/w;
    } else {
      window_opt.ow2 = ow2;
      window_opt.beta = beta;
      shift_t0 = 1;
      scale_factor = oh;
    }

    int idx;
    int idx_from[3];

    // compute index range and centering
    if(is_odd(p)) {
      for(int j=0; j<3; j++)
        {
          idx = (int) round(x[j]*oh);
          idx_from[j] = idx - p_half;
          t0[j] = x[j]*oh - idx + shift_t0*p_half;
        }
    }
    else {
      for(int j=0; j<3; j++)
        {
          idx = (int) floor(x[j]*oh);
          idx_from[j] = idx - (p_half-1);
          t0[j] = x[j]*oh - idx + shift_t0*(p_half-1);
        }
    }

    // compute window function and its derivative
    window(t0, p, &window_opt, z2_0, z2_1, z2_2);
    window_deriv(t0, p, &window_opt, zf_0, zf_1, zf_2);

    // extra scale factors applied to derivative
    for (int i=0; i<p; i++)
    {
        zf_0[i] *= scale_factor;
        zf_1[i] *= scale_factor;
        zf_2[i] *= scale_factor;
    }

    // save some flops by multiplying one vector with q
    for(int i=0; i<p; i++)
      z2_0[i] *= q;
    // FIXME/TODO: should we multiply also zf_0 with q?
    //             Or should we just remove this completely?

    return __IDX3_RMAJ(idx_from[0]+p_half,
                       idx_from[1]+p_half,
                       idx_from[2]+p_half,
                       params->npdims[1], params->npdims[2]);
}
#elif TWO_PERIODIC
static
int kaiser_expansion_2p_force(const double x[3], const double q,
                              const SE_FGG_params* params,
                              double z2_0[P_MAX],
                              double z2_1[P_MAX],
                              double z2_2[P_MAX],
                              double zf_0[P_MAX],
                              double zf_1[P_MAX],
                              double zf_2[P_MAX])
{
    // unpack params
    const int p = params->P;
    const int p_half = params->P_half;
    const double h = params->h;
    const double oh = 1./h;
    const double w = params->P/2.;
    const double ow2 = 1./(w*w);
    const double beta = params->beta;
    const double a = params->a;
    double t0[3];
    int shift_t0;
    double scale_factor;
    KaiserlikeWindow window = params->window;
    KaiserlikeWindow window_deriv = params->window_deriv;
    SE_window_params window_opt;
    if (params->use_polynomial_window) {
      window_opt.polynomial_degree = params->polynomial_degree;
      shift_t0 = 0;
      scale_factor = oh/w;
    } else {
      window_opt.ow2 = ow2;
      window_opt.beta = beta;
      shift_t0 = 1;
      scale_factor = oh;
    }

    int idx;
    int idx_from[3];

    // compute index range and centering
    if (is_odd(p)) {
      idx = (int) round(x[0]*oh);
      idx_from[0] = idx - p_half;
      t0[0] = x[0]*oh - idx + shift_t0*p_half;

      idx = (int) round(x[1]*oh);
      idx_from[1] = idx - p_half;
      t0[1] = x[1]*oh - idx + shift_t0*p_half;

      idx = (int) round((x[2]-(a+h/2))*oh);
      idx_from[2] = idx - p_half;
      t0[2] = (x[2]-(a+h/2))*oh - idx + shift_t0*p_half;
    }
    else {
      idx = (int) floor(x[0]*oh);
      idx_from[0] = idx - (p_half-1);
      t0[0] = x[0]*oh - idx + shift_t0*(p_half-1);

      idx = (int) floor(x[1]*oh);
      idx_from[1] = idx - (p_half-1);
      t0[1] = x[1]*oh - idx + shift_t0*(p_half-1);

      idx = (int) floor((x[2]-(a+h/2))*oh);
      idx_from[2] = idx - (p_half-1);
      t0[2] = (x[2]-(a+h/2))*oh - idx + shift_t0*(p_half-1);
    }

    // compute window function and its derivative
    window(t0, p, &window_opt, z2_0, z2_1, z2_2);
    window_deriv(t0, p, &window_opt, zf_0, zf_1, zf_2);

    // extra scale factors applied to derivative
    for (int i=0; i<p; i++)
    {
        zf_0[i] *= scale_factor;
        zf_1[i] *= scale_factor;
        zf_2[i] *= scale_factor;
    }

    // save some flops by multiplying one vector with q
    for (int i=0; i<p; i++)
      z2_0[i] *= q;
    // FIXME/TODO: should we multiply also zf_0 with q?
    //             Or should we just remove this completely?

    return __IDX3_RMAJ(idx_from[0]+p_half,
                       idx_from[1]+p_half,
                       idx_from[2],
                       params->npdims[1], params->npdims[2]);
}
#elif ONE_PERIODIC
static
int kaiser_expansion_1p_force(const double x[3], const double q,
                              const SE_FGG_params* params,
                              double z2_0[P_MAX],
                              double z2_1[P_MAX],
                              double z2_2[P_MAX],
                              double zf_0[P_MAX],
                              double zf_1[P_MAX],
                              double zf_2[P_MAX])
{
    // unpack params
    const int p = params->P;
    const int p_half = params->P_half;
    const double h = params->h;
    const double oh = 1./h;
    const double w = params->P/2.;
    const double ow2 = 1./(w*w);
    const double beta = params->beta;
    const double a = params->a;
    const double b = params->b;
    double t0[3];
    int shift_t0;
    double scale_factor;
    KaiserlikeWindow window = params->window;
    KaiserlikeWindow window_deriv = params->window_deriv;
    SE_window_params window_opt;
    if (params->use_polynomial_window) {
      window_opt.polynomial_degree = params->polynomial_degree;
      shift_t0 = 0;
      scale_factor = oh/w;
    } else {
      window_opt.ow2 = ow2;
      window_opt.beta = beta;
      shift_t0 = 1;
      scale_factor = oh;
    }

    int idx;
    int idx_from[3];

    // compute index range and centering
    if (is_odd(p)) {
      idx = (int) round(x[0]*oh);
      idx_from[0] = idx - p_half;
      t0[0] = x[0]*oh - idx + shift_t0*p_half;

      idx = (int) round((x[1]-(a+h/2))*oh);
      idx_from[1] = idx - p_half;
      t0[1] = (x[1]-(a+h/2))*oh - idx + shift_t0*p_half;

      idx = (int) round((x[2]-(b+h/2))*oh);
      idx_from[2] = idx - p_half;
      t0[2] = (x[2]-(b+h/2))*oh - idx + shift_t0*p_half;
    }
    else {
      idx = (int) floor(x[0]*oh);
      idx_from[0] = idx - (p_half-1);
      t0[0] = x[0]*oh - idx + shift_t0*(p_half-1);

      idx = (int) floor((x[1]-(a+h/2))*oh);
      idx_from[1] = idx - (p_half-1);
      t0[1] = (x[1]-(a+h/2))*oh - idx + shift_t0*(p_half-1);

      idx = (int) floor((x[2]-(b+h/2))*oh);
      idx_from[2] = idx - (p_half-1);
      t0[2] = (x[2]-(b+h/2))*oh - idx + shift_t0*(p_half-1);
    }

    // compute window function and its derivative
    window(t0, p, &window_opt, z2_0, z2_1, z2_2);
    window_deriv(t0, p, &window_opt, zf_0, zf_1, zf_2);

    // extra scale factors applied to derivative
    for (int i=0; i<p; i++)
    {
        zf_0[i] *= scale_factor;
        zf_1[i] *= scale_factor;
        zf_2[i] *= scale_factor;
    }

    // save some flops by multiplying one vector with q
    for (int i=0; i<p; i++)
      z2_0[i] *= q;
    // FIXME/TODO: should we multiply also zf_0 with q?
    //             Or should we just remove this completely?

    return __IDX3_RMAJ(idx_from[0]+p_half,
                       idx_from[1],
                       idx_from[2],
                       params->npdims[1], params->npdims[2]);
}
#endif

// -----------------------------------------------------------------------------
void SE_FKG_expand_all_force(SE_FGG_work* work,
                             const SE_state* st,
                             const SE_FGG_params* params)
{
    double xn[3] MEM_ALIGNED;
    const int N = params->N;
    const int Peval = params->P_eval_window;

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int n=0; n<N; n++)
    {
        // compute index and expansion vectors
        xn[0] = st->x[n]; xn[1] = st->x[n+N]; xn[2] = st->x[n+2*N];

        *(work->idx+n) = __FKG_EXPA_FORCE(xn,1,params,
                                          work->zx+n*Peval,
                                          work->zy+n*Peval,
                                          work->zz+n*Peval,
                                          work->zfx+n*Peval,
                                          work->zfy+n*Peval,
                                          work->zfz+n*Peval);
    }
}

// -----------------------------------------------------------------------------
// vanilla grid gather to calculate forces
void SE_FKG_int_force(double* restrict force,
                      const SE_FGG_work* work,
                      const SE_state* st,
                      const SE_FGG_params* params)
{
    double z2_0[P_MAX] MEM_ALIGNED;
    double z2_1[P_MAX] MEM_ALIGNED;
    double z2_2[P_MAX] MEM_ALIGNED;
    // to calculate forces:
    double zf_0[P_MAX] MEM_ALIGNED;
    double zf_1[P_MAX] MEM_ALIGNED;
    double zf_2[P_MAX] MEM_ALIGNED;

    // unpack params
    const double* restrict H = work->H;
    const int p = params->P;
    const int N = params->N;
    const double h = params->h;
    const double h3 = h*h*h;

    double xm[3];
    int i,j,k,idx;
    double force_m[3], aij, bij, cij, Hz;
#ifdef CALC_ENERGY
    double phi_m;
#endif

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++) {
        xm[0] = st->x[m]; xm[1] = st->x[m+N]; xm[2] = st->x[m+2*N];

        idx = __FKG_EXPA_FORCE(xm, 1, params, z2_0, z2_1, z2_2, zf_0, zf_1, zf_2);

        force_m[0] = 0; force_m[1] = 0; force_m[2] = 0;
#ifdef CALC_ENERGY
        phi_m = 0;
#endif

        for(i = 0; i<p; i++)
        {
            for(j = 0; j<p; j++)
            {
                aij = zf_0[i] * z2_1[j];
                bij = z2_0[i] * zf_1[j];
                cij = z2_0[i] * z2_1[j];
                for(k = 0; k<p; k++)
                {
                    Hz = H[idx] * z2_2[k];
#ifdef CALC_ENERGY
                    phi_m += Hz * cij;
#endif
                    force_m[0] += Hz * aij;
                    force_m[1] += Hz * bij;
                    force_m[2] += H[idx] * zf_2[k] * cij;
                    idx++;
                }
                idx += incrj;
            }
            idx += incri;
        }
#ifdef CALC_ENERGY
        st->phi[m] = h3*phi_m;
#endif
        force[m    ] = h3*force_m[0];
        force[m+  N] = h3*force_m[1];
        force[m+2*N] = h3*force_m[2];
    }
}

// -----------------------------------------------------------------------------
void SE_FKG_int_split_force(double* restrict force,
                            SE_state* st,
                            const SE_FGG_work* work,
                            const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    const double* restrict zfx = work->zfx;
    const double* restrict zfy = work->zfy;
    const double* restrict zfz = work->zfz;

    const int p = params->P;
    const int N = params->N;
    const double h = params->h;
    const double h3 = h*h*h;

    int i,j,k,idx,idx_zz;
    double force_m[3], aij, bij, cij, Hz;
#ifdef CALC_ENERGY
    double phi_m;
#endif

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
        idx = work->idx[m];
        force_m[0] = 0; force_m[1] = 0; force_m[2] = 0;
#ifdef CALC_ENERGY
        phi_m = 0;
#endif

        for(i = 0; i<p; i++)
        {
            for(j = 0; j<p; j++)
            {
                aij = zfx[m*p+i] * zy[m*p+j];
                bij = zx[m*p+i] * zfy[m*p+j];
                cij = zx[m*p+i] * zy[m*p+j];
                idx_zz = m*p;
                for(k = 0; k<p; k++)
                {
                    Hz = H[idx] * zz[idx_zz];
#ifdef CALC_ENERGY
                    phi_m += Hz * cij;
#endif
                    force_m[0] += Hz * aij;
                    force_m[1] += Hz * bij;
                    force_m[2] += H[idx] * zfz[idx_zz] * cij;
                    idx++; idx_zz++;
                }
                idx += incrj;
            }
            idx += incri;
        }
#ifdef CALC_ENERGY
        st->phi[m] = h3*phi_m;
#endif
        force[m    ] = h3*force_m[0];
        force[m+  N] = h3*force_m[1];
        force[m+2*N] = h3*force_m[2];
    }
}

#endif // end FORCE
