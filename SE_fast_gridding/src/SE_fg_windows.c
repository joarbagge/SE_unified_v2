#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include "SE_fg_windows.h"

/**
 * Exponential semicircle window
 */
static inline double _basic_expsemicirc(const double x, const double ow2,
                                        const double beta) {
  const double t = sqrt(1. - x*x*ow2);
  return exp(beta*(t-1));
}

void SE_fg_window_expsemicirc(const double t0[3], const int P,
                              const SE_window_params * opt,
                              double out1[P_MAX],
                              double out2[P_MAX],
                              double out3[P_MAX]) {
  const double ow2 = opt->ow2;
  const double beta = opt->beta;
  for (int i=0; i<P; i++) {
    out1[i] = _basic_expsemicirc(t0[0]-i, ow2, beta);
    out2[i] = _basic_expsemicirc(t0[1]-i, ow2, beta);
    out3[i] = _basic_expsemicirc(t0[2]-i, ow2, beta);
  }
}

/**
 * Exact Kaiser-Bessel window
 */
static inline double _basic_kaiser_exact(const double x, const double ow2,
                                         const double beta) {
  const double t = sqrt(1. - x*x*ow2);
  return gsl_sf_bessel_I0(beta*t) / gsl_sf_bessel_I0(beta);
}

static inline double _basic_kaiser_exact_deriv(const double x, const double ow2,
                                               const double beta) {
  const double t = sqrt(1. - x*x*ow2);
  if (t < __DBL_EPSILON__)
  {
    return 0; // limit
  }
  else
  {
    return -beta*x*ow2 * gsl_sf_bessel_I1(beta*t) / (gsl_sf_bessel_I0(beta) * t);
  }
}

void SE_fg_window_kaiser_exact(const double t0[3], const int P,
                               const SE_window_params * opt,
                               double out1[P_MAX],
                               double out2[P_MAX],
                               double out3[P_MAX]) {
  const double ow2 = opt->ow2;
  const double beta = opt->beta;
  for (int i=0; i<P; i++) {
    out1[i] = _basic_kaiser_exact(t0[0]-i, ow2, beta);
    out2[i] = _basic_kaiser_exact(t0[1]-i, ow2, beta);
    out3[i] = _basic_kaiser_exact(t0[2]-i, ow2, beta);
  }
}

void SE_fg_window_kaiser_exact_deriv(const double t0[3], const int P,
                                     const SE_window_params * opt,
                                     double out1[P_MAX],
                                     double out2[P_MAX],
                                     double out3[P_MAX]) {
  const double ow2 = opt->ow2;
  const double beta = opt->beta;
  for (int i=0; i<P; i++) {
    out1[i] = _basic_kaiser_exact_deriv(t0[0]-i, ow2, beta);
    out2[i] = _basic_kaiser_exact_deriv(t0[1]-i, ow2, beta);
    out3[i] = _basic_kaiser_exact_deriv(t0[2]-i, ow2, beta);
  }
}

/**
 * Polynomial approximation of Kaiser-Bessel window
 */
static inline void _basic_kaiser_poly(const double t0, const int P,
                                      const int degree,
                                      double out[P_MAX]) {
  const double z = 2*t0 - 1; // scale so that local grid offset z is in [-1,1]
  //printf("_basic_kaiser_poly: z=%g, P=%d, degree=%d\n", z, P, degree);
  // insert generated code which expects z, P, degree and writes to out
#include "gen_basic_kaiser_poly.c"
}

static inline void _basic_kaiser_poly_deriv(const double t0, const int P,
                                            const int degree,
                                            double out[P_MAX]) {
  const double z = 2*t0 - 1; // scale so that local grid offset z is in [-1,1]
  //printf("_basic_kaiser_poly: z=%g, P=%d, degree=%d\n", z, P, degree);
  // insert generated code which expects z, P, degree and writes to out
#include "gen_basic_kaiser_poly_deriv.c"
}

void SE_fg_window_kaiser_poly(const double t0[3], const int P,
                              const SE_window_params * opt,
                              double out1[P_MAX],
                              double out2[P_MAX],
                              double out3[P_MAX]) {
  const int degree = opt->polynomial_degree;
  _basic_kaiser_poly(t0[0], P, degree, out1);
  _basic_kaiser_poly(t0[1], P, degree, out2);
  _basic_kaiser_poly(t0[2], P, degree, out3);
}

void SE_fg_window_kaiser_poly_deriv(const double t0[3], const int P,
                                    const SE_window_params * opt,
                                    double out1[P_MAX],
                                    double out2[P_MAX],
                                    double out3[P_MAX]) {
  const int degree = opt->polynomial_degree;
  _basic_kaiser_poly_deriv(t0[0], P, degree, out1);
  _basic_kaiser_poly_deriv(t0[1], P, degree, out2);
  _basic_kaiser_poly_deriv(t0[2], P, degree, out3);
}
