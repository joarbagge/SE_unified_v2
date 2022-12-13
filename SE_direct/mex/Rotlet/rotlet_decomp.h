#ifndef SE_ROTLET_DECOMP_H
#define SE_ROTLET_DECOMP_H

#include <math.h>

inline double dot(const double * a, const double * b)
{
    return ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

/* Computes contribution of u += A cross(t, r) */
void real_part_rotlet_contribution(double u[3],
                                   const double r[3],
                                   const double t[3],
                                   const double xi)
{
    const double R2 = dot(r, r);
    const double R = sqrt(R2);
    const double A = (2*xi*exp(-xi*xi*R2)/sqrt(PI) + erfc(xi*R)/R) / R2;
    u[0] += A * (t[1]*r[2] - t[2]*r[1]);
    u[1] += A * (t[2]*r[0] - t[0]*r[2]);
    u[2] += A * (t[0]*r[1] - t[1]*r[0]);
}

/* Computes contribution of u += B cross(f, k) */
void fourier_part_rotlet_contribution(double u[3],
                                      const double k[3],
                                      const double f[3],
                                      const double xi)
{
    const double k2 = dot(k, k);
    const double B = -4*PI * exp(-k2/(4*xi*xi)) / k2;
    u[0] += B * (f[1]*k[2] - f[2]*k[1]);
    u[1] += B * (f[2]*k[0] - f[0]*k[2]);
    u[2] += B * (f[0]*k[1] - f[1]*k[0]);
}

#endif
