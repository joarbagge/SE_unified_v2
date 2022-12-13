#ifndef SE_STRESSLET_HASIMOTO_DECOMP_H
#define SE_STRESSLET_HASIMOTO_DECOMP_H

#include <math.h>

inline double dot(const double * a, const double * b)
{
    return ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

/* Computes contribution of u_j += TR_jlm(r,xi) q_l n_m */
void real_part_stresslet_contribution(double u[3],
                                      const double r[3],
                                      const double q[3],
                                      const double n[3],
                                      const double xi)
{
    const double r_dot_q = dot(r, q);
    const double r_dot_n = dot(r, n);
    const double q_dot_n = dot(q, n);
    const double R2 = dot(r, r);
    const double R = sqrt(R2);
    const double xi2 = xi*xi;
    const double c = xi2*R2;
    const double expc = exp(-c);

    const double A = -2/(R2*R2) * (3*erfc(xi*R)/R + 2*xi/sqrt(PI)*(3+2*c)*expc);
    const double B = 4*xi*xi2/sqrt(PI) * expc;
    const double Kr = A * r_dot_q * r_dot_n + B * q_dot_n;
    const double Kq = B * r_dot_n;
    const double Kn = B * r_dot_q;
    for (int j=0; j<3; j++)
    {
        u[j] += Kr*r[j] + Kq*q[j] + Kn*n[j];
    }
}

/* Computes contribution of u_j += -i TF_jlm(k,xi) f_lm */
void fourier_part_stresslet_contribution(double u[3],
                                         const double k[3],
                                         const double f[3][3],
                                         const double xi)
{
    const double k2 = dot(k, k);
    const double xi2 = xi*xi;
    const double c = k2/(4*xi2);
    const double gammaB = exp(-c)*(1+c) * (-8*PI)/(k2*k2);
    double du[] = {0, 0, 0};
    for (int j=0; j<3; j++)
    {
        for (int m=0; m<3; m++)
        {
            du[j] += -(k[m]*f[j][m] + k[m]*f[m][j] + k[j]*f[m][m])*k2;
            for (int l=0; l<3; l++)
                du[j] += 2*k[j]*k[l]*k[m]*f[l][m];
        }
        u[j] += gammaB * du[j];
    }
}

#endif
