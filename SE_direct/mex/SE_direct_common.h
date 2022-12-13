#ifndef SE_DIRECT_COMMON_H
#define SE_DIRECT_COMMON_H

#include <stddef.h>

#define PI 3.141592653589793

typedef struct
{
    double box[3]; // size of periodic cell [L1,L2,L3]
    double xi; // Ewald decomposition parameter
    double rc; // cutoff radius
    size_t layers; // controls the number of replications of the
                   // periodic cell (real space), or the maximum
                   // wavenumber (Fourier space)
    size_t N_eval; // number of source points to evaluate at
    size_t * eval_idx; // indices of source points to evaluate at
    size_t N_ext; // number of external points to evaluate at
    double * eval_ext_x; // external points to evaluate at
    double stokeslet_k0_constant; // constant in the 2D biharmonic Green's function
} ewald_opts;

#endif
