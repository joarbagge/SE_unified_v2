#include "mex.h"
#include "../../SE3P/mex/mex_compat.h"
#include "math.h"
#include "../../SE3P/mex/cell_list.h"

#ifdef INTEL_MKL
#include "mkl.h"
#endif

// CODE_VARIANT can be 1 or 2
// (Variant 1 seems to be faster, but can only do 1 layer)
#define CODE_VARIANT 1

#define X   prhs[0] // source locations (Nx3)
#define F   prhs[1] // source strengths (Nx1)
#define RC  prhs[2] // cutoff radius
#define XI  prhs[3] // Ewald decomposition parameter
#if CODE_VARIANT == 1
#define BOX prhs[4] // domain size (1x3)
#endif
#if CODE_VARIANT == 2
#define P   prhs[4] // number of periodic wrap layers
#define BOX prhs[5] // domain size (1x3)
#endif

#define OUT plhs[0]  // Output (Nx1 or N×3 depending on GRADIENT)

#ifndef VERBOSE
#define VERBOSE 0
#endif

#ifdef _OPENMP
#define CRITICAL _Pragma("omp critical")
#else
#define CRITICAL
#endif

#define PI 3.141592653589793

inline double norm2(double * a)
{
    return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
}

// Transpose matrix
static void transpose(const double* restrict in, double* restrict out, const int N)
{
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<3; j++)
        {
            out[i*3+j] = in[i+j*N];
        }
    }
}

// Compute buffer stuff
#define BUF_SIZE 256
typedef struct {
    int n;
    int idx_t[BUF_SIZE];
    double rvec[3*BUF_SIZE];
    double r2[BUF_SIZE];
    double a;
} ComputeBuffer;

static void buffer_push(ComputeBuffer* buffer, int idx_t, double* rvec, double r2)
{
    int n = buffer->n;
    buffer->idx_t[n] = idx_t;
    for(int i=0; i<3; i++)
    {
        buffer->rvec[3*n + i] = rvec[i];
    }
    buffer->r2[n] = r2;
    buffer->n = n + 1;
}

static void empty_buffer(ComputeBuffer* buffer,
                         const double* restrict x,
                         const double* restrict f,
                         double* restrict u,
                         int idx_s,
                         double xi,
                         const int N)
{
    int nbuf = buffer->n;
    double fs, ft;
#ifndef GRADIENT
    double us = 0;
#else
    double us[] = {0, 0, 0};
    double xi2 = xi*xi;
#endif
    fs = f[idx_s];

    // Do what we can to help the compiler vectorize exp and erfc, if possible
#ifdef GRADIENT
    const double* restrict rvec = buffer->rvec;
    const double a = buffer->a;
#endif
    const double* restrict r2 = buffer->r2;
    double c1[BUF_SIZE];
// TODO: The following section needs to be modified if INTEL_MKL
// and GRADIENT are to work simultaneously.
#ifdef INTEL_MKL
    double r[BUF_SIZE];
    double xir[BUF_SIZE];
    double xi2r2[BUF_SIZE];
#pragma ivdep
    for (int n=0; n<nbuf; n++)
    {
        r[n] = sqrt(r2[n]);
        xir[n] = xi*r[n];
        xi2r2[n] = -xi2*r2[n];
    }
    double erfc_vec[BUF_SIZE];
    double exp_vec[BUF_SIZE];
    vdErfc(nbuf, xir, erfc_vec);
    vdExp(nbuf, xi2r2, exp_vec);
#pragma ivdep
    for (int n=0; n<nbuf; n++)
    {
        double xiexp = xi*exp_vec[n];
        c1[n] = erfc_vec[n] / r[n];
    }
#else
    for (int n=0; n<nbuf; n++)
    {
        double r = sqrt(r2[n]);
#ifndef GRADIENT
        c1[n] = erfc(xi*r) / r;
#else
        c1[n] = -( erfc(xi*r) / r + a*exp(-xi2*r2[n]) )/r2[n];
#endif
    }
#endif
    // Compute interactions
#ifdef INTEL_MKL
#pragma ivdep
#endif
    for (int n=0; n<nbuf; n++)
    {
        int idx_t = buffer->idx_t[n];
        ft = f[idx_t];
#ifndef GRADIENT
        u[idx_t] += fs*c1[n];
        us += ft*c1[n];
#else
        // rvec points *from* t *to* s, sign is important here
        u[idx_t*3    ] -= fs*c1[n] * rvec[n*3    ];
        u[idx_t*3 + 1] -= fs*c1[n] * rvec[n*3 + 1];
        u[idx_t*3 + 2] -= fs*c1[n] * rvec[n*3 + 2];
        us[0] += ft*c1[n] * rvec[n*3    ];
        us[1] += ft*c1[n] * rvec[n*3 + 1];
        us[2] += ft*c1[n] * rvec[n*3 + 2];
#endif
    }
#ifndef GRADIENT
    u[idx_s] += us;
#else
    u[idx_s*3    ] += us[0];
    u[idx_s*3 + 1] += us[1];
    u[idx_s*3 + 2] += us[2];
#endif
    buffer->n = 0;
}

// Entry point
void
mexFunction( int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[] )
{
    // Input
    const int N = mxGetM(X);
    const double xi = mxGetScalar(XI);
    const double rc = mxGetScalar(RC);
    const double rc2 = rc*rc;
#if CODE_VARIANT == 2
    const double p = mxGetScalar(P);
#endif
    const double* restrict x_in = mxGetPr(X);
    const double* restrict f_tmp = mxGetPr(F);
    const double* restrict box = mxGetPr(BOX);

    // Transpose input for better memory access
    double* restrict x_tmp = __MALLOC(3*N*sizeof(double));
    transpose(x_in, x_tmp, N);

    // Setup cell list variables
    double rn[3];
    int ncell[3];
    int* restrict cell_list;
    int* restrict cell_idx;
    int px[27] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    int py[27] = {-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    int pz[27] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    // Build cell list
    build_cell_list(x_tmp, N, box, rc, rn, ncell, &cell_list, &cell_idx);
    if (cell_list == NULL)
    {
        mexErrMsgIdAndTxt("SE:RealCutCell", "Could not build cell list (is rc too large?)");
    }

    // Sort the arrays "x" and "f" according to cell list for
    // faster memory access
    double* restrict x = __MALLOC(3*N*sizeof(double));
    double* restrict f = __MALLOC(  N*sizeof(double));
    for (int n=0; n<N; n++)
    {
        for (int i=0; i<3; i++)
        {
            x[3*n+i] = x_tmp[3*cell_list[n]+i];
        }
        f[n] = f_tmp[cell_list[n]];
    }
    __FREE(x_tmp);

    // Output
#ifndef GRADIENT
#define OUT_DIM_2 1
    OUT = mxCreateDoubleMatrix(N, OUT_DIM_2, mxREAL);
    double* restrict pot_out = mxGetPr(OUT);
    double* restrict pot_tmp = __CALLOC(OUT_DIM_2*N, sizeof(double));
#else
#define OUT_DIM_2 3
    OUT = mxCreateDoubleMatrix(N, OUT_DIM_2, mxREAL);
    double* restrict grad_out = mxGetPr(OUT);
    double* restrict grad_tmp = __CALLOC(OUT_DIM_2*N, sizeof(double));
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    { // Begin parallel section
        // Create a thread-local compute buffer
        ComputeBuffer buffer;
        buffer.a = 2*xi/sqrt(PI);
        // Setup local output
        double* restrict u;
        CRITICAL {
            u = __CALLOC(OUT_DIM_2*N, sizeof(double));
        }

        // Main loop
// TODO: schedule(dynamic) didn't work for some reason...
#ifdef _OPENMP
#pragma omp for nowait
#endif
        // Loop over all points (work-shared)
        for (int idx_s=0; idx_s<N; idx_s++)
        {
            double xs[3];
            int home_cell[3], icell[3];
            for(int i=0; i<3; i++)
            {
                // Source point
                xs[i] = x[idx_s*3 + i];
                // Determine home cell
                home_cell[i] = xs[i]/rn[i];
            }

            // Loop over neighbouring cells (including home cell)
            buffer.n = 0;
            for(int ip=0; ip<27; ip++)
            {
                // Get neighbour cell
                icell[0] = home_cell[0] + px[ip];
                icell[1] = home_cell[1] + py[ip];
                icell[2] = home_cell[2] + pz[ip];
#if CODE_VARIANT == 1
                // Periodic wrap in periodic directions
                double pshift[3] = {0,0,0};
                bool cell_out_of_bounds = false;
                for (int j=0; j<2; j++)
                {
                    // (Could do this with mod)
                    if (icell[j] >= ncell[j])
                    {
                        icell[j] = 0;
                        pshift[j] = box[j];
                    }
                    else if (icell[j] < 0)
                    {
                        icell[j] = ncell[j] - 1;
                        pshift[j] = -box[j];
                    }
                }
                if (icell[2] >= ncell[2] || icell[2] < 0)
                {
                    cell_out_of_bounds = true;
                }
                if (cell_out_of_bounds)
                {
                    continue;
                }
#endif
#if CODE_VARIANT == 2
                // Stop at boundaries
                int inside = 1;
                for(int j=0; j<3; j++)
                {
                    if (icell[j] < 0 || icell[j] == ncell[j])
                        inside = 0;
                }
                if (!inside)
                    continue;
#endif
                int icell_idx =
                    icell[0] +
                    icell[1]*ncell[0] +
                    icell[2]*ncell[1]*ncell[0];
                // Extract indices (for cell_list) of the particles in this cell
                int cell_a = cell_idx[icell_idx];
                int cell_b = cell_idx[icell_idx+1];
                // Loop over particles in this cell
                for (int idx_t=cell_a; idx_t<cell_b; idx_t++) // particle index (for x, f)
                {
                    if (idx_s >= idx_t) continue;
                    double rvec[3];
#if CODE_VARIANT == 2
                    // Periodic wrap
                    for (int j0=-p; j0<=p; j0++)
                    {
                        for (int j1=-p; j1<=p; j1++)
                        {
                            double pshift[] = {j0*box[0], j1*box[1], 0};
#endif
                            for (int i=0; i<3; i++)
                                rvec[i] = xs[i] - x[idx_t*3 + i] - pshift[i];
                            double r2 = norm2(rvec);
                            if (r2 > rc2) continue;
                            buffer_push(&buffer, idx_t, rvec, r2);
                            if (buffer.n == BUF_SIZE)
                                empty_buffer(&buffer, x, f, u, idx_s, xi, N);
#if CODE_VARIANT == 2
                        }
                    }
#endif
                } // End loop over particles in this cell
            } // End loop over neighbouring cells
            empty_buffer(&buffer, x, f, u, idx_s, xi, N);
        } // End loop over all points
        // Collect results
        CRITICAL {
            for (int n=0; n<N; n++)
            {
#ifndef GRADIENT
                pot_tmp[n] += u[n];
#else
                grad_tmp[3*n  ] += u[3*n  ];
                grad_tmp[3*n+1] += u[3*n+1];
                grad_tmp[3*n+2] += u[3*n+2];
#endif
            }
        }
        // free/malloc not thread safe under MEX
        CRITICAL {
            __FREE(u);
        }
    } // End of parallel section
    __FREE(x);
    __FREE(f);
#ifndef GRADIENT
    // Sort u to agree with input order
    for (int n=0; n<N; n++)
    {
        pot_out[cell_list[n]] = pot_tmp[n];
    }
    __FREE(pot_tmp);
#else
    // Sort u to agree with input order, and transpose
    for (int n=0; n<N; n++)
    {
        for (int i=0; i<3; i++)
        {
            grad_out[cell_list[n]+i*N] = grad_tmp[3*n+i];
        }
    }
    __FREE(grad_tmp);
#endif
    __FREE(cell_list);
    __FREE(cell_idx);
}
