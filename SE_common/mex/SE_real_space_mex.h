typedef struct
{
    int periodicity;
    double box[3];
    double xi;
    double rc;
} Options;

typedef struct
{
    double (*fun)(const double); // function pointer
    double (*taylor)(const double, const double, const int); // function pointer
    double threshold;
    int terms;
} Smoothing;

typedef
void (*kernel_fun_t)(double* restrict, const double, const double, const double* restrict,
         const double, const double, const double* restrict, const double* restrict,
         const double, const double, const double, const double, const double,
         const Smoothing*);

int ALL_SHIFTS[][3] = {{0,0,0}, // 0P
                       {-1,0,0}, {1,0,0}, // 1P
                       {0,-1,0}, {-1,-1,0}, {1,-1,0}, {0,1,0}, {-1,1,0}, {1,1,0}, // 2P
                       {0,0,-1}, {-1,0,-1}, {1,0,-1}, {0,-1,-1}, {-1,-1,-1}, {1,-1,-1},
                       {0,1,-1}, {-1,1,-1}, {1,1,-1},
                       {0,0,1}, {-1,0,1}, {1,0,1}, {0,-1,1}, {-1,-1,1}, {1,-1,1},
                       {0,1,1}, {-1,1,1}, {1,1,1}}; // 3P

double inline dot3(const double* restrict x, const double* restrict y)
{
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

void load_options(const mxArray* s, Options* restrict opt);
int get_num_shifts(const int periodicity);

void core_computation(double* restrict pot, kernel_fun_t kernel_fun, bool double_layer,
                      const Options* restrict opt, const double* restrict targets,
                      size_t Nt, const double* restrict sources, size_t Ns,
                      const double* restrict normals, const double* restrict density,
                      const double* restrict delta, size_t numel_delta,
                      const Smoothing* restrict smooth, int num_shifts);

void stokes_sl_kernel(double* restrict pot, const double xi, const double xi2,
                      const double* restrict rvec, const double r, const double r2,
                      const double* restrict n, const double* restrict f,
                      const double rdotf, const double delt,
                      const double A, const double B, const double C,
                      const Smoothing* smooth);
void stokes_dl_kernel(double* restrict pot, const double xi, const double xi2,
                      const double* restrict rvec, const double r, const double r2,
                      const double* restrict n, const double* restrict f,
                      const double rdotf, const double delt,
                      const double A, const double B, const double C,
                      const Smoothing* smooth);

double smoothfun1(const double x);
double smoothfun1_taylor(const double d, const double rd, const int terms);
double get_smoothfun1_threshold(const int terms);
double smoothfun2(const double x);
double smoothfun2_taylor(const double d, const double rd, const int terms);
double get_smoothfun2_threshold(const int terms);
double smoothfun3_full(const double x);
double smoothfun3_full_taylor(const double d, const double rd, const int terms);
double get_smoothfun3_full_threshold(const int terms);
double smoothfun3_simple(const double x);
double smoothfun3_simple_taylor(const double d, const double rd, const int terms);
double get_smoothfun3_simple_threshold(const int terms);

// vim:set shiftwidth=4 softtabstop=4:
