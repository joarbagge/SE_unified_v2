#ifndef __SE_FG_WINDOWS_H
#define __SE_FG_WINDOWS_H

// Maximal size of window support (defined to help the compiler)
#define P_MAX 32

// Window parameters
typedef struct {
  // For exact Kaiser-like windows
  double ow2;
  double beta;
  // For polynomial windows
  int polynomial_degree;
} SE_window_params;

// Signature of Kaiser-like windows
typedef void (*KaiserlikeWindow)(const double t0[3], const int P,
                                 const SE_window_params * opt,
                                 double out1[P_MAX],
                                 double out2[P_MAX],
                                 double out3[P_MAX]);

#endif
