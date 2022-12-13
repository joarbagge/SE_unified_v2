add_definitions(-DMATLAB_MEX_FILE)
# add_definitions(-DVERBOSE)

## 3P
# 3P real
matlab_add_mex(
  NAME SE3P_Laplace_direct_real_mex
  SRC SE3P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE3P_Laplace_direct_real_mex PUBLIC THREE_PERIODIC EWALD_REAL)
matlab_add_mex(
  NAME SE3P_Laplace_direct_real_grad_mex
  SRC SE3P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE3P_Laplace_direct_real_grad_mex PUBLIC THREE_PERIODIC EWALD_REAL GRADIENT)

# 3P real_cut
matlab_add_mex(
  NAME SE3P_Laplace_direct_real_cut_mex
  SRC SE3P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE3P_Laplace_direct_real_cut_mex PUBLIC THREE_PERIODIC EWALD_REAL CUTOFF)
matlab_add_mex(
  NAME SE3P_Laplace_direct_real_cut_grad_mex
  SRC SE3P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE3P_Laplace_direct_real_cut_grad_mex PUBLIC THREE_PERIODIC EWALD_REAL CUTOFF GRADIENT)

# 3P real_cut_ext
matlab_add_mex(
  NAME SE3P_Laplace_direct_real_cut_ext_mex
  SRC SE3P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE3P_Laplace_direct_real_cut_ext_mex PUBLIC THREE_PERIODIC EWALD_REAL CUTOFF EXTERNAL)
matlab_add_mex(
  NAME SE3P_Laplace_direct_real_cut_ext_grad_mex
  SRC SE3P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE3P_Laplace_direct_real_cut_ext_grad_mex PUBLIC THREE_PERIODIC EWALD_REAL CUTOFF EXTERNAL GRADIENT)

# 3P fourier
matlab_add_mex(
  NAME SE3P_Laplace_direct_fourier_mex
  SRC SE3P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE3P_Laplace_direct_fourier_mex PUBLIC THREE_PERIODIC EWALD_FOURIER)
matlab_add_mex(
  NAME SE3P_Laplace_direct_fourier_grad_mex
  SRC SE3P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE3P_Laplace_direct_fourier_grad_mex PUBLIC THREE_PERIODIC EWALD_FOURIER GRADIENT)

# 3P fourier_ext
matlab_add_mex(
  NAME SE3P_Laplace_direct_fourier_ext_mex
  SRC SE3P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE3P_Laplace_direct_fourier_ext_mex PUBLIC THREE_PERIODIC EWALD_FOURIER EXTERNAL)
matlab_add_mex(
  NAME SE3P_Laplace_direct_fourier_ext_grad_mex
  SRC SE3P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE3P_Laplace_direct_fourier_ext_grad_mex PUBLIC THREE_PERIODIC EWALD_FOURIER EXTERNAL GRADIENT)

# 3P self
matlab_add_mex(
  NAME SE3P_Laplace_direct_self_mex
  SRC SE3P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE3P_Laplace_direct_self_mex PUBLIC THREE_PERIODIC EWALD_SELF)

## 2P
# 2P real
matlab_add_mex(
  NAME SE2P_Laplace_direct_real_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_real_mex PUBLIC TWO_PERIODIC EWALD_REAL)
matlab_add_mex(
  NAME SE2P_Laplace_direct_real_grad_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_real_grad_mex PUBLIC TWO_PERIODIC EWALD_REAL GRADIENT)

# 2P real_cut
matlab_add_mex(
  NAME SE2P_Laplace_direct_real_cut_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_real_cut_mex PUBLIC TWO_PERIODIC EWALD_REAL CUTOFF)
matlab_add_mex(
  NAME SE2P_Laplace_direct_real_cut_grad_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_real_cut_grad_mex PUBLIC TWO_PERIODIC EWALD_REAL CUTOFF GRADIENT)

# 2P real_cut_ext
matlab_add_mex(
  NAME SE2P_Laplace_direct_real_cut_ext_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_real_cut_ext_mex PUBLIC TWO_PERIODIC EWALD_REAL CUTOFF EXTERNAL)
matlab_add_mex(
  NAME SE2P_Laplace_direct_real_cut_ext_grad_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_real_cut_ext_grad_mex PUBLIC TWO_PERIODIC EWALD_REAL CUTOFF EXTERNAL GRADIENT)

# 2P fourier
matlab_add_mex(
  NAME SE2P_Laplace_direct_fourier_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_fourier_mex PUBLIC TWO_PERIODIC EWALD_FOURIER)
matlab_add_mex(
  NAME SE2P_Laplace_direct_fourier_grad_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_fourier_grad_mex PUBLIC TWO_PERIODIC EWALD_FOURIER GRADIENT)

# 2P fourier_ext
matlab_add_mex(
  NAME SE2P_Laplace_direct_fourier_ext_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_fourier_ext_mex PUBLIC TWO_PERIODIC EWALD_FOURIER EXTERNAL)
matlab_add_mex(
  NAME SE2P_Laplace_direct_fourier_ext_grad_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_fourier_ext_grad_mex PUBLIC TWO_PERIODIC EWALD_FOURIER EXTERNAL GRADIENT)

# 2P fourier_k0
matlab_add_mex(
  NAME SE2P_Laplace_direct_fourier_k0_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_fourier_k0_mex PUBLIC TWO_PERIODIC EWALD_FOURIER_K0)
matlab_add_mex(
  NAME SE2P_Laplace_direct_fourier_k0_grad_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_fourier_k0_grad_mex PUBLIC TWO_PERIODIC EWALD_FOURIER_K0 GRADIENT)

# 2P fourier_k0_ext
matlab_add_mex(
  NAME SE2P_Laplace_direct_fourier_k0_ext_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_fourier_k0_ext_mex PUBLIC TWO_PERIODIC EWALD_FOURIER_K0 EXTERNAL)
matlab_add_mex(
  NAME SE2P_Laplace_direct_fourier_k0_ext_grad_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_fourier_k0_ext_grad_mex PUBLIC TWO_PERIODIC EWALD_FOURIER_K0 EXTERNAL GRADIENT)

# 2P self
matlab_add_mex(
  NAME SE2P_Laplace_direct_self_mex
  SRC SE2P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE2P_Laplace_direct_self_mex PUBLIC TWO_PERIODIC EWALD_SELF)

## 1P
# 1P full
matlab_add_mex(
  NAME SE1P_Laplace_direct_full_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_full_mex PUBLIC ONE_PERIODIC EWALD_FULL)
matlab_add_mex(
  NAME SE1P_Laplace_direct_full_grad_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_full_grad_mex PUBLIC ONE_PERIODIC EWALD_FULL GRADIENT)

# 1P full_ext
matlab_add_mex(
  NAME SE1P_Laplace_direct_full_ext_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_full_ext_mex PUBLIC ONE_PERIODIC EWALD_FULL EXTERNAL)
matlab_add_mex(
  NAME SE1P_Laplace_direct_full_ext_grad_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_full_ext_grad_mex PUBLIC ONE_PERIODIC EWALD_FULL EXTERNAL GRADIENT)

# 1P real
matlab_add_mex(
  NAME SE1P_Laplace_direct_real_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_real_mex PUBLIC ONE_PERIODIC EWALD_REAL)
matlab_add_mex(
  NAME SE1P_Laplace_direct_real_grad_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_real_grad_mex PUBLIC ONE_PERIODIC EWALD_REAL GRADIENT)

# 1P real_cut
matlab_add_mex(
  NAME SE1P_Laplace_direct_real_cut_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_real_cut_mex PUBLIC ONE_PERIODIC EWALD_REAL CUTOFF)
matlab_add_mex(
  NAME SE1P_Laplace_direct_real_cut_grad_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_real_cut_grad_mex PUBLIC ONE_PERIODIC EWALD_REAL CUTOFF GRADIENT)

# 1P real_cut_ext
matlab_add_mex(
  NAME SE1P_Laplace_direct_real_cut_ext_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_real_cut_ext_mex PUBLIC ONE_PERIODIC EWALD_REAL CUTOFF EXTERNAL)
matlab_add_mex(
  NAME SE1P_Laplace_direct_real_cut_ext_grad_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_real_cut_ext_grad_mex PUBLIC ONE_PERIODIC EWALD_REAL CUTOFF EXTERNAL GRADIENT)

# 1P fourier
matlab_add_mex(
  NAME SE1P_Laplace_direct_fourier_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_fourier_mex PUBLIC ONE_PERIODIC EWALD_FOURIER)
target_link_libraries(SE1P_Laplace_direct_fourier_mex incomplete_bessel_K gsl gslcblas)
matlab_add_mex(
  NAME SE1P_Laplace_direct_fourier_grad_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_fourier_grad_mex PUBLIC ONE_PERIODIC EWALD_FOURIER GRADIENT)
target_link_libraries(SE1P_Laplace_direct_fourier_grad_mex incomplete_bessel_K gsl gslcblas)

# 1P fourier_ext
matlab_add_mex(
  NAME SE1P_Laplace_direct_fourier_ext_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_fourier_ext_mex PUBLIC ONE_PERIODIC EWALD_FOURIER EXTERNAL)
target_link_libraries(SE1P_Laplace_direct_fourier_ext_mex incomplete_bessel_K gsl gslcblas)
matlab_add_mex(
  NAME SE1P_Laplace_direct_fourier_ext_grad_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_fourier_ext_grad_mex PUBLIC ONE_PERIODIC EWALD_FOURIER EXTERNAL GRADIENT)
target_link_libraries(SE1P_Laplace_direct_fourier_ext_grad_mex incomplete_bessel_K gsl gslcblas)

# 1P fourier_k0
matlab_add_mex(
  NAME SE1P_Laplace_direct_fourier_k0_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_fourier_k0_mex PUBLIC ONE_PERIODIC EWALD_FOURIER_K0)
target_link_libraries(SE1P_Laplace_direct_fourier_k0_mex incomplete_bessel_K gsl gslcblas)
matlab_add_mex(
  NAME SE1P_Laplace_direct_fourier_k0_grad_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_fourier_k0_grad_mex PUBLIC ONE_PERIODIC EWALD_FOURIER_K0 GRADIENT)
target_link_libraries(SE1P_Laplace_direct_fourier_k0_grad_mex incomplete_bessel_K gsl gslcblas)

# 1P fourier_k0_ext
matlab_add_mex(
  NAME SE1P_Laplace_direct_fourier_k0_ext_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_fourier_k0_ext_mex PUBLIC ONE_PERIODIC EWALD_FOURIER_K0 EXTERNAL)
target_link_libraries(SE1P_Laplace_direct_fourier_k0_ext_mex incomplete_bessel_K gsl gslcblas)
matlab_add_mex(
  NAME SE1P_Laplace_direct_fourier_k0_ext_grad_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_fourier_k0_ext_grad_mex PUBLIC ONE_PERIODIC EWALD_FOURIER_K0 EXTERNAL GRADIENT)
target_link_libraries(SE1P_Laplace_direct_fourier_k0_ext_grad_mex incomplete_bessel_K gsl gslcblas)

# 1P self
matlab_add_mex(
  NAME SE1P_Laplace_direct_self_mex
  SRC SE1P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE1P_Laplace_direct_self_mex PUBLIC ONE_PERIODIC EWALD_SELF)

## 0P
# 0P full
matlab_add_mex(
  NAME SE0P_Laplace_direct_full_mex
  SRC SE0P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE0P_Laplace_direct_full_mex PUBLIC ZERO_PERIODIC EWALD_FULL)
matlab_add_mex(
  NAME SE0P_Laplace_direct_full_grad_mex
  SRC SE0P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE0P_Laplace_direct_full_grad_mex PUBLIC ZERO_PERIODIC EWALD_FULL GRADIENT)

# 0P full_ext
matlab_add_mex(
  NAME SE0P_Laplace_direct_full_ext_mex
  SRC SE0P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE0P_Laplace_direct_full_ext_mex PUBLIC ZERO_PERIODIC EWALD_FULL EXTERNAL)
matlab_add_mex(
  NAME SE0P_Laplace_direct_full_ext_grad_mex
  SRC SE0P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE0P_Laplace_direct_full_ext_grad_mex PUBLIC ZERO_PERIODIC EWALD_FULL EXTERNAL GRADIENT)

# 0P real
matlab_add_mex(
  NAME SE0P_Laplace_direct_real_mex
  SRC SE0P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE0P_Laplace_direct_real_mex PUBLIC ZERO_PERIODIC EWALD_REAL)
matlab_add_mex(
  NAME SE0P_Laplace_direct_real_grad_mex
  SRC SE0P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE0P_Laplace_direct_real_grad_mex PUBLIC ZERO_PERIODIC EWALD_REAL GRADIENT)

# 0P real_cut
matlab_add_mex(
  NAME SE0P_Laplace_direct_real_cut_mex
  SRC SE0P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE0P_Laplace_direct_real_cut_mex PUBLIC ZERO_PERIODIC EWALD_REAL CUTOFF)
matlab_add_mex(
  NAME SE0P_Laplace_direct_real_cut_grad_mex
  SRC SE0P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE0P_Laplace_direct_real_cut_grad_mex PUBLIC ZERO_PERIODIC EWALD_REAL CUTOFF GRADIENT)

# 0P real_cut_ext
matlab_add_mex(
  NAME SE0P_Laplace_direct_real_cut_ext_mex
  SRC SE0P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE0P_Laplace_direct_real_cut_ext_mex PUBLIC ZERO_PERIODIC EWALD_REAL CUTOFF EXTERNAL)
matlab_add_mex(
  NAME SE0P_Laplace_direct_real_cut_ext_grad_mex
  SRC SE0P_Laplace_direct.c SE_Laplace_direct_mex.c
)
target_compile_definitions(SE0P_Laplace_direct_real_cut_ext_grad_mex PUBLIC ZERO_PERIODIC EWALD_REAL CUTOFF EXTERNAL GRADIENT)
