add_definitions(-DMATLAB_MEX_FILE)
# add_definitions(-DVERBOSE)

## 3P
# 3P real
matlab_add_mex(
  NAME SE3P_Rotlet_direct_real_mex
  SRC SE3P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE3P_Rotlet_direct_real_mex PUBLIC THREE_PERIODIC EWALD_REAL)

# 3P real_cut
matlab_add_mex(
  NAME SE3P_Rotlet_direct_real_cut_mex
  SRC SE3P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE3P_Rotlet_direct_real_cut_mex PUBLIC THREE_PERIODIC EWALD_REAL CUTOFF)

# 3P real_cut_ext
matlab_add_mex(
  NAME SE3P_Rotlet_direct_real_cut_ext_mex
  SRC SE3P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE3P_Rotlet_direct_real_cut_ext_mex PUBLIC THREE_PERIODIC EWALD_REAL CUTOFF EXTERNAL)

# 3P fourier
matlab_add_mex(
  NAME SE3P_Rotlet_direct_fourier_mex
  SRC SE3P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE3P_Rotlet_direct_fourier_mex PUBLIC THREE_PERIODIC EWALD_FOURIER)

# 3P fourier_ext
matlab_add_mex(
  NAME SE3P_Rotlet_direct_fourier_ext_mex
  SRC SE3P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE3P_Rotlet_direct_fourier_ext_mex PUBLIC THREE_PERIODIC EWALD_FOURIER EXTERNAL)

## 2P
# 2P real
matlab_add_mex(
  NAME SE2P_Rotlet_direct_real_mex
  SRC SE2P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE2P_Rotlet_direct_real_mex PUBLIC TWO_PERIODIC EWALD_REAL)

# 2P real_cut
matlab_add_mex(
  NAME SE2P_Rotlet_direct_real_cut_mex
  SRC SE2P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE2P_Rotlet_direct_real_cut_mex PUBLIC TWO_PERIODIC EWALD_REAL CUTOFF)

# 2P real_cut_ext
matlab_add_mex(
  NAME SE2P_Rotlet_direct_real_cut_ext_mex
  SRC SE2P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE2P_Rotlet_direct_real_cut_ext_mex PUBLIC TWO_PERIODIC EWALD_REAL CUTOFF EXTERNAL)

# 2P fourier
matlab_add_mex(
  NAME SE2P_Rotlet_direct_fourier_mex
  SRC SE2P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE2P_Rotlet_direct_fourier_mex PUBLIC TWO_PERIODIC EWALD_FOURIER)

# 2P fourier_ext
matlab_add_mex(
  NAME SE2P_Rotlet_direct_fourier_ext_mex
  SRC SE2P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE2P_Rotlet_direct_fourier_ext_mex PUBLIC TWO_PERIODIC EWALD_FOURIER EXTERNAL)

# 2P fourier_k0
matlab_add_mex(
  NAME SE2P_Rotlet_direct_fourier_k0_mex
  SRC SE2P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE2P_Rotlet_direct_fourier_k0_mex PUBLIC TWO_PERIODIC EWALD_FOURIER_K0)

# 2P fourier_k0_ext
matlab_add_mex(
  NAME SE2P_Rotlet_direct_fourier_k0_ext_mex
  SRC SE2P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE2P_Rotlet_direct_fourier_k0_ext_mex PUBLIC TWO_PERIODIC EWALD_FOURIER_K0 EXTERNAL)

## 1P
# 1P full
matlab_add_mex(
  NAME SE1P_Rotlet_direct_full_mex
  SRC SE1P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE1P_Rotlet_direct_full_mex PUBLIC ONE_PERIODIC EWALD_FULL)

# 1P full_ext
matlab_add_mex(
  NAME SE1P_Rotlet_direct_full_ext_mex
  SRC SE1P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE1P_Rotlet_direct_full_ext_mex PUBLIC ONE_PERIODIC EWALD_FULL EXTERNAL)

# 1P real
matlab_add_mex(
  NAME SE1P_Rotlet_direct_real_mex
  SRC SE1P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE1P_Rotlet_direct_real_mex PUBLIC ONE_PERIODIC EWALD_REAL)

# 1P real_cut
matlab_add_mex(
  NAME SE1P_Rotlet_direct_real_cut_mex
  SRC SE1P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE1P_Rotlet_direct_real_cut_mex PUBLIC ONE_PERIODIC EWALD_REAL CUTOFF)

# 1P real_cut_ext
matlab_add_mex(
  NAME SE1P_Rotlet_direct_real_cut_ext_mex
  SRC SE1P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE1P_Rotlet_direct_real_cut_ext_mex PUBLIC ONE_PERIODIC EWALD_REAL CUTOFF EXTERNAL)

# 1P fourier
matlab_add_mex(
  NAME SE1P_Rotlet_direct_fourier_mex
  SRC SE1P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE1P_Rotlet_direct_fourier_mex PUBLIC ONE_PERIODIC EWALD_FOURIER)
target_link_libraries(SE1P_Rotlet_direct_fourier_mex incomplete_bessel_K gsl gslcblas)

# 1P fourier_ext
matlab_add_mex(
  NAME SE1P_Rotlet_direct_fourier_ext_mex
  SRC SE1P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE1P_Rotlet_direct_fourier_ext_mex PUBLIC ONE_PERIODIC EWALD_FOURIER EXTERNAL)
target_link_libraries(SE1P_Rotlet_direct_fourier_ext_mex incomplete_bessel_K gsl gslcblas)

# 1P fourier_k0
matlab_add_mex(
  NAME SE1P_Rotlet_direct_fourier_k0_mex
  SRC SE1P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE1P_Rotlet_direct_fourier_k0_mex PUBLIC ONE_PERIODIC EWALD_FOURIER_K0)
target_link_libraries(SE1P_Rotlet_direct_fourier_k0_mex incomplete_bessel_K gsl gslcblas)

# 1P fourier_k0_ext
matlab_add_mex(
  NAME SE1P_Rotlet_direct_fourier_k0_ext_mex
  SRC SE1P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE1P_Rotlet_direct_fourier_k0_ext_mex PUBLIC ONE_PERIODIC EWALD_FOURIER_K0 EXTERNAL)
target_link_libraries(SE1P_Rotlet_direct_fourier_k0_ext_mex incomplete_bessel_K gsl gslcblas)

## 0P
# 0P full
matlab_add_mex(
  NAME SE0P_Rotlet_direct_full_mex
  SRC SE0P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE0P_Rotlet_direct_full_mex PUBLIC ZERO_PERIODIC EWALD_FULL)

# 0P full_ext
matlab_add_mex(
  NAME SE0P_Rotlet_direct_full_ext_mex
  SRC SE0P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE0P_Rotlet_direct_full_ext_mex PUBLIC ZERO_PERIODIC EWALD_FULL EXTERNAL)

# 0P real
matlab_add_mex(
  NAME SE0P_Rotlet_direct_real_mex
  SRC SE0P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE0P_Rotlet_direct_real_mex PUBLIC ZERO_PERIODIC EWALD_REAL)

# 0P real_cut
matlab_add_mex(
  NAME SE0P_Rotlet_direct_real_cut_mex
  SRC SE0P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE0P_Rotlet_direct_real_cut_mex PUBLIC ZERO_PERIODIC EWALD_REAL CUTOFF)

# 0P real_cut_ext
matlab_add_mex(
  NAME SE0P_Rotlet_direct_real_cut_ext_mex
  SRC SE0P_Rotlet_direct.c SE_Rotlet_direct_mex.c
)
target_compile_definitions(SE0P_Rotlet_direct_real_cut_ext_mex PUBLIC ZERO_PERIODIC EWALD_REAL CUTOFF EXTERNAL)
