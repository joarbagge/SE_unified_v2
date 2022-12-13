## Fast real-space code using cell list
if(MKL_FOUND)
  # Build SE0P_Laplace realspace using cell list and MKL
  matlab_add_mex(
    NAME SE0P_Laplace_real_cut_cell_mkl_mex
    SRC SE0P_Laplace_real_cut_cell.c ../../SE3P/mex/cell_list.c
    )
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE0P_Laplace_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE0P_Laplace_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_compile_definitions(SE0P_Laplace_real_cut_cell_mkl_mex PUBLIC INTEL_MKL)

  # Build SE0P_Rotlet realspace using cell list and MKL
  matlab_add_mex(
    NAME SE0P_Rotlet_real_cut_cell_mkl_mex
    SRC SE0P_Rotlet_real_cut_cell.c ../../SE3P/mex/cell_list.c
    )
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE0P_Rotlet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE0P_Rotlet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_compile_definitions(SE0P_Rotlet_real_cut_cell_mkl_mex PUBLIC INTEL_MKL)

  # Build SE0P_Stokeslet realspace using cell list and MKL
  matlab_add_mex(
    NAME SE0P_Stokeslet_real_cut_cell_mkl_mex
    SRC SE0P_Stokeslet_real_cut_cell.c ../../SE3P/mex/cell_list.c
    )
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE0P_Stokeslet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE0P_Stokeslet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_compile_definitions(SE0P_Stokeslet_real_cut_cell_mkl_mex PUBLIC INTEL_MKL)

  # Build SE0P_Stresslet realspace using cell list and MKL
  matlab_add_mex(
    NAME SE0P_Stresslet_real_cut_cell_mkl_mex
    SRC SE0P_Stresslet_real_cut_cell.c ../../SE3P/mex/cell_list.c
    )
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE0P_Stresslet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE0P_Stresslet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_compile_definitions(SE0P_Stresslet_real_cut_cell_mkl_mex PUBLIC INTEL_MKL)
endif()

## Standard, non-MKL, code
# Laplace
matlab_add_mex(
  NAME SE0P_Laplace_real_cut_cell_mex
  SRC SE0P_Laplace_real_cut_cell.c ../../SE3P/mex/cell_list.c
  )
matlab_add_mex(
  NAME SE0P_Laplace_real_cut_cell_grad_mex
  SRC SE0P_Laplace_real_cut_cell.c ../../SE3P/mex/cell_list.c
  )
target_compile_definitions(SE0P_Laplace_real_cut_cell_grad_mex PUBLIC GRADIENT)

# Rotlet
matlab_add_mex(
  NAME SE0P_Rotlet_real_cut_cell_mex
  SRC SE0P_Rotlet_real_cut_cell.c ../../SE3P/mex/cell_list.c
  )

# Stokeslet
matlab_add_mex(
  NAME SE0P_Stokeslet_real_cut_cell_mex
  SRC SE0P_Stokeslet_real_cut_cell.c ../../SE3P/mex/cell_list.c
  )

# Stresslet
matlab_add_mex(
  NAME SE0P_Stresslet_real_cut_cell_mex
  SRC SE0P_Stresslet_real_cut_cell.c ../../SE3P/mex/cell_list.c
  )

## Kernel FFT precomputation
# Harmonic (Laplace)
matlab_add_mex(
  NAME SE0P_precompute_harmonic_kernel_fft_mex
  SRC SE0P_precompute_kernel_fft_mex.c
  )
target_compile_definitions(SE0P_precompute_harmonic_kernel_fft_mex PUBLIC HARMONIC_KERNEL)
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  matlab_add_mex(
    NAME SE0P_precompute_harmonic_kernel_fft_quad_mex
    SRC SE0P_precompute_kernel_fft_quad_mex.c
    )
  target_compile_definitions(SE0P_precompute_harmonic_kernel_fft_quad_mex PUBLIC HARMONIC_KERNEL)
endif()
if(MPFR_FOUND)
  matlab_add_mex(
    NAME SE0P_precompute_harmonic_kernel_fft_mpfr_mex
    SRC SE0P_precompute_kernel_fft_mpfr_mex.c
    )
  target_compile_definitions(SE0P_precompute_harmonic_kernel_fft_mpfr_mex PUBLIC HARMONIC_KERNEL)
  target_link_libraries(SE0P_precompute_harmonic_kernel_fft_mpfr_mex ${MPFR_LIBRARIES})
endif()

# Biharmonic (Stokeslet)
matlab_add_mex(
  NAME SE0P_precompute_biharmonic_kernel_fft_mex
  SRC SE0P_precompute_kernel_fft_mex.c
  )
target_compile_definitions(SE0P_precompute_biharmonic_kernel_fft_mex PUBLIC BIHARMONIC_KERNEL)
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  matlab_add_mex(
    NAME SE0P_precompute_biharmonic_kernel_fft_quad_mex
    SRC SE0P_precompute_kernel_fft_quad_mex.c
    )
  target_compile_definitions(SE0P_precompute_biharmonic_kernel_fft_quad_mex PUBLIC BIHARMONIC_KERNEL)
endif()
if(MPFR_FOUND)
  matlab_add_mex(
    NAME SE0P_precompute_biharmonic_kernel_fft_mpfr_mex
    SRC SE0P_precompute_kernel_fft_mpfr_mex.c
    )
  target_compile_definitions(SE0P_precompute_biharmonic_kernel_fft_mpfr_mex PUBLIC BIHARMONIC_KERNEL)
  target_link_libraries(SE0P_precompute_biharmonic_kernel_fft_mpfr_mex ${MPFR_LIBRARIES})
endif()
