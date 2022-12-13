## Fast real-space code using cell list
if(MKL_FOUND)
  # Build SE3P_Laplace realspace using cell list and MKL
  matlab_add_mex(
    NAME SE3P_Laplace_real_cut_cell_mkl_mex
    SRC SE3P_Laplace_real_cut_cell.c cell_list.c
    )
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE3P_Laplace_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE3P_Laplace_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_compile_definitions(SE3P_Laplace_real_cut_cell_mkl_mex PUBLIC INTEL_MKL)

  # Build SE3P_Rotlet realspace using cell list and MKL
  matlab_add_mex(
    NAME SE3P_Rotlet_real_cut_cell_mkl_mex
    SRC SE3P_Rotlet_real_cut_cell.c cell_list.c
    )
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE3P_Rotlet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE3P_Rotlet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_compile_definitions(SE3P_Rotlet_real_cut_cell_mkl_mex PUBLIC INTEL_MKL)

  # Build SE3P_Stokeslet realspace using cell list and MKL
  matlab_add_mex(
    NAME SE3P_Stokeslet_real_cut_cell_mkl_mex
    SRC SE3P_Stokeslet_real_cut_cell.c cell_list.c
    )
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE3P_Stokeslet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE3P_Stokeslet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_compile_definitions(SE3P_Stokeslet_real_cut_cell_mkl_mex PUBLIC INTEL_MKL)

  # Build SE3P_Stresslet realspace using cell list and MKL
  matlab_add_mex(
    NAME SE3P_Stresslet_real_cut_cell_mkl_mex
    SRC SE3P_Stresslet_real_cut_cell.c cell_list.c
    )
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE3P_Stresslet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE3P_Stresslet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_compile_definitions(SE3P_Stresslet_real_cut_cell_mkl_mex PUBLIC INTEL_MKL)
endif()

## Standard, non-MKL, code
# Laplace
matlab_add_mex(
  NAME SE3P_Laplace_real_cut_cell_mex
  SRC SE3P_Laplace_real_cut_cell.c cell_list.c
  )
matlab_add_mex(
  NAME SE3P_Laplace_real_cut_cell_grad_mex
  SRC SE3P_Laplace_real_cut_cell.c cell_list.c
  )
target_compile_definitions(SE3P_Laplace_real_cut_cell_grad_mex PUBLIC GRADIENT)

# Rotlet
matlab_add_mex(
  NAME SE3P_Rotlet_real_cut_cell_mex
  SRC SE3P_Rotlet_real_cut_cell.c cell_list.c
  )

# Stokeslet
matlab_add_mex(
  NAME SE3P_Stokeslet_real_cut_cell_mex
  SRC SE3P_Stokeslet_real_cut_cell.c cell_list.c
  )

# Stresslet
matlab_add_mex(
  NAME SE3P_Stresslet_real_cut_cell_mex
  SRC SE3P_Stresslet_real_cut_cell.c cell_list.c
  )
