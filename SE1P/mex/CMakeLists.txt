## Fast real-space code using cell list
if(MKL_FOUND)
  # Build SE1P_Laplace realspace using cell list and MKL
  matlab_add_mex(
    NAME SE1P_Laplace_real_cut_cell_mkl_mex
    SRC SE1P_Laplace_real_cut_cell.c ../../SE3P/mex/cell_list.c
    )
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE1P_Laplace_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE1P_Laplace_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_compile_definitions(SE1P_Laplace_real_cut_cell_mkl_mex PUBLIC INTEL_MKL)

  # Build SE1P_Rotlet realspace using cell list and MKL
  matlab_add_mex(
    NAME SE1P_Rotlet_real_cut_cell_mkl_mex
    SRC SE1P_Rotlet_real_cut_cell.c ../../SE3P/mex/cell_list.c
    )
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE1P_Rotlet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE1P_Rotlet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_compile_definitions(SE1P_Rotlet_real_cut_cell_mkl_mex PUBLIC INTEL_MKL)

  # Build SE1P_Stokeslet realspace using cell list and MKL
  matlab_add_mex(
    NAME SE1P_Stokeslet_real_cut_cell_mkl_mex
    SRC SE1P_Stokeslet_real_cut_cell.c ../../SE3P/mex/cell_list.c
    )
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE1P_Stokeslet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE1P_Stokeslet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_compile_definitions(SE1P_Stokeslet_real_cut_cell_mkl_mex PUBLIC INTEL_MKL)

  # Build SE1P_Stresslet realspace using cell list and MKL
  matlab_add_mex(
    NAME SE1P_Stresslet_real_cut_cell_mkl_mex
    SRC SE1P_Stresslet_real_cut_cell.c ../../SE3P/mex/cell_list.c
    )
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE1P_Stresslet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE1P_Stresslet_real_cut_cell_mkl_mex ${MKL_STATIC})
  target_compile_definitions(SE1P_Stresslet_real_cut_cell_mkl_mex PUBLIC INTEL_MKL)
endif()

## Standard, non-MKL, code
# Laplace
matlab_add_mex(
  NAME SE1P_Laplace_real_cut_cell_mex
  SRC SE1P_Laplace_real_cut_cell.c ../../SE3P/mex/cell_list.c
  )
matlab_add_mex(
  NAME SE1P_Laplace_real_cut_cell_grad_mex
  SRC SE1P_Laplace_real_cut_cell.c ../../SE3P/mex/cell_list.c
  )
target_compile_definitions(SE1P_Laplace_real_cut_cell_grad_mex PUBLIC GRADIENT)

# Rotlet
matlab_add_mex(
  NAME SE1P_Rotlet_real_cut_cell_mex
  SRC SE1P_Rotlet_real_cut_cell.c ../../SE3P/mex/cell_list.c
  )

# Stokeslet
matlab_add_mex(
  NAME SE1P_Stokeslet_real_cut_cell_mex
  SRC SE1P_Stokeslet_real_cut_cell.c ../../SE3P/mex/cell_list.c
  )

# Stresslet
matlab_add_mex(
  NAME SE1P_Stresslet_real_cut_cell_mex
  SRC SE1P_Stresslet_real_cut_cell.c ../../SE3P/mex/cell_list.c
  )
