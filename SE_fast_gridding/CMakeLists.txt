## FGG MEX libraries
# 3D
add_library(
  fgg
  src/SE_fgg.c src/SE_fgg.h src/SE_fgg_MEX_params.c
)
target_compile_definitions(fgg PUBLIC THREE_PERIODIC)
target_link_libraries(fgg gsl gslcblas)

add_library(
  fgg_thrd
  src/SE_fgg.c src/SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_thrd PUBLIC THREE_PERIODIC FGG_THRD)
target_link_libraries(fgg_thrd gsl gslcblas)

add_library(
  fgg_force
  src/SE_fgg.c src/SE_fgg.h src/SE_fgg_force.c src/SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_force PUBLIC THREE_PERIODIC FORCE)
target_link_libraries(fgg_force gsl gslcblas)

## FGG MEX libraries

# 2D
add_library(
  fgg_2p
  src/SE_fgg.c src/SE_fgg.h src/SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_2p PUBLIC TWO_PERIODIC)
target_link_libraries(fgg_2p gsl gslcblas)

add_library(
  fgg_thrd_2p
  src/SE_fgg.c src/SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_thrd_2p PUBLIC TWO_PERIODIC FGG_THRD)
target_link_libraries(fgg_thrd_2p gsl gslcblas)

add_library(
  fgg_force_2p
  src/SE_fgg.c src/SE_fgg.h src/SE_fgg_force.c src/SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_force_2p PUBLIC TWO_PERIODIC FORCE)
target_link_libraries(fgg_force_2p gsl gslcblas)

## FGG MEX libraries

# 1D
add_library(
  fgg_1p
  src/SE_fgg.c src/SE_fgg.h src/SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_1p PUBLIC ONE_PERIODIC)
target_link_libraries(fgg_1p gsl gslcblas)

add_library(
  fgg_thrd_1p
  src/SE_fgg.c src/SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_thrd_1p PUBLIC ONE_PERIODIC FGG_THRD)
target_link_libraries(fgg_thrd_1p gsl gslcblas)

add_library(
  fgg_force_1p
  src/SE_fgg.c src/SE_fgg.h src/SE_fgg_MEX_params.c
)
target_compile_definitions(fgg_force_1p PUBLIC ONE_PERIODIC FORCE)
target_link_libraries(fgg_force_1p gsl gslcblas)


## MEX functions 3P

matlab_add_mex(
  NAME SE_fg_grid_mex
  SRC mex/SE_fg_grid_mex.c
)
target_link_libraries(SE_fg_grid_mex fgg)

matlab_add_mex(
  NAME SE_fg_int_mex
  SRC mex/SE_fg_int_mex.c
)
target_link_libraries(SE_fg_int_mex fgg)

matlab_add_mex(
  NAME SE_fgg_expand_all_mex
  SRC mex/SE_fgg_expand_all_mex.c  src/SE_fgg.c src/SE_fgg_MEX_params.c
)
target_link_libraries(SE_fgg_expand_all_mex fgg gsl gslcblas)

matlab_add_mex(
  NAME SE_fgg_base_gaussian_mex
  SRC mex/SE_fgg_base_gaussian_mex.c
)
target_link_libraries(SE_fgg_base_gaussian_mex fgg)

matlab_add_mex(
  NAME SE_fg_grid_split_mex
  SRC mex/SE_fg_grid_split_mex.c
)
target_link_libraries(SE_fg_grid_split_mex fgg)

matlab_add_mex(
  NAME SE_fg_int_split_mex
  SRC mex/SE_fg_int_split_mex.c
)
target_link_libraries(SE_fg_int_split_mex fgg)

matlab_add_mex(
  NAME SE_fg_grid_thrd_mex
  SRC mex/SE_fg_grid_thrd_mex.c
)
target_link_libraries(SE_fg_grid_thrd_mex fgg_thrd)

matlab_add_mex(
  NAME SE_fg_grid_split_thrd_mex
  SRC mex/SE_fg_grid_split_thrd_mex.c
)
target_link_libraries(SE_fg_grid_split_thrd_mex fgg_thrd)

matlab_add_mex(
  NAME SE_fg_extend_fcn_mex
  SRC mex/SE_fg_extend_fcn_mex.c
)
target_link_libraries(SE_fg_extend_fcn_mex fgg)

matlab_add_mex(
  NAME SE_fg_int_extended_split_mex
  SRC mex/SE_fg_int_extended_split_mex.c
)
target_link_libraries(SE_fg_int_extended_split_mex fgg)

## MEX functions 2P
matlab_add_mex(
  NAME SE_fg_grid_mex_2p
  SRC mex/SE_fg_grid_mex.c
)
target_link_libraries(SE_fg_grid_mex_2p fgg_2p)

matlab_add_mex(
  NAME SE_fg_int_mex_2p
  SRC mex/SE_fg_int_mex.c
)
target_link_libraries(SE_fg_int_mex_2p fgg_2p)

matlab_add_mex(
  NAME SE_fgg_expand_all_mex_2p
  SRC mex/SE_fgg_expand_all_mex.c  src/SE_fgg.c src/SE_fgg_MEX_params.c
)
target_link_libraries(SE_fgg_expand_all_mex_2p fgg_2p gsl gslcblas)

matlab_add_mex(
  NAME SE_fgg_base_gaussian_mex_2p
  SRC mex/SE_fgg_base_gaussian_mex.c
)
target_link_libraries(SE_fgg_base_gaussian_mex_2p fgg_2p)

matlab_add_mex(
  NAME SE_fg_grid_split_mex_2p
  SRC mex/SE_fg_grid_split_mex.c
)
target_link_libraries(SE_fg_grid_split_mex_2p fgg_2p)

matlab_add_mex(
  NAME SE_fg_int_split_mex_2p
  SRC mex/SE_fg_int_split_mex.c
)
target_link_libraries(SE_fg_int_split_mex_2p fgg_2p)

matlab_add_mex(
  NAME SE_fg_grid_split_thrd_mex_2p
  SRC mex/SE_fg_grid_split_thrd_mex.c
)
target_link_libraries(SE_fg_grid_split_thrd_mex_2p fgg_thrd_2p)

matlab_add_mex(
  NAME SE_fg_extend_fcn_mex_2p
  SRC mex/SE_fg_extend_fcn_mex.c
)
target_link_libraries(SE_fg_extend_fcn_mex_2p fgg_2p)

matlab_add_mex(
  NAME SE_fg_int_extended_split_mex_2p
  SRC mex/SE_fg_int_extended_split_mex.c
)
target_link_libraries(SE_fg_int_extended_split_mex_2p fgg_2p)

## MEX functions 1P
matlab_add_mex(
  NAME SE_fg_grid_mex_1p
  SRC mex/SE_fg_grid_mex.c
)
target_link_libraries(SE_fg_grid_mex_1p fgg_1p)

matlab_add_mex(
  NAME SE_fg_int_mex_1p
  SRC mex/SE_fg_int_mex.c
)
target_link_libraries(SE_fg_int_mex_1p fgg_1p)

matlab_add_mex(
  NAME SE_fgg_expand_all_mex_1p
  SRC mex/SE_fgg_expand_all_mex.c  src/SE_fgg.c src/SE_fgg_MEX_params.c
)
target_link_libraries(SE_fgg_expand_all_mex_1p fgg_1p gsl gslcblas)

matlab_add_mex(
  NAME SE_fgg_base_gaussian_mex_1p
  SRC mex/SE_fgg_base_gaussian_mex.c
)
target_link_libraries(SE_fgg_base_gaussian_mex_1p fgg_1p)

matlab_add_mex(
  NAME SE_fg_grid_split_mex_1p
  SRC mex/SE_fg_grid_split_mex.c
)
target_link_libraries(SE_fg_grid_split_mex_1p fgg_1p)

matlab_add_mex(
  NAME SE_fg_int_split_mex_1p
  SRC mex/SE_fg_int_split_mex.c
)
target_link_libraries(SE_fg_int_split_mex_1p fgg_1p)

matlab_add_mex(
  NAME SE_fg_grid_split_thrd_mex_1p
  SRC mex/SE_fg_grid_split_thrd_mex.c
)
target_link_libraries(SE_fg_grid_split_thrd_mex_1p fgg_thrd_1p)

matlab_add_mex(
  NAME SE_fg_extend_fcn_mex_1p
  SRC mex/SE_fg_extend_fcn_mex.c
)
target_link_libraries(SE_fg_extend_fcn_mex_1p fgg_1p)

matlab_add_mex(
  NAME SE_fg_int_extended_split_mex_1p
  SRC mex/SE_fg_int_extended_split_mex.c
)
target_link_libraries(SE_fg_int_extended_split_mex_1p fgg_1p)

# MEX functions 3P FORCE
matlab_add_mex(
  NAME SE_fgg_expand_all_force_mex
  SRC mex/SE_fgg_expand_all_force_mex.c  src/SE_fgg.c src/SE_fgg_MEX_params.c
)
target_link_libraries(SE_fgg_expand_all_force_mex fgg_force gsl gslcblas)

matlab_add_mex(
  NAME SE_fg_int_force_mex
  SRC mex/SE_fg_int_force_mex.c
)
target_link_libraries(SE_fg_int_force_mex fgg_force)

matlab_add_mex(
  NAME SE_fg_int_split_force_mex
  SRC mex/SE_fg_int_split_force_mex.c
)
target_link_libraries(SE_fg_int_split_force_mex fgg_force)

# MEX functions 2P FORCE
matlab_add_mex(
  NAME SE_fgg_expand_all_force_mex_2p
  SRC mex/SE_fgg_expand_all_force_mex.c src/SE_fgg.c src/SE_fgg_MEX_params.c
)
target_link_libraries(SE_fgg_expand_all_force_mex_2p fgg_force_2p gsl gslcblas)

matlab_add_mex(
  NAME SE_fg_int_force_mex_2p
  SRC mex/SE_fg_int_force_mex.c
)
target_link_libraries(SE_fg_int_force_mex_2p fgg_force_2p)

matlab_add_mex(
  NAME SE_fg_int_split_force_mex_2p
  SRC mex/SE_fg_int_split_force_mex.c
)
target_link_libraries(SE_fg_int_split_force_mex_2p fgg_force_2p)

# MEX functions 1P FORCE
matlab_add_mex(
  NAME SE_fgg_expand_all_force_mex_1p
  SRC mex/SE_fgg_expand_all_force_mex.c  src/SE_fgg.c src/SE_fgg_MEX_params.c
)
target_link_libraries(SE_fgg_expand_all_force_mex_1p fgg_force_1p gsl gslcblas)

matlab_add_mex(
  NAME SE_fg_int_force_mex_1p
  SRC mex/SE_fg_int_force_mex.c
)
target_link_libraries(SE_fg_int_force_mex_1p fgg_force_1p)

matlab_add_mex(
  NAME SE_fg_int_split_force_mex_1p
  SRC mex/SE_fg_int_split_force_mex.c
)
target_link_libraries(SE_fg_int_split_force_mex_1p fgg_force_1p)


############## KAISER MEX FUNCTIONS
## FGG MEX libraries
# 3D
add_library(
  fkg
  src/SE_fgg.c src/SE_fkg.h src/SE_fgg_MEX_params.c
)
target_compile_definitions(fkg PUBLIC THREE_PERIODIC KAISER)
target_link_libraries(fkg gsl gslcblas)

add_library(
  fkg_thrd
  src/SE_fgg.c src/SE_fgg_MEX_params.c
)
target_compile_definitions(fkg_thrd PUBLIC THREE_PERIODIC KAISER FGG_THRD)
target_link_libraries(fkg_thrd gsl gslcblas)

add_library(
  fkg_force
  src/SE_fgg.c src/SE_fkg.h src/SE_fkg_force.c src/SE_fgg_MEX_params.c
)
target_compile_definitions(fkg_force PUBLIC THREE_PERIODIC KAISER FORCE)
target_link_libraries(fkg_force gsl gslcblas)

## FGG MEX libraries

# 2D
add_library(
  fkg_2p
  src/SE_fgg.c src/SE_fkg.h src/SE_fgg_MEX_params.c
)
target_compile_definitions(fkg_2p PUBLIC TWO_PERIODIC KAISER)
target_link_libraries(fkg_2p gsl gslcblas)

add_library(
  fkg_thrd_2p
  src/SE_fgg.c src/SE_fgg_MEX_params.c
)
target_compile_definitions(fkg_thrd_2p PUBLIC TWO_PERIODIC KAISER FGG_THRD)
target_link_libraries(fkg_thrd_2p gsl gslcblas)

add_library(
  fkg_force_2p
  src/SE_fgg.c src/SE_fkg.h src/SE_fkg_force.c src/SE_fgg_MEX_params.c
)
target_compile_definitions(fkg_force_2p PUBLIC TWO_PERIODIC KAISER FORCE)
target_link_libraries(fkg_force_2p gsl gslcblas)

## FGG MEX libraries

# 1D
add_library(
  fkg_1p
  src/SE_fgg.c src/SE_fkg.h src/SE_fgg_MEX_params.c
)
target_compile_definitions(fkg_1p PUBLIC ONE_PERIODIC KAISER)
target_link_libraries(fkg_1p gsl gslcblas)

add_library(
  fkg_thrd_1p
  src/SE_fgg.c src/SE_fgg_MEX_params.c
)
target_compile_definitions(fkg_thrd_1p PUBLIC ONE_PERIODIC KAISER FGG_THRD)
target_link_libraries(fkg_thrd_1p gsl gslcblas)

add_library(
  fkg_force_1p
  src/SE_fgg.c src/SE_fkg.h src/SE_fkg_force.c src/SE_fgg_MEX_params.c
)
target_compile_definitions(fkg_force_1p PUBLIC ONE_PERIODIC KAISER FORCE)
target_link_libraries(fkg_force_1p gsl gslcblas)

#  KAISER MEX functions 3P
matlab_add_mex(
  NAME SE_fg_grid_kaiser_mex
  SRC mex/SE_fg_grid_kaiser_mex.c
)
target_link_libraries(SE_fg_grid_kaiser_mex fkg_thrd)

matlab_add_mex(
  NAME SE_fg_int_kaiser_mex
  SRC mex/SE_fg_int_kaiser_mex.c
)
target_link_libraries(SE_fg_int_kaiser_mex fkg)
matlab_add_mex(
  NAME SE_fg_int_kaiser_force_mex
  SRC mex/SE_fg_int_kaiser_force_mex.c
)
target_link_libraries(SE_fg_int_kaiser_force_mex fkg_force)

matlab_add_mex(
  NAME SE_fg_grid_split_kaiser_mex
  SRC mex/SE_fg_grid_split_kaiser_mex.c
)
target_link_libraries(SE_fg_grid_split_kaiser_mex fkg_thrd)

matlab_add_mex(
  NAME SE_fg_int_split_kaiser_mex
  SRC mex/SE_fg_int_split_kaiser_mex.c
)
target_link_libraries(SE_fg_int_split_kaiser_mex fkg)
matlab_add_mex(
  NAME SE_fg_int_split_kaiser_force_mex
  SRC mex/SE_fg_int_split_kaiser_force_mex.c
)
target_link_libraries(SE_fg_int_split_kaiser_force_mex fkg_force)

matlab_add_mex(
  NAME SE_fg_extend_fcn_kaiser_mex
  SRC mex/SE_fg_extend_fcn_kaiser_mex.c
)
target_link_libraries(SE_fg_extend_fcn_kaiser_mex fkg)

matlab_add_mex(
  NAME SE_fkg_expand_all_mex
  SRC mex/SE_fkg_expand_all_mex.c
)
target_link_libraries(SE_fkg_expand_all_mex fkg)
matlab_add_mex(
  NAME SE_fkg_expand_all_force_mex
  SRC mex/SE_fkg_expand_all_force_mex.c
)
target_link_libraries(SE_fkg_expand_all_force_mex fkg_force)

#  KAISER MEX functions 2P
matlab_add_mex(
  NAME SE_fg_grid_kaiser_mex_2p
  SRC mex/SE_fg_grid_kaiser_mex.c
)
target_link_libraries(SE_fg_grid_kaiser_mex_2p fkg_thrd_2p)

matlab_add_mex(
  NAME SE_fg_int_kaiser_mex_2p
  SRC mex/SE_fg_int_kaiser_mex.c
)
target_link_libraries(SE_fg_int_kaiser_mex_2p fkg_2p)
matlab_add_mex(
  NAME SE_fg_int_kaiser_force_mex_2p
  SRC mex/SE_fg_int_kaiser_force_mex.c
)
target_link_libraries(SE_fg_int_kaiser_force_mex_2p fkg_force_2p)

matlab_add_mex(
  NAME SE_fg_grid_split_kaiser_mex_2p
  SRC mex/SE_fg_grid_split_kaiser_mex.c
)
target_link_libraries(SE_fg_grid_split_kaiser_mex_2p fkg_thrd_2p)

matlab_add_mex(
  NAME SE_fg_int_split_kaiser_mex_2p
  SRC mex/SE_fg_int_split_kaiser_mex.c
)
target_link_libraries(SE_fg_int_split_kaiser_mex_2p fkg_2p)
matlab_add_mex(
  NAME SE_fg_int_split_kaiser_force_mex_2p
  SRC mex/SE_fg_int_split_kaiser_force_mex.c
)
target_link_libraries(SE_fg_int_split_kaiser_force_mex_2p fkg_force_2p)

matlab_add_mex(
  NAME SE_fkg_expand_all_mex_2p
  SRC mex/SE_fkg_expand_all_mex.c
)
target_link_libraries(SE_fkg_expand_all_mex_2p fkg_2p)
matlab_add_mex(
  NAME SE_fkg_expand_all_force_mex_2p
  SRC mex/SE_fkg_expand_all_force_mex.c
)
target_link_libraries(SE_fkg_expand_all_force_mex_2p fkg_force_2p)

#  KAISER MEX functions 1P
matlab_add_mex(
  NAME SE_fg_grid_kaiser_mex_1p
  SRC mex/SE_fg_grid_kaiser_mex.c
)
target_link_libraries(SE_fg_grid_kaiser_mex_1p fkg_thrd_1p)

matlab_add_mex(
  NAME SE_fg_int_kaiser_mex_1p
  SRC mex/SE_fg_int_kaiser_mex.c
)
target_link_libraries(SE_fg_int_kaiser_mex_1p fkg_1p)
matlab_add_mex(
  NAME SE_fg_int_kaiser_force_mex_1p
  SRC mex/SE_fg_int_kaiser_force_mex.c
)
target_link_libraries(SE_fg_int_kaiser_force_mex_1p fkg_force_1p)

matlab_add_mex(
  NAME SE_fg_grid_split_kaiser_mex_1p
  SRC mex/SE_fg_grid_split_kaiser_mex.c
)
target_link_libraries(SE_fg_grid_split_kaiser_mex_1p fkg_thrd_1p)

matlab_add_mex(
  NAME SE_fg_int_split_kaiser_mex_1p
  SRC mex/SE_fg_int_split_kaiser_mex.c
)
target_link_libraries(SE_fg_int_split_kaiser_mex_1p fkg_1p)
matlab_add_mex(
  NAME SE_fg_int_split_kaiser_force_mex_1p
  SRC mex/SE_fg_int_split_kaiser_force_mex.c
)
target_link_libraries(SE_fg_int_split_kaiser_force_mex_1p fkg_force_1p)

matlab_add_mex(
  NAME SE_fkg_expand_all_mex_1p
  SRC mex/SE_fkg_expand_all_mex.c
)
target_link_libraries(SE_fkg_expand_all_mex_1p fkg_1p)
matlab_add_mex(
  NAME SE_fkg_expand_all_force_mex_1p
  SRC mex/SE_fkg_expand_all_force_mex.c
)
target_link_libraries(SE_fkg_expand_all_force_mex_1p fkg_force_1p)
