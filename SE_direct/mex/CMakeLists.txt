## Libraries used by direct summation code

add_library(
  incomplete_bessel_K
  incomplete_bessel_K.c)

## Subdirectories with more MEX files
add_subdirectory("Laplace")
add_subdirectory("Rotlet")
add_subdirectory("Stokeslet")
add_subdirectory("Stresslet")
