# Research papers

The two most recent papers, describing the framework that is used in this package to treat all periodic cases in a unified way, are:

1. **(Electrostatics)** D. S. Shamshirgar, J. Bagge and A.-K.
  Tornberg, Fast Ewald summation for electrostatic potentials
  with arbitrary periodicity, *J Chem Phys* **154**(16), 164109
  (2021); https://doi.org/10.1063/5.0044895
2. **(Stokes flow)** J. Bagge and A.-K. Tornberg,
  Fast Ewald summation for Stokes flow with arbitrary periodicity,
  *J Comput Phys* **493**, 112473 (2023); https://doi.org/10.1016/j.jcp.2023.112473

There is a number of papers describing the development of the
Spectral Ewald method, listed below.

**Electrostatics:**

3. D. Lindbo and A.-K. Tornberg, Spectral accuracy in fast Ewald-based methods
  for particle simulations, *J Comput Phys* **230**(24), 8744–8761
  (2011), https://doi.org/10.1016/j.jcp.2011.08.022
4. D. Lindbo and A.-K. Tornberg, Fast and spectrally accurate Ewald summation
  for 2-periodic electrostatic systems, *J Chem Phys* **136**(16), 164111
  (2012); https://doi.org/10.1063/1.4704177
5. A.-K. Tornberg, The Ewald sums for singly, doubly and triply periodic
  electrostatic systems, *Adv Comput Math* **42**, 227–248
  (2016); https://doi.org/10.1007/s10444-015-9422-3
6. D. S. Shamshirgar and A.-K. Tornberg, The Spectral Ewald method for singly
  periodic domains, *J Comput Phys* **347**, 341–366 (2017);
  https://doi.org/10.1016/j.jcp.2017.07.001

**Stokes flow:**

7. D. Lindbo and A.-K. Tornberg, Spectrally accurate fast summation for
  periodic Stokes potentials, *J Comput Phys* **229**(23), 8994–9010
  (2010); https://doi.org/10.1016/j.jcp.2010.08.026
8. D. Lindbo and A.-K. Tornberg, Fast and spectrally accurate summation of
  2-periodic Stokes potentials, *Preprint* (2011); http://arxiv.org/abs/1111.1815
9. L. af Klinteberg and A.-K. Tornberg, Fast Ewald summation for Stokesian particle
  suspensions, *Int J Numer Meth Fluids* **76**, 669–698 (2014);
  https://doi.org/10.1002/fld.3953
10. L. af Klinteberg, Ewald summation for the rotlet singularity of
  Stokes flow, *Preprint* (2016); https://arxiv.org/abs/1603.07467
11. L. af Klinteberg, D. S. Shamshirgar and A.-K. Tornberg, Fast Ewald summation
  for free-space Stokes potentials, *Res Math Sci* **4**, 1 (2017);
  https://doi.org/10.1186/s40687-016-0092-7

## An overview of the papers

The Spectral Ewald method was first formulated in the triply
periodic case for the stokeslet (7.) and harmonic (3.) kernels.
It was subsequently extended to the doubly periodic case for the
stokeslet (8.) and harmonic (4.) kernels.

The method was extended to the triply periodic stresslet (9.)
and rotlet (10.) kernels, which are used in addition to the
stokeslet in boundary integral formulations of Stokes flow.

Efforts to unify the formulation started in the electrostatic
case (5.) and was applied to the singly periodic harmonic kernel
(6.). Around the same time, new ideas by [Vico et al.][https://doi.org/10.1016/j.jcp.2016.07.028]
were used to extend the Spectral Ewald method to the free-space
case for the Stokes kernels (11.).

The unification of all periodic cases was concluded first for the
harmonic kernel (1.), and then for the Stokes kernels (2.). These
papers also introduced a new window function based on the
Kaiser–Bessel window, which replaces the Gaussian window function
previously used in the method.
