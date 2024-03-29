# Exercises in scientific computing

These are various exercises in scientific computing.

## Contents

### Sequential (`seq`)

* `basics`
	* `ieee754`
		* [`doubles_density`](src/seq/basics/ieee754/doubles_density/): Relative density of double floating-point numbers
		* [`info`](src/seq/basics/ieee754/info/): Extracting information from IEEE 754 floating-point numbers
		* [`ulps`](src/seq/basics/ieee754/ulps/): Number of full floating-point intervals between two numbers
	* `summation`
		* [`exp_taylor_series`](src/seq/basics/summation/exp_taylor_series/): Exponent calculation via Taylor series summation
		* [`kahan_summation`](src/seq/basics/summation/kahan_summation/): Kahan summation algorithm
	* [`shift_and_subtract_division`](src/seq/basics/shift_and_subtract_division/): Shift-and-subtract division division algorithm
* `lin`
	* `laplace_2d_fdm`
		* [`relaxation`](src/seq/lin/laplace_2d_fdm/relaxation/): Laplace BVP solution using finite-difference discretization and the Jacobi, the Gauss-Seidel, and the SOR iterative methods
	* [`poisson_2d_fem`](src/seq/lin/poisson_2d_fem/): Poisson BVP solution using mixed finite-element P0+P1 method
	* [`given_qr_factorization`](src/seq/lin/given_qr_factorization/): QR factorization using Givens rotation
* `eigen`
	* [`jacobi`](src/seq/eigen/jacobi/): Matrix diagonalization using the Jacobi eigenvalue algorithm
* `random`
	* [`ising_model`](src/seq/random/ising_model/): 2D Ising model simulation with Metropolis Monte-Carlo algorithm
	* [`uniform_distr_on_sphere`](src/seq/random/uniform_distr_on_sphere/): Uniform distribution on sphere
* `wavelet`
	* [`multiresolution_analysis`](src/seq/wavelet/multiresolution_analysis/): Wavelets, discrete wavelet transform

### MPI (`mpi`)

* `mat_mat_multiplication.cpp`: Matrix-matrix multiplication

<!--| 20	| 4.4		| LU factorization					| LU factorization without pivoting, MPI						|-->

## Results

Generated images can be found in the corresponding directories inside the `figs` directory. Video animations can be found on the following [YouTube channel](https://www.youtube.com/channel/UCvaVjVoG0KRS9TGZaBZxfnQ).

## How to build

Set `MKLROOT` environment variable to point to the MKL installation directory,
and be sure that your CMake version is >= 3.13. Then:

```sh
git clone --recursive https://github.com/eugnsp/sci_comp_exercises.git
cd sci_comp_exercises
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE .. && make
```

C++17/C++20 compiler is required. Tested with GCC 8.3.0.

## External dependencies

* [Intel MKL](https://software.intel.com/en-us/mkl)
* [`eslib` library](https://github.com/eugnsp/eslib)
* MPI library, e.g. [Intel MPI](https://software.intel.com/en-us/mpi-library)

## License

This code is distributed under GNU General Public License v3.0.
