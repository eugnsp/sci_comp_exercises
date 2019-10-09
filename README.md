# Exercises in scientific computing

These are various exercises in scientific computing.

## Contents

### Sequential (`/seq`)

| Directory / File						| Description 															|
|:--------------------------------------|:----------------------------------------------------------------------|
| `kahan_summation.cpp`					| Kahan summation algorithm												|
| `shift_and_subtract_division.cpp`		| Shift-and-subtract division division algorithm						|
| `/lin/laplace/iterative`				| Laplace BVP, Jacobi, Gauss-Seidel, and SOR iterative methods			|
| `/lin/given_qr_factorization`			| QR factorization, Givens rotation										|
| `/eigen/jacobi.cpp`					| Matrix diagonalization, Jacobi eigenvalue algorithm					|
| `/random/ising_model`					| 2D Ising model simulation with Metropolis Monte-Carlo algorithm		|
| `/random/uniform_distr_on_sphere`		| Uniform distribution on sphere										|
| `/wavelet/multiresolution_analysis`	| Wavelets, discrete wavelet transform									|

### MPI (`/mpi`)

| Directory / File						| Description 															|
|:--------------------------------------|:----------------------------------------------------------------------|
| `mat_mat_multiplication.cpp`			| Matrix-matrix multiplication											|

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

C++17 compiler is required. Tested with GCC 8.3.0 and Clang 7.0.

## External dependencies

* [Intel MKL](https://software.intel.com/en-us/mkl)
* [`esl` library](https://github.com/eugnsp/esl)
* MPI library, e.g. [Intel MPI](https://software.intel.com/en-us/mpi-library)

## License

This code is distributed under GNU General Public License v3.0.
