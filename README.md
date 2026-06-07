# hompack

[![Test](https://github.com/HugoMVale/hompack/actions/workflows/test_gcc.yml/badge.svg)](https://github.com/HugoMVale/hompack/actions)
[![codecov](https://codecov.io/gh/HugoMVale/hompack/graph/badge.svg?token=iaVkdNeKJj)](https://codecov.io/gh/HugoMVale/hompack)
[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)


## Description

`hompack` is a package for finding zeros or fixed points of nonlinear systems using globally convergent probability-one homotopy algorithms.

## History

The first version of the library, named HOMPACK, was originally released in 1987 [^1], and the last "official" update, named HOMPACK90, dates from 1997 [^2].

`hompack` is (will be) a modernization of the HOMPACK90 code [^3], intended to make the library easier to use and maintain. The main changes include:

* [ ] Conversion from fixed-form (`.f`) to free-form (`.f90`).
* [ ] Conversion from upper case to lower case.
* [x] Modularization.
* [ ] Removal of `DATA` statements, labeled do loops, and (most) `goto`s.
* [ ] Addition of `intent(in/out)` to all procedures.
* [x] Addition of explicit interfaces to BLAS/LAPACK routines.
* [ ] Implementation of a C API.
* [x] Automatic code documentation with FORD.
* [ ] Python bindings, available in the companion repo [hompack-python](https://github.com/HugoMVale/hompack-python).

|    Version    | Year |   Standard   |
|:-------------:|:----:|:------------:|
|  hompack      | 202x | Fortran 2018 |
|  HOMPACK90    | 1997 |  FORTRAN 90  |
|  HOMPACK      | 1987 |  FORTRAN 77  |


## Build instructions

### Dependencies

`hompack` depends on a number of functions from BLAS and LAPACK.

The build configuration files provided with the code (see further below) assume a suitable BLAS/LAPACK library is locally installed. Alternatively, the subset of required functions is available in [./external](./external).


### With fpm

The easiest way to build/test the code and run the examples is by means of [`fpm`](https://fpm.fortran-lang.org/).

To build the library, do:

```sh
fpm build --profile release
```

To run the tests, do:

```sh
fpm test --profile release
```

To run the provided examples, do:

```sh
fpm run --example "example_name" --profile release
```

### With meson

First, setup the build:

```sh
meson setup builddir -Dbuild_tests=true
```

To build the libraries, do:

```sh
meson compile -C builddir
```

To run the tests, do:

```sh
meson test -C builddir
```

## Licence

* The original `hompack` code is public domain.
* Modifications introduced in this project are covered under the MIT license.


### References

[^1]: Layne T. Watson, Stephen C. Billups, and Alexander P. Morgan. 1987. Algorithm 652: HOMPACK: a suite of codes for globally convergent homotopy algorithms. ACM Trans. Math. Softw. 13, 3 (Sept. 1987), 281–310. https://doi.org/10.1145/29380.214343

[^2]: Layne T. Watson, Maria Sosonkina, Robert C. Melville, Alexander P. Morgan, and Homer F. Walker. 1997. Algorithm 777: HOMPACK90: a suite of Fortran 90 codes for globally convergent homotopy algorithms. ACM Trans. Math. Softw. 23, 4 (Dec. 1997), 514–549. https://doi.org/10.1145/279232.279235

[^3]: Original source code from [Netlib](https://www.netlib.org/hompack/).
