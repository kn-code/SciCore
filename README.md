_SciCore_ is a C++20 scientific computing library. The main focus lies on computations with matrix-valued data.
For example, given a function \f$ f: \mathbb{R} \rightarrow \mathbb{C}^{n\times m} \f$ one can interpolate \f$f\f$ using Chebyshev polynomials, integrate it with various methods or use it in the solution of
(integro)-differential equations.

## Dependencies

This library uses some C++-20 features and therefore needs a not too old compiler. The following compilers were used successfully in the build process:

* gcc 11.4.0

Some functionalities require _blas_ and _lapacke_ packages, which are on Ubuntu-type systems most easily installed with

```console
apt install libopenblas-dev liblapacke-dev
```

The following dependencies will be automatically downloaded in the build process:

* [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) for linear algebra and representation of matrix/vector types
* [cereal](https://uscilab.github.io/cereal/) for serialization
* [Taskflow](https://github.com/taskflow/taskflow) for multithreading capabilities
* [GoogleTest](https://github.com/google/googletest) for tests

## Building

The minimal steps to build and install the library are

```console
$ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
$ cmake --build build --target install
```

A custom installation path can be chosen when configuring the project via

```console
$ cmake -S . -B build  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install
```

Tests can be compiled and run with

```console
$ cmake --build build --target build_tests
$ cmake --build build --target test
```

## License

This library is licensed under the [Mozilla Public License (MPL) version 2.0](https://www.mozilla.org/en-US/MPL/2.0/FAQ/).