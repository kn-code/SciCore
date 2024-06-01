//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

///
/// \file   Definitions.h
///
/// \brief  Basic type definitions and concepts used throughout.
///

#ifndef SCICORE_DEFINITIONS_H
#define SCICORE_DEFINITIONS_H

#include <cmath>
#include <complex>
#include <numbers>
#include <type_traits>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "SciCore_export.h"

namespace SciCore
{

//
// Basic types
//

/// Real number type.
using Real = double;

/// Complex number type.
using Complex = std::complex<Real>;

/// Dense real vector with variable size.
using RealVector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

/// Dense integer vector with variable size.
using IntVector = Eigen::Matrix<int, Eigen::Dynamic, 1>;

/// Dense real rowvector with variable size.
using RealRowVector = Eigen::Matrix<Real, 1, Eigen::Dynamic>;

/// Dense complex vector with variable size.
using Vector = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;

/// Dense complex rowvector with variable size.
using RowVector = Eigen::Matrix<Complex, 1, Eigen::Dynamic>;

/// Dense real matrix with variable number of rows and columns.
using RealMatrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

/// Dense complex matrix with variable number of rows and columns.
using Matrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;

/// Sparse real matrix with variable number of rows and columns.
using SparseRealMatrix = Eigen::SparseMatrix<Real>;

/// Sparse complex matrix with variable number of rows and columns.
using SparseMatrix = Eigen::SparseMatrix<Complex>;

/// Dense real vector with fixed size _n_.
template <int n>
using StaticRealVector = Eigen::Matrix<Real, n, 1>;

/// Dense integer vector with fixed size _n_.
template <int n>
using StaticIntVector = Eigen::Matrix<int, n, 1>;

/// Dense real rowvector with fixed size _n_.
template <int n>
using StaticRealRowVector = Eigen::Matrix<Real, 1, n>;

/// Dense complex vector with fixed size _n_.
template <int n>
using StaticVector = Eigen::Matrix<Complex, n, 1>;

/// Dense complex rowvector with fixed size _n_.
template <int n>
using StaticRowVector = Eigen::Matrix<Complex, 1, n>;

/// Dense real matrix with fixed number of rows _m_ and columns _n_.
template <int m, int n>
using StaticRealMatrix = Eigen::Matrix<Real, m, n>;

/// Dense complex matrix with fixed number of rows _m_ and columns _n_.
template <int m, int n>
using StaticMatrix = Eigen::Matrix<Complex, m, n>;

//
// Basic traits
//

// Type trait determining whether _T_ is a scalar type.
template <typename T>
struct SCICORE_EXPORT IsScalar
{
    constexpr static bool value = std::is_same_v<std::decay_t<T>, Real> || std::is_same_v<std::decay_t<T>, Complex>;
};

// Type trait determining whether _T_ is a dense matrix type. Use as IsDenseMatrix<int>::value
template <typename T>
struct SCICORE_EXPORT IsDenseMatrix : std::is_base_of<Eigen::DenseBase<std::decay_t<T>>, std::decay_t<T>>
{
};

// Use like IsSparseMatrix<int>::value
template <typename T>
struct SCICORE_EXPORT IsSparseMatrix : std::is_base_of<Eigen::SparseMatrixBase<std::decay_t<T>>, std::decay_t<T>>
{
};

template <typename T>
struct SCICORE_EXPORT IsMatrix
{
    constexpr static bool value = IsDenseMatrix<std::decay_t<T>>::value || IsSparseMatrix<std::decay_t<T>>::value;
};

template <typename T>
struct SCICORE_EXPORT IsMatrixOrScalar
{
    constexpr static bool value = IsMatrix<std::decay_t<T>>::value || std::is_same_v<std::decay_t<T>, Real> ||
                                  std::is_same_v<std::decay_t<T>, Complex>;
};

template <typename T>
struct SCICORE_EXPORT IsRealMatrix
{
    constexpr static bool value =
        IsMatrix<T>::value && (std::is_same_v<typename T::Scalar, float> || std::is_same_v<typename T::Scalar, double>);
};

template <typename T>
struct SCICORE_EXPORT IsComplexMatrix
{
    constexpr static bool value = IsMatrix<T>::value && std::is_same_v<typename T::Scalar, Complex>;
};

template <typename T>
struct SCICORE_EXPORT IsStdVector : std::false_type
{
};

template <typename T, typename A>
struct SCICORE_EXPORT IsStdVector<std::vector<T, A>> : std::true_type
{
};

template <typename T>
inline constexpr bool IsStdVector_v = IsStdVector<T>::value;

template <typename T>
struct SCICORE_EXPORT GetScalarType
{
    using type = typename T::Scalar;
};

template <>
struct SCICORE_EXPORT GetScalarType<float>
{
    using type = float;
};

template <>
struct SCICORE_EXPORT GetScalarType<double>
{
    using type = double;
};

template <>
struct SCICORE_EXPORT GetScalarType<std::complex<float>>
{
    using type = std::complex<float>;
};

template <>
struct SCICORE_EXPORT GetScalarType<std::complex<double>>
{
    using type = std::complex<double>;
};

//
// Basic concepts
//

template <typename T>
concept ScalarType = IsScalar<T>::value;

template <typename T>
concept MatrixOrScalarType = IsMatrixOrScalar<T>::value;

template <typename T>
concept MatrixType = IsMatrix<T>::value;

template <typename T>
concept DenseMatrixType = IsDenseMatrix<T>::value;

template <typename T>
concept SparseMatrixType = IsSparseMatrix<T>::value;

} // namespace SciCore

#endif // SCICORE_DEFINITIONS_H
