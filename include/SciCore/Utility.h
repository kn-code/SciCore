//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

///
/// \file   Utility.h
///
/// \brief  Small utility functions used throughout.
///

#ifndef SCICORE_UTILITY_H
#define SCICORE_UTILITY_H

#include <cmath>
#include <type_traits>
#include <vector>

#include "Definitions.h"
#include "SciCore_export.h"

SCICORE_EXPORT inline std::complex<SciCore::Real> operator+(SciCore::Complex z, int n)
{
    return (z + static_cast<SciCore::Real>(n));
}

SCICORE_EXPORT inline SciCore::Complex operator-(SciCore::Complex z, int n)
{
    return (z - static_cast<SciCore::Real>(n));
}

SCICORE_EXPORT inline SciCore::Complex operator*(SciCore::Complex z, int n)
{
    return (z * static_cast<SciCore::Real>(n));
}

SCICORE_EXPORT inline SciCore::Complex operator/(SciCore::Complex z, int n)
{
    return (z / static_cast<SciCore::Real>(n));
}

SCICORE_EXPORT inline SciCore::Complex operator+(int n, SciCore::Complex z)
{
    return (static_cast<SciCore::Real>(n) + z);
}

SCICORE_EXPORT inline SciCore::Complex operator-(int n, SciCore::Complex z)
{
    return (static_cast<SciCore::Real>(n) - z);
}

SCICORE_EXPORT inline SciCore::Complex operator*(int n, SciCore::Complex z)
{
    return (static_cast<SciCore::Real>(n) * z);
}

SCICORE_EXPORT inline SciCore::Complex operator/(int n, SciCore::Complex z)
{
    return (static_cast<SciCore::Real>(n) / z);
}

SCICORE_EXPORT inline SciCore::Complex operator+(SciCore::Complex z, unsigned int n)
{
    return (z + static_cast<SciCore::Real>(n));
}

SCICORE_EXPORT inline SciCore::Complex operator-(SciCore::Complex z, unsigned int n)
{
    return (z - static_cast<SciCore::Real>(n));
}

SCICORE_EXPORT inline SciCore::Complex operator*(SciCore::Complex z, unsigned int n)
{
    return (z * static_cast<SciCore::Real>(n));
}

SCICORE_EXPORT inline SciCore::Complex operator/(SciCore::Complex z, unsigned int n)
{
    return (z / static_cast<SciCore::Real>(n));
}

SCICORE_EXPORT inline SciCore::Complex operator+(unsigned int n, SciCore::Complex z)
{
    return (static_cast<SciCore::Real>(n) + z);
}

SCICORE_EXPORT inline SciCore::Complex operator-(unsigned int n, SciCore::Complex z)
{
    return (static_cast<SciCore::Real>(n) - z);
}

SCICORE_EXPORT inline SciCore::Complex operator*(unsigned int n, SciCore::Complex z)
{
    return (static_cast<SciCore::Real>(n) * z);
}

SCICORE_EXPORT inline SciCore::Complex operator/(unsigned int n, SciCore::Complex z)
{
    return (static_cast<SciCore::Real>(n) / z);
}

SCICORE_EXPORT inline bool operator==(SciCore::Complex z, int n)
{
    return (z == static_cast<SciCore::Real>(n));
}

SCICORE_EXPORT inline bool operator!=(SciCore::Complex z, int n)
{
    return (z != static_cast<SciCore::Real>(n));
}

SCICORE_EXPORT inline bool operator==(int n, SciCore::Complex z)
{
    return (static_cast<SciCore::Real>(n) == z);
}

SCICORE_EXPORT inline bool operator!=(int n, SciCore::Complex z)
{
    return (static_cast<SciCore::Real>(n) != z);
}

SCICORE_EXPORT inline bool operator==(SciCore::Complex z, unsigned int n)
{
    return (z == static_cast<SciCore::Real>(n));
}

SCICORE_EXPORT inline bool operator!=(SciCore::Complex z, unsigned int n)
{
    return (z != static_cast<SciCore::Real>(n));
}

SCICORE_EXPORT inline bool operator==(unsigned int n, SciCore::Complex z)
{
    return (static_cast<SciCore::Real>(n) == z);
}

SCICORE_EXPORT inline bool operator!=(unsigned int n, SciCore::Complex z)
{
    return (static_cast<SciCore::Real>(n) != z);
}

namespace SciCore
{

///
/// \brief      Holds the result of a computation together with an error estimate.
///
/// \copybrief  ResultWithError
///
/// \tparam     T   Result type of a computation.
///
/// \headerfile Utility.h <SciCore/Utility.h>
///
template <typename T>
struct SCICORE_EXPORT ResultWithError
{
    /// \brief Result of the computation.
    T value;

    /// \brief Estimate of the absolute error.
    Real absErr;

    ResultWithError& operator=(const ResultWithError<std::decay_t<T>>& other)
    {
        value  = other.value;
        absErr = other.absErr;
        return *this;
    }

    ResultWithError& operator=(const ResultWithError<std::decay_t<T>&> other)
    {
        value  = other.value;
        absErr = other.absErr;
        return *this;
    }

    operator T() const&
    {
        return value;
    }

    operator T() && noexcept
    {
        return std::move(value);
    }
};

template <typename T>
SCICORE_EXPORT ResultWithError<T&> tie(T& value, Real absErr) noexcept
{
    return ResultWithError<T&>{value, absErr};
}

template <ScalarType T>
SCICORE_EXPORT bool isFinite(T x) noexcept
{
    if constexpr (std::is_same_v<std::decay_t<T>, Real> == true)
    {
        return std::isfinite(x);
    }
    else // Complex
    {
        return std::isfinite(x.real()) && std::isfinite(x.imag());
    }
}

template <DenseMatrixType MatrixT>
SCICORE_EXPORT bool isFinite(const MatrixT& x) noexcept
{
    return (x.array().isFinite().all() == true);
}

template <SparseMatrixType MatrixT>
SCICORE_EXPORT bool isFinite(const MatrixT& x) noexcept
{
    for (int k = 0; k < x.outerSize(); ++k)
    {
        for (Eigen::InnerIterator<std::decay_t<MatrixT>> it(x, k); it; ++it)
        {
            if (SciCore::isFinite(it.value()) == false)
            {
                return false;
            }
        }
    }

    return true;
}

template <typename T, typename A>
SCICORE_EXPORT bool isFinite(const std::vector<T, A>& x) noexcept
{
    for (size_t i = 0; i < x.size(); ++i)
    {
        if (isFinite(x[i]) == false)
        {
            return false;
        }
    }

    return true;
}

template <ScalarType T>
SCICORE_EXPORT Real maxNorm(T x) noexcept
{
    if constexpr (std::is_same_v<std::decay_t<T>, Real> == true)
    {
        return std::abs(x);
    }
    else // Complex
    {
        return std::max(std::abs(x.real()), std::abs(x.imag()));
    }
}

template <typename Derived>
SCICORE_EXPORT Real maxNorm(const Eigen::DenseBase<Derived>& x) noexcept
{
    return x.derived().template lpNorm<Eigen::Infinity>();
}

template <SparseMatrixType MatrixT>
SCICORE_EXPORT Real maxNorm(const MatrixT& x) noexcept
{
    Real normMax = 0;
    for (int k = 0; k < x.outerSize(); ++k)
    {
        for (Eigen::InnerIterator<std::decay_t<MatrixT>> it(x, k); it; ++it)
        {
            Real norm = std::abs(it.value());
            if (norm > normMax)
            {
                normMax = norm;
            }
        }
    }
    return normMax;
}

template <typename T, typename A>
SCICORE_EXPORT Real maxNorm(const std::vector<T, A>& vec) noexcept
{
    Real normMax = 0;
    for (const auto& x : vec)
    {
        Real norm = maxNorm(x);
        if (norm > normMax)
        {
            normMax = norm;
        }
    }
    return normMax;
}

SCICORE_EXPORT inline Real cwiseAbs(Real x) noexcept
{
    return maxNorm(x);
}

SCICORE_EXPORT inline Real cwiseAbs(Complex x) noexcept
{
    return maxNorm(x);
}

template <typename Derived>
SCICORE_EXPORT auto cwiseAbs(const Eigen::DenseBase<Derived>& x) noexcept
{
    return x.derived().cwiseAbs();
}

///
/// @brief      Returns the relative error between _testValue_ and _trueValue_.
///
/// Returns the relative error between _testValue_ and _trueValue_.
/// If _trueValue != 0_, then
///  \f[
///     \frac{|x_\text{test} - x_\text{true}|}{|x_\text{true}|}
///  \f]
/// is returned. Otherwise, _std::numeric_limits<Real>::max()_ is returned.
///
SCICORE_EXPORT Real relError(Real testValue, Real trueValue) noexcept;

///
/// @brief      Returns the relative error between the complex numbers _testValue_ and _trueValue_.
///
/// Returns the relative error between the complex numbers _testValue_ and _trueValue_.
/// If _trueValue != 0_, then
///  \f[
///     \frac{|x_\text{test} - x_\text{true}|}{|x_\text{true}|}
///  \f]
/// is returned. Otherwise, _std::numeric_limits<Real>::max()_ is returned.
///
SCICORE_EXPORT Real relError(Complex testValue, Complex trueValue) noexcept;

///
/// @brief      Returns the relative error between the matrices _testValue_ and _trueValue_.
///
/// Returns the relative error between the matrices _testValue_ and _trueValue_.
/// First, for each non-zero matrix element of _trueValue_ the relative error with the corresponding matrix element
/// of _testValue_ is computed. Out of all these computed relative errors the maximum is returned.
///
template <DenseMatrixType MatrixT>
SCICORE_EXPORT Real relError(const MatrixT& testValue, const MatrixT& trueValue) noexcept
{
    assert(testValue.rows() == trueValue.rows());
    assert(testValue.cols() == trueValue.cols());

    Real maxRelError = 0;
    for (int i = 0; i < testValue.rows(); ++i)
    {
        for (int j = 0; j < testValue.cols(); ++j)
        {
            if (trueValue(i, j) != Real(0))
            {
                Real relErr_ij = relError(testValue(i, j), trueValue(i, j));
                if (relErr_ij > maxRelError)
                {
                    maxRelError = relErr_ij;
                }
            }
        }
    }

    return maxRelError;
}

template <MatrixOrScalarType T>
SCICORE_EXPORT T cwiseQuotient(const T& A, const T& B, Real undefined)
{
    if constexpr (IsScalar<T>::value == true)
    {
        if (A == Real(0) && B == Real(0))
        {
            return undefined;
        }
        else
        {
            return A / B;
        }
    }
    else
    {
        int rows = A.rows();
        int cols = A.cols();
        assert(rows == B.rows());
        assert(cols == B.cols());

        T returnValue(rows, cols);
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                if (A(i, j) == Real(0) && B(i, j) == Real(0))
                {
                    returnValue(i, j) = undefined;
                }
                else
                {
                    returnValue(i, j) = A(i, j) / B(i, j);
                }
            }
        }
        return returnValue;
    }
}

template <MatrixOrScalarType T>
SCICORE_EXPORT bool isZero(const T& A, Real epsAbs = 0)
{
    return (maxNorm(A) <= epsAbs);
};

SCICORE_EXPORT constexpr inline Real truncToZero(Real x, Real precision)
{
    return (std::abs(x) >= precision) * x;
}

SCICORE_EXPORT inline Complex truncToZero(Complex x, Real precision)
{
    return Complex(truncToZero(x.real(), precision), truncToZero(x.imag(), precision));
}

template <DenseMatrixType MatrixT>
SCICORE_EXPORT MatrixT truncToZero(MatrixT&& x, Real precision)
{
    for (int i = 0; i < x.rows(); ++i)
    {
        for (int j = 0; j < x.cols(); ++j)
        {
            x(i, j) = truncToZero(x(i, j), precision);
        }
    }

    return x;
}

template <DenseMatrixType MatrixT>
SCICORE_EXPORT MatrixT truncToZero(const MatrixT& x, Real precision)
{
    MatrixT copy(x);
    return truncToZero(std::move(copy), precision);
}

SCICORE_EXPORT size_t nextPowerOf2(size_t n);

} // namespace SciCore

#endif // SCICORE_UTILITY_H
