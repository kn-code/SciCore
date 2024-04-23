//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

///
/// \file   Cheb.h
///
/// \brief  Chebyshev interpolation in one dimension.
///

#ifndef SCICORE_CHEB_H
#define SCICORE_CHEB_H

#include <concepts>

#include "DCT.h"
#include "Definitions.h"
#include "Serialization.h"
#include "Utility.h"

namespace SciCore
{

///
///  \defgroup Interpolation Chebyshev Interpolation
///
///  \brief Contains classes and methods concerning the Chebyshev interpolation of functions.
///
///  This module contains classes and methods concerning the Chebyshev interpolation of functions.
///
///  \{
///

///
///  \brief Computes _n_ Chebyshev nodes between _a_ and _b_.
///
///  The function computes _n_ Chebyshev nodes between _a_ and _b_.
///
///  \param  a          Lower interval point.
///  \param  b          Upper interval point.
///  \param  n          Number of Chebyshev nodes.
///  \param  nodes      A pointer to an array containing at least _n_ elements into which the result is written.
///
void chebNodes(Real a, Real b, int n, Real* nodes) noexcept;

///
///  \brief Computes _n_ Chebyshev nodes between _a_ and _b_ and returns the result.
///
///  The function computes _n_ Chebyshev nodes between _a_ and _b_ and returns the result.
///
///  \param  a          Lower interval point.
///  \param  b          Upper interval point.
///  \param  n          Number of Chebyshev nodes.
///
RealVector chebNodes(Real a, Real b, int n);

namespace Detail
{

template <MatrixOrScalarType T>
T chebyshevClenshawRecurrence(const T* c, int length, Real x) noexcept
{
    assert(length > 0);

    if (length < 2)
    {
        // if (length == 0)
        // {
        //     return 0;
        // }
        return Real(0.5) *  c[0];
    }

    T b1 = c[length - 1];
    T b2, tmp;

    // Initialize b2 to be the same dimensions as c, but all set to zero.
    if constexpr (IsDenseMatrix<T>::value == true)
    {
        b2 = T::Zero(c[0].rows(), c[0].cols());
    }
    else if constexpr (IsSparseMatrix<T>::value == true)
    {
        b2 = T(c[0].rows(), c[0].cols());
        b2.setZero();
    }
    else
    {
        b2 = 0.0;
    }

    for (int j = length - 2; j >= 1; --j)
    {
        tmp = 2 * x * b1 - b2 + c[j];

        b2 = b1;
        b1 = tmp;
    }

    return x * b1 - b2 + Real(0.5) * c[0];
}

// This is supposed to be more stable for x \approx 1 or x \approx -1.
// In these cases the normal chebyshevClenshawRecurrence can become numerically unstable (example?).
// See J. Oliver, An error analysis of the modified Clenshaw method for evaluating Chebyshev and Fourier series, J.
// Inst. Math. Appl., 20 (1977), pp. 379-391.
//? Example where this really happens?
template <MatrixOrScalarType T>
T modifiedChebyshevClenshawRecurrence(const T* c, int length, Real x) noexcept
{
    assert(length > 0);

    if (length < 2)
    {
        // if (length == 0)
        // {
        //     return 0;
        // }
        return c[0] / Real(2);
    }

    T b1 = c[length - 1];
    T d  = c[length - 1];
    T b2;

    // Initialize b2 to be the same dimensions as c, but all set to zero.
    if constexpr (IsDenseMatrix<T>::value == true)
    {
        b2 = T::Zero(c[0].rows(), c[0].cols());
    }
    else if constexpr (IsSparseMatrix<T>::value == true)
    {
        b2 = T(c[0].rows(), c[0].cols());
        b2.setZero();
    }
    else
    {
        b2 = 0.0;
    }

    if (x >= 0)
    {
        for (int j = length - 2; j >= 1; --j)
        {
            d  += (2 * (x - 1)) * b1 + c[j];
            b2  = b1;
            b1 += d;
        }
    }
    else
    {
        for (int j = length - 2; j >= 1; --j)
        {
            d  = (2 * (x + 1)) * b1 - d + c[j];
            b2 = b1;
            b1 = d - b2;
        }
    }

    return x * b1 - b2 + c[0] / Real(2);
}

template <MatrixOrScalarType T>
int chopChebSeriesAbsolute(const T* c, int N, Real epsAbs)
{
    // Enforce minimum size
    constexpr int minSize = 4;
    if (N < minSize)
    {
        return N;
    }

    RealVector envelope(N);
    envelope[N - 1] = maxNorm(c[N - 1]);
    for (int j = N - 2; j >= 0; --j)
    {
        envelope[j] = std::max(envelope[j + 1], maxNorm(c[j]));
    }

    if (envelope[0] == 0)
    {
        return 1;
    }

    for (int k = 0; k + 3 < N; ++k)
    {
        if (envelope[k] < epsAbs)
        {
            return k + 3;
        }
    }

    return N;
}

// http://arxiv.org/abs/1512.01803
template <MatrixOrScalarType T>
int chopChebSeriesRelative(const T* c, int N, Real epsRel)
{
    assert(N > 0);

    // Enforce minimum size
    constexpr int minSize = 15;
    if (N < minSize)
    {
        return N;
    }

    // Step 1: Create envelope
    RealVector envelope(N);
    if constexpr (IsScalar<T>::value == true)
    {
        envelope[N - 1] = std::abs(c[N - 1]);
        for (int j = N - 2; j >= 0; --j)
        {
            envelope[j] = std::max(envelope[j + 1], std::abs(c[j]));
        }

        if (envelope[0] == 0)
        {
            return 1;
        }
        envelope /= envelope[0];
    }
    else
    {
        std::vector<RealMatrix> envelopeMatrix(N);
        envelopeMatrix[N - 1] = c[N - 1].cwiseAbs();
        for (int j = N - 2; j >= 0; --j)
        {
            envelopeMatrix[j] = envelopeMatrix[j + 1].cwiseMax(c[j].cwiseAbs());
        }

        if (isZero(envelopeMatrix[0]) == true)
        {
            return 1;
        }

        for (int i = 1; i < N; ++i)
        {
            envelopeMatrix[i] = cwiseQuotient(envelopeMatrix[i], envelopeMatrix[0], 0);
        }
        envelopeMatrix[0] = cwiseQuotient(envelopeMatrix[0], envelopeMatrix[0], 0);

        for (int i = 0; i < N; ++i)
        {
            envelope[i] = maxNorm(envelopeMatrix[i]);
        }
    }

    // Step 2: Scan envelope for plateau point.
    // Plateau is located at c[plateauPoint+1], c[plateauPoint+2], ... c[N-1]
    int plateauPoint = -1;
    int j2           = -1;
    for (int j = 1; j < N; ++j)
    {
        j2 = std::round(1.25 * j + 5.25);
        if (j2 > N - 1)
        {
            // This is not in original algorithm:
            // Check if envelope[k] <= epsRel for any k=j...N-3 ?
            // In this case we would always have r <= 0, and therefore
            // the plateauPoint would be set to k - 1.
            // Since j2 is at this point in the code undefined, we skip
            // step 3 and return immidiately.
            for (int k = j; k + 3 < N; ++k)
            {
                if (envelope[k] < epsRel)
                {
                    return k + 3;
                }
            }
            return N;
        }

        Real e1 = envelope[j];
        Real e2 = envelope[j2];
        Real r  = 3 * (1 - std::log(e1) / std::log(epsRel));

        if ((e1 == 0) || (e2 / e1 > r))
        {
            plateauPoint = j - 1;
            break;
        }
    }

    // Step 3: Decide precise cutoff point, see paper
    if (envelope[plateauPoint] == 0)
    {
        return plateauPoint + 1;
    }

    Real tol = std::pow(epsRel, 7. / 6.);
    int j3   = -1;
    for (Real x : envelope)
    {
        j3 += static_cast<int>(x >= tol);
    }
    if (j3 < j2)
    {
        j2           = j3 + 1;
        envelope[j2] = tol;
    }

    RealVector cc  = envelope.head(j2 + 1).array().log10();
    cc            += RealVector::LinSpaced(j2 + 1, 0, Real(-1. / 3.) * std::log10(epsRel));

    int d = -1;
    cc.minCoeff(&d);

    return std::max(d, 1);
}

} // namespace Detail

///
/// \brief Computes Chebyshev interpolations for scalar- and matrix-valued functions of a real parameter.
///
/// The class computes interpolations at the zeros of the Chebyshev polynomials for scalar- and matrix-valued functions
/// of a real parameter. Chebyshev interpolations are a powerful tool to interpolate functions \f$f(x)\f$ of a single
/// real variable. This is especially useful if \f$f(x)\f$ is hard to compute: instead of recomputing for different
/// arguments \f$x\f$ we only compute once and then interpolate. The basic usage is as follows:
///
/// \snippet tests/ChebTest.cpp Basic cheb usage
///
/// \headerfile Cheb.h <SciCore/Cheb.h>
///
template <MatrixOrScalarType T>
class Cheb
{
  public:
    ///
    /// \brief The return type of the interpolated function (Real, Complex, Matrix, ...)
    ///
    using ElementType = T;

    ///
    /// \brief Constructs an empty Chebyshev object.
    ///
    Cheb() noexcept
    {
    }

    ///
    /// \brief  Creates a Chebyshev interpolation of the function \f$f(x)\f$ with \f$x\in[a,b]\f$ using \f$n\f$ points.
    ///
    /// Creates a Chebyshev interpolation of the function \f$f(x)\f$ with \f$x\in[a,b]\f$ using \f$n\f$ points.
    /// The function \f$f(x)\f$ can either return a scalar (real or complex) or a matrix (dense or sparse). Note that
    /// \f$n\f$ is required to be a power of 2.
    ///
    /// During the computation a _buffer_ of \f$n\f$ matrices is needed, which can be supplied as an additional argument
    /// (to prevent unnecessary allocations). This is useful if several Chebyshev interpolations have to be computed in
    /// sequence. If _buffer_ is _NULL_ it will be allocated automatically.
    ///
    /// \param  f          The function for which the Chebyshev approximation is computed.
    /// \param  a          Lower interval point.
    /// \param  b          Upper interval point.
    /// \param  n          Order of the Chebyshev approximation.
    /// \param  buffer     A buffer of _n_ objects of type _T_.
    ///
    /// \tparam FunctionT  The type of function. Has to take a \a Real and return a \a T.
    ///
    template <typename FunctionT>
        requires std::invocable<FunctionT, Real>
    Cheb(FunctionT&& f, Real a, Real b, int n, T* buffer = nullptr) : _coeffs((size_t)n), _a(a), _b(b)
    {
        assert(n > 0);

        if (buffer != nullptr)
        {
            _computeFixedN(std::forward<FunctionT>(f), buffer);
        }
        else
        {
            std::vector<T> newBuffer((size_t)n);
            _computeFixedN(std::forward<FunctionT>(f), newBuffer.data());
        }
    }

    Cheb(const T* f, Real a, Real b, int n, T* buffer = nullptr) : _coeffs((size_t)n), _a(a), _b(b)
    {
        assert(n > 0);

        if (buffer != nullptr)
        {
            _computeFixedN(f, n, buffer);
        }
        else
        {
            std::vector<T> newBuffer((size_t)n);
            _computeFixedN(f, n, newBuffer.data());
        }
    }

    ///
    /// \brief      Computes a Chebyshev interpolation of \f$f(x)\f$ using an adaptive algorithm to determine the number
    /// of needed points.
    ///
    /// Computes a Chebyshev interpolation of \f$f(x)\f$ using an adaptive algorithm to determine the number
    /// of needed points.
    /// At first \f$f\f$ is interpolated using _nStart_ points. Then the error is estimated. If the error is too big,
    /// then a new interpolation using _3*nStart_ points is constructed. It is used that the Chebyshev points of order
    /// _n_ and _3*n_ are nested, which means that all previous computations of \f$f\f$ can be reused. The number of
    /// points is tripled at most _nTrip_ times.
    ///
    /// \param      f               The function for which the Chebyshev approximation is computed.
    /// \param      a               Lower interval point.
    /// \param      b               Upper interval point.
    /// \param      epsAbs          Absolute error goal.
    /// \param      epsRel          Relative error goal.
    /// \param      nStart          Order of the initial Chebyshev approximation.
    /// \param      nTrip           The maximum number of triplings. Can currently at most be equal to 2.
    /// \param      ok              Set to _true_ if error goal was achieved, otherwise set to _false_.
    /// \param      buffer          A buffer containing \f$n_{\text{Start}} 3^{n_{\text{Trip}}}\f$ elements.
    ///                             If set to _NULL_ then the buffer is allocated automatically.
    ///
    /// \tparam     FunctionT  The type of function. Has to take a \a Real and return a \a T.
    ///
    template <typename FunctionT>
        requires std::invocable<FunctionT, Real>
    Cheb(
        FunctionT&& f,
        Real a,
        Real b,
        Real epsAbs,
        Real epsRel,
        int nStart,
        int nTrip = 2,
        bool* ok  = nullptr,
        T* buffer = nullptr)
        : _a(a), _b(b)
    {
        assert(nStart >= 2);
        assert(nTrip >= 0);

        bool errorOk;
        if (buffer == nullptr)
        {
            std::vector<T> newBuffer((size_t)(nStart * pow(3, nTrip)));
            errorOk =
                _computeAdaptive(std::forward<FunctionT>(f), epsAbs, epsRel, nStart, nTrip, newBuffer.data());
        }
        else
        {
            errorOk = _computeAdaptive(std::forward<FunctionT>(f), epsAbs, epsRel, nStart, nTrip, buffer);
        }

        if (ok != nullptr)
        {
            *ok = errorOk;
        }
    }

    ///
    /// \brief Copy constructor.
    ///
    Cheb(const Cheb& other) noexcept : _coeffs(other._coeffs), _a(other._a), _b(other._b)
    {
    }

    ///
    /// \brief Move constructor.
    ///
    Cheb(Cheb&& other) noexcept : _coeffs(std::move(other._coeffs)), _a(other._a), _b(other._b)
    {
    }

    ///
    /// \brief Copy assignment operator.
    ///
    Cheb& operator=(const Cheb& other) noexcept
    {
        _coeffs = other._coeffs;
        _a      = other._a;
        _b      = other._b;
        return *this;
    }

    ///
    /// \brief Move assignment operator.
    ///
    Cheb& operator=(Cheb&& other) noexcept
    {
        _coeffs = std::move(other._coeffs);
        _a      = other._a;
        _b      = other._b;
        return *this;
    }

    ///
    /// \brief Evaluates the Chebyshev series at \f$x\f$ and returns the result.
    ///
    T operator()(Real x) const
    {
        int n = std::ssize(_coeffs);
        if (n <= 0)
        {
            throw std::runtime_error("Cheb::operator(): Can't evaluate empty Cheb object.");
        }

        Real z = (Real(2) * x - _a - _b) / (_b - _a);
        if (z < Real(-0.95) || z > Real(0.95))
        {
            return Detail::modifiedChebyshevClenshawRecurrence(_coeffs.data(), n, z);
        }
        else
        {
            return Detail::chebyshevClenshawRecurrence(_coeffs.data(), n, z);
        }
    }

    ///
    /// \brief Evaluates the Chebyshev series at at all points in \a x. Only works if \a ElementType is \a Real.
    ///
    RealVector operator()(const RealVector& x) requires std::is_same_v<std::decay_t<T>, Real>
    {
        RealVector returnValue(x.size());
        for (int i = 0; i < x.size(); ++i)
        {
            returnValue[i] = this->operator()(x[i]);
        }
        return returnValue;
    }

    ///
    /// \brief Evaluates the Chebyshev series at at all points in \a x. Only works if \a ElementType is \a Complex.
    ///
    Vector operator()(const RealVector& x) requires std::is_same_v<std::decay_t<T>, Complex>
    {
        Vector returnValue(x.size());
        for (int i = 0; i < x.size(); ++i)
        {
            returnValue[i] = this->operator()(x[i]);
        }
        return returnValue;
    }

    ///
    /// \brief      Reduces the number of used Chebyshev coefficients with the goal to achieve a given error goal.
    ///
    /// \copybrief  chopCoefficients(Real epsAbs, Real epsRel)
    /// Reducing the number of used coefficients can increase the evaluation speed. Returns _true_ if the error goal was likely
    /// met, and _false_ otherwise. The used algorithm to decide the cutoff point is described in
    /// <a href="http://arxiv.org/abs/1512.01803"> J. L. Aurentz and L. N. Trefethen, "Chopping a Chebyshev series"</a>.
    ///
    /// \param epsAbs       The absolute error goal.
    /// \param epsRel       The relative error goal.
    ///
    bool chopCoefficients(Real epsAbs, Real epsRel)
    {
        int n       = std::ssize(_coeffs);
        int chopAbs = ((epsAbs > 0) ? Detail::chopChebSeriesAbsolute(_coeffs.data(), n, epsAbs) : n);
        int chopRel = ((epsRel > 0) ? Detail::chopChebSeriesRelative(_coeffs.data(), n, epsRel) : n);
        int chop    = std::min(chopAbs, chopRel);

        if (chop >= n)
        {
            return false;
        }
        else
        {
            if (n > 4)
            {
                _coeffs.resize(std::max(chop, 4));
            }

            return true;
        }
    }

    ///
    /// \brief  Computes the derivative of the represented function and returns it as a new Cheb object.
    ///
    /// \copybrief diff
    ///
    /// \warning Computing the derivative is numerically not very stable,
    ///          see also <a href="https://www.sciencedirect.com/science/article/pii/0021999192902743">
    ///          On the errors incurred calculating derivatives using Chebyshev polynomials.</a>
    ///
    Cheb diff() const
    {
        int n = std::ssize(_coeffs);

        Cheb returnValue;
        returnValue._a = _a;
        returnValue._b = _b;
        returnValue._coeffs.resize((size_t)n);

        Real con = 2 / (_b - _a);

        if (n > 0)
        {
            if constexpr (IsDenseMatrix<T>::value == true)
            {
                returnValue._coeffs[n - 1] = T::Zero(_coeffs[0].rows(), _coeffs[0].cols());
            }
            else if constexpr (IsSparseMatrix<T>::value == true)
            {
                returnValue._coeffs[n - 1] = T(_coeffs[0].rows(), _coeffs[0].cols());
                returnValue._coeffs[n - 1].setZero();
            }
            else // is Real or Complex
            {
                returnValue._coeffs[n - 1] = 0;
            }
        }

        if (n > 1)
        {
            returnValue._coeffs[n - 2] = Real(2 * (n - 1)) * _coeffs[n - 1];

            for (int i = n; i >= 3; --i)
            {
                returnValue._coeffs[i - 3] = returnValue._coeffs[i - 1] + Real(2 * (i - 2)) * _coeffs[i - 2];
            }

            for (int i = 0; i < n; ++i)
            {
                returnValue._coeffs[i] *= con;
            }
        }

        return returnValue;
    }

    ///
    /// \brief  Returns the integral of the represented function.
    ///
    /// \copybrief integrate
    /// If the object represents the function \f$f(x)\f$, then the antiderivative \f$ F(x) \f$, defined by
    /// \f[
    ///     F(x) = \int_{a}^x f(x') dx',
    /// \f]
    /// is returned where \a a corresponds to lowerLimit().
    ///
    Cheb integrate() const
    {
        int n = std::ssize(_coeffs);
        Cheb returnValue;
        returnValue._a = _a;
        returnValue._b = _b;
        returnValue._coeffs.resize((size_t)n);

        Real con = (_b - _a) / 4;

        if (n == 1)
        {
            if constexpr (IsMatrix<T>::value == true)
            {
                returnValue._coeffs[0].setZero();
            }
            else // is Real or Complex
            {
                returnValue._coeffs[0] = Real(0);
            }
        }
        else if (n == 2)
        {
            returnValue._coeffs[1] = con * _coeffs[0];
            returnValue._coeffs[0] = Real(2) * returnValue._coeffs[1];
        }
        else
        {
            T sum;
            if constexpr (IsDenseMatrix<T>::value == true)
            {
                sum = T::Zero(_coeffs[0].rows(), _coeffs[0].cols());
            }
            else if constexpr (IsSparseMatrix<T>::value == true)
            {
                sum = T(_coeffs[0].rows(), _coeffs[0].cols());
                sum.setZero();
            }
            else // is Real or Complex
            {
                sum = 0.0;
            }

            Real fac = 1.0;
            for (int i = 1; i <= n - 2; ++i)
            {
                returnValue._coeffs[i]  = (con / ((Real)i)) * (_coeffs[i - 1] - _coeffs[i + 1]);
                sum                    += fac * returnValue._coeffs[i];
                fac                     = -fac;
            }
            returnValue._coeffs[n - 1]  = (con / (n - 1)) * _coeffs[n - 2];
            sum                        += fac * returnValue._coeffs[n - 1];
            returnValue._coeffs[0]      = Real(2) * sum;
        }

        return returnValue;
    }

    ///
    /// \brief  Adds the function represented by \a other to this object.
    ///
    /// \copybrief Cheb::operator+=
    /// It is necessary that this object and \a other have the same lower limit and upper limit,
    /// otherwise the operation is undefined.
    ///
    /// \param  other  The object that is added.
    ///
    Cheb& operator+=(const Cheb& other) noexcept
    {
        int sizeMe    = std::ssize(_coeffs);
        int sizeOther = std::ssize(other._coeffs);

        if (sizeMe == 0)
        {
            return operator=(other);
        }

        if (sizeOther == sizeMe)
        {
            for (int i = 0; i < sizeMe; ++i)
            {
                _coeffs[i] += other._coeffs[i];
            }
        }
        else if (sizeOther > sizeMe)
        {
            _coeffs.resize(sizeOther);
            for (int i = 0; i < sizeMe; ++i)
            {
                _coeffs[i] += other._coeffs[i];
            }
        }
        else
        {
            for (int i = 0; i < sizeOther; ++i)
            {
                _coeffs[i] += other._coeffs[i];
            }
        }

        return *this;
    }

    ///
    /// \brief  Adds the constant _constant_ to the function represented by this object.
    ///
    Cheb& operator+=(const T& constant)
    {
        int n = std::ssize(_coeffs);
        if (n == 0)
        {
            throw std::runtime_error("Cheb::operator+=(const T& ): Can't use on empty Cheb object.");
        }
        else
        {
            _coeffs[0] += Real(2) * constant;
        }

        return *this;
    }

    ///
    /// \brief  Multiplies the represented function by a scalar _a_.
    ///
    Cheb& operator*=(Real a) noexcept
    {
        for (auto& c : _coeffs)
        {
            c *= a;
        }
        return *this;
    }

    ///
    /// \brief Test if two objects are the same.
    ///
    bool operator==(const Cheb& other) const noexcept
    {
        return _coeffs == other._coeffs && _a == other._a && _b == other._b;
    }

    ///
    /// \brief Test if two objects are not the same.
    ///
    bool operator!=(const Cheb& other) const noexcept
    {
        return !(*this == other);
    }

    ///
    /// \brief  Returns the Chebyshev coefficients of the approximation.
    ///
    std::vector<T>& coefficients() noexcept
    {
        return _coeffs;
    }

    ///
    /// \brief  Returns the Chebyshev coefficients of the approximation.
    ///
    const std::vector<T>& coefficients() const noexcept
    {
        return _coeffs;
    }

    ///
    /// \brief  Returns the lower limit of the approximation.
    ///
    Real lowerLimit() const noexcept
    {
        return _a;
    }

    ///
    /// \brief  Returns the upper limit of the approximation.
    ///
    Real upperLimit() const noexcept
    {
        return _b;
    }

    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(CEREAL_NVP(_coeffs), CEREAL_NVP(_a), CEREAL_NVP(_b));
    }

  private:
    std::vector<T> _coeffs;
    Real _a = -1;
    Real _b = 1;

    // Computes the first n Chebyshev coefficients of the function f,
    // where n is _coeffs.size() and buffer contains n elements.
    template <typename FunctionT>
        requires std::invocable<FunctionT, Real>
    void _computeFixedN(FunctionT&& f, T* buffer) noexcept
    {
        Real bma = (_b - _a) / 2;
        Real bpa = (_b + _a) / 2;

        int n      = std::ssize(_coeffs);
        Real inv_n = 1 / static_cast<Real>(n);

        for (int j = n - 1; j >= 0; --j) // This order makes sure that f is sampled at increasing values
        {
            Real y     = std::cos(std::numbers::pi_v<Real> * (j + Real(0.5)) * inv_n);
            _coeffs[j] = (2 * inv_n) * f(std::fma(y, bma, bpa));
        }

        dct(_coeffs.data(), n, buffer);
    }

    void _computeFixedN(const T* f, int n, T* buffer)
    {
        Real inv_n = 1 / static_cast<Real>(n);

        for (int j = n - 1; j >= 0; --j)
        {
            _coeffs[j] = (2 * inv_n) * f[n - j - 1];
        }

        dct(_coeffs.data(), n, buffer);
    }

    // Buffer needs to contain nStart*(nTrip+1) elements
    template <typename FunctionT>
    bool _computeAdaptive(FunctionT&& f, Real epsAbs, Real epsRel, int nStart, int nTrip, T* buffer)
    {
        assert(nStart >= 2);
        assert(nTrip >= 0);

        Real bma = (_b - _a) / 2;
        Real bpa = (_b + _a) / 2;

        //
        // Compute first approximation
        //
        int n = nStart;
        std::vector<T> newCoeff(n);
        _coeffs    = std::vector<T>(n);
        Real inv_n = 1 / static_cast<Real>(n);
        for (int j = 0; j < n; ++j)
        {
            Real y      = std::cos(std::numbers::pi_v<Real> * (j + Real(0.5)) * inv_n);
            newCoeff[j] = (2 * inv_n) * f(std::fma(y, bma, bpa));
        }
        _coeffs = newCoeff;
        dct(_coeffs.data(), n, buffer);

        if (chopCoefficients(epsAbs, epsRel) == true)
        {
            return true;
        }

        //
        // Else if initial approximation is not good enough:
        // improve by tripling number of points
        //
        for (int iterTripling = 0; iterTripling < nTrip; ++iterTripling)
        {
            _coeffs  = newCoeff; // Restore original data points
            n       *= 3;
            inv_n    = 1 / static_cast<Real>(n);
            newCoeff.resize(n);
            Real y              = std::cos(std::numbers::pi_v<Real> / 2 * inv_n);
            newCoeff[0]         = (2 * inv_n) * f(std::fma(y, bma, bpa));
            size_t oldDataIndex = 0;
            for (int j = 1; j < n; ++j)
            {
                if ((j - 1) % 3 == 0)
                {
                    newCoeff[j] = _coeffs[oldDataIndex] / Real(3);
                    ++oldDataIndex;
                }
                else
                {
                    y           = std::cos(std::numbers::pi_v<Real> * (j + Real(0.5)) * inv_n);
                    newCoeff[j] = (2 * inv_n) * f(std::fma(y, bma, bpa));
                }
            }

            _coeffs = newCoeff;
            dct(_coeffs.data(), n, buffer);

            if (chopCoefficients(epsAbs, epsRel) == true)
            {
                return true;
            }
        }

        return false;
    }
};

template <typename FunctionT>
Cheb(FunctionT&& f, Real a, Real b, int n) -> Cheb<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
Cheb(FunctionT&& f, Real a, Real b, Real epsAbs, Real epsRel, int nStart)
    -> Cheb<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
Cheb(FunctionT&& f, Real a, Real b, Real epsAbs, Real epsRel, int nStart, int nTrip)
    -> Cheb<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
Cheb(FunctionT&& f, Real a, Real b, Real epsAbs, Real epsRel, int nStart, int nTrip, bool* ok)
    -> Cheb<std::invoke_result_t<FunctionT, Real>>;

#ifndef SCICORE_DONT_PRECOMPILE_TEMPLATES
extern template class Cheb<Real>;
extern template class Cheb<Complex>;
extern template class Cheb<RealVector>;
extern template class Cheb<Vector>;
extern template class Cheb<RealMatrix>;
extern template class Cheb<Matrix>;
#endif // SCICORE_DONT_PRECOMPILE_TEMPLATES

/// \} // end of Interpolation

} // namespace SciCore

#endif // SCICORE_CHEB_H
