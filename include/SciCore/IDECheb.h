//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

///
/// \file   IDECheb.h
///
/// \brief  Solution of integro-differential equations using Chebyshev methods.
///

#ifndef SCICORE_IDE_CHEB_H
#define SCICORE_IDE_CHEB_H

#include <Eigen/LU>

#include "ChebAdaptive.h"
#include "Definitions.h"
#include "Integration.h"

namespace SciCore
{

namespace Detail
{

//
// \brief  Calculates the convolution matrix as presented in
//         "SPECTRAL APPROXIMATION OF CONVOLUTION OPERATORS, KUAN XU and ANA F. LOUREIRO".
//
// \param  chebCoefficients   The Chebyshev coefficients of the memory kernel.
// \param  tMax_minus_t0
//
// \tparam VectorT            A \a Real or \a Complex vector type.
//
template <MatrixOrScalarType T>
Eigen::Matrix<typename GetScalarType<T>::type, Eigen::Dynamic, Eigen::Dynamic> convolutionMatrix(
    const std::vector<T>& chebCoefficients,
    Real tMax_minus_t0)
{
    int n = chebCoefficients.size();
    assert(n >= 4);

    int rows = 1;
    int cols = 1;
    if constexpr (IsMatrix<T>::value == true)
    {
        rows = chebCoefficients[0].rows();
        cols = chebCoefficients[0].cols();
    }

    using ElementType   = typename GetScalarType<T>::type;
    using DynamicMatrix = Eigen::Matrix<ElementType, Eigen::Dynamic, Eigen::Dynamic>;
    DynamicMatrix convMat(rows * n, cols * n);
    convMat.setZero();

    auto getValue = [&convMat, rows,
                     cols](int r, int c) -> std::conditional_t<IsScalar<T>::value, T&, Eigen::Block<DynamicMatrix>>
    {
        if constexpr (IsScalar<T>::value == true)
        {
            return convMat(r, c);
        }
        else
        {
            return convMat.block(r * rows, c * cols, rows, cols);
        }
    };

    //
    // Compute zeroth column (3.1)
    // Carful: we have different convention than the paper
    //         such that a_0 in the paper is chebCoefficients[0]/2.
    //

    getValue(1, 0) = (chebCoefficients[0] - chebCoefficients[2]) / 2.0;
    for (int k = 2; k < n - 1; ++k)
    {
        getValue(k, 0) = (chebCoefficients[k - 1] - chebCoefficients[k + 1]) / (2.0 * k);
    }
    getValue(n - 1, 0) = chebCoefficients[n - 2] / (2.0 * (n - 1));

    Real pm = 1; // plus or minus one
    for (int j = 1; j < n; ++j)
    {
        getValue(0, 0) += pm * getValue(j, 0);
        pm             *= -1;
    }

    //
    // Compute first column (3.2a) on and below main diagonal
    //
    getValue(1, 1) = -getValue(1, 0) + getValue(0, 0) - getValue(2, 0) / 2.0;
    for (int k = 2; k < n - 1; ++k)
    {
        getValue(k, 1) = -getValue(k, 0) + getValue(k - 1, 0) / (2.0 * k) - getValue(k + 1, 0) / (2.0 * k);
    }
    getValue(n - 1, 1) = -getValue(n - 1, 0) + getValue(n - 2, 0) / (2.0 * (n - 1));

    //
    // Compute second column (3.2b) on and below main diagonal
    //
    for (int k = 2; k < n - 1; ++k)
    {
        getValue(k, 2) = getValue(k, 0) + getValue(k - 1, 1) * Real(2.0 / k) - getValue(k + 1, 1) * Real(2.0 / k);
    }
    getValue(n - 1, 2) = getValue(n - 1, 0) + getValue(n - 2, 1) * (2.0 / (n - 1));

    //
    // Compute all other elements (3.2c) on and below main diagonal
    //
    for (int m = 2; m + 1 < n - 1; ++m)
    {
        for (int k = m + 1; k < n - 1; ++k)
        {
            pm = (m % 2 == 0) ? 1 : -1;
            getValue(k, m + 1) =
                2 * pm / Real(m - 1) * getValue(k, 0) + Real(m + 1) / Real(m - 1) * getValue(k, m - 1) +
                Real(m + 1) / Real(k) * getValue(k - 1, m) - Real(m + 1) / Real(k) * getValue(k + 1, m);
        }
        pm                     = (m % 2 == 0) ? 1 : -1;
        getValue(n - 1, m + 1) = 2 * pm / Real(m - 1) * getValue(n - 1, 0) +
                                 Real(m + 1) / Real(m - 1) * getValue(n - 1, m - 1) +
                                 Real(m + 1) / Real(n - 1) * getValue(n - 2, m);
    }
    pm                     = ((n - 2) % 2 == 0) ? 1 : -1;
    getValue(n - 1, n - 1) = 2 * pm / Real(n - 3) * getValue(n - 1, 0) +
                             Real(n - 1) / Real(n - 3) * getValue(n - 1, n - 3) +
                             Real(n - 1) / Real(n - 1) * getValue(n - 2, n - 2);

    //
    // Compute elements above the main diagonal using (3.4) from bottom to top
    //
    pm                     = ((n - 1) % 2 == 0) ? 1 : -1;
    getValue(n - 2, n - 1) = -2 * (n - 1) * pm / Real((n - 1) * (n - 1) - 1) * getValue(n - 1, 0) -
                             (n - 1) / Real(n - 2) * getValue(n - 1, n - 2);
    for (int k = n - 2; k > 1; --k)
    {
        for (int m = k; m < n - 1; ++m)
        {
            pm                 = (m % 2 == 0) ? 1 : -1;
            getValue(k - 1, m) = -2 * k * pm / Real(m * m - 1) * getValue(k, 0) - k / Real(m - 1) * getValue(k, m - 1) +
                                 k / Real(m + 1) * getValue(k, m + 1) + getValue(k + 1, m);
        }
        getValue(k - 1, n - 1) = -2 * k * pm / Real((n - 1) * (n - 1) - 1) * getValue(k, 0) -
                                 k / Real(n - 2) * getValue(k, n - 2) + getValue(k + 1, n - 1);
    }

    for (int m = 2; m < n - 1; ++m)
    {
        pm             = (m % 2 == 0) ? 1 : -1;
        getValue(0, m) = 0.5 * (-2 * pm / Real(m * m - 1) * getValue(1, 0) - 1 / Real(m - 1) * getValue(1, m - 1) +
                                1 / Real(m + 1) * getValue(1, m + 1) + getValue(2, m));
    }
    getValue(0, n - 1) = 0.5 * (-2 * pm / Real((n - 1) * (n - 1) - 1) * getValue(1, 0) -
                                1 / Real(n - 2) * getValue(1, n - 2) + getValue(2, n - 1));

    pm = 1;
    for (int j = 1; j < n; ++j)
    {
        getValue(0, 1) += pm * getValue(j, 1);
        pm             *= -1;
    }

    for (int i = 1; i < n; ++i)
    {
        getValue(i, 0) /= Real(2);
        getValue(0, i) *= Real(2);
    }

    convMat *= tMax_minus_t0 / Real(2);

    return convMat;
}

template <ScalarType ElementType, MatrixOrScalarType T>
void process(
    Eigen::Matrix<ElementType, Eigen::Dynamic, Eigen::Dynamic>& convolutionMatrix,
    const T& L,
    Real tMax_minus_t0)
{
    int dimT = 1;
    if constexpr (IsMatrix<T>::value == true)
    {
        dimT = L.rows();
        assert(L.cols() == dimT);
    }

    using DynamicMatrix = Eigen::Matrix<ElementType, Eigen::Dynamic, Eigen::Dynamic>;
    auto getValue       = [&convolutionMatrix,
                     dimT](int r, int c) -> std::conditional_t<IsScalar<T>::value, T&, Eigen::Block<DynamicMatrix>>
    {
        if constexpr (IsScalar<T>::value == true)
        {
            return convolutionMatrix(r, c);
        }
        else
        {
            return convolutionMatrix.block(r * dimT, c * dimT, dimT, dimT);
        }
    };

    int n = convolutionMatrix.rows();
    assert(convolutionMatrix.cols() == n);
    if constexpr (IsMatrix<T>::value == true)
    {
        n /= L.rows();
        assert(n * L.rows() == convolutionMatrix.rows());
    }
    assert(n >= 2);

    Real factor     = 2 / tMax_minus_t0;
    int startValue  = 2;
    int startColumn = 1;
    for (int row = 0; row < n; ++row)
    {
        int value = startValue;
        for (int col = startColumn; col < n; col += 2)
        {
            if constexpr (IsScalar<T>::value == true)
            {
                getValue(row, col) -= factor * value;
            }
            else
            {
                getValue(row, col) -= factor * value * T::Identity(dimT, dimT);
            }
            value += 4;
        }
        startValue += 2;
        ++startColumn;
    }

    // Add Liouvillian condition
    for (int row = 0; row < n; ++row)
    {
        getValue(row, row) += L;
    }

    // Write initial condition
    if constexpr (IsScalar<T>::value == true)
    {
        getValue(n - 1, 0) = Real(0.5);
    }
    else
    {
        getValue(n - 1, 0) = Real(0.5) * T::Identity(dimT, dimT);
    }
    factor = -1;
    for (int col = 1; col < n; ++col)
    {
        if constexpr (IsScalar<T>::value == true)
        {
            getValue(n - 1, col) = factor;
        }
        else
        {
            getValue(n - 1, col) = factor * T::Identity(dimT, dimT);
        }
        factor *= -1;
    }
}

} // namespace Detail

// \int_{t_0}^t ds k(t-s) f(s)
template <MatrixOrScalarType T, MatrixOrScalarType T1, MatrixOrScalarType T2>
Cheb<T> convolution(const Cheb<T1>& k, const Cheb<T2>& f)
{
    int N = k.coefficients().size();
    if ((int)f.coefficients().size() != N)
    {
        throw std::runtime_error(
            "convolution: currently only implemented if arguments have same number of Chebyshev coefficients.");
    }

    if (N < 4)
    {
        throw std::runtime_error("convolution: kernel needs at least 4 Chebyshev coefficients.");
    }

    if (k.lowerLimit() != 0)
    {
        throw std::runtime_error("convolution: kernel must have lower limit zero.");
    }

    if (std::abs(k.upperLimit() - (f.upperLimit() - f.lowerLimit())) > 10 * std::numeric_limits<Real>::epsilon())
    {
        throw std::runtime_error("convolution: arguments have inconsistent limits.");
    }

    int rows = 1;
    int cols = 1;
    if constexpr (IsMatrix<T>::value == true)
    {
        rows = k.coefficients()[0].rows();
        cols = k.coefficients()[0].cols();
    }

    auto convMatrix = Detail::convolutionMatrix(k.coefficients(), f.upperLimit() - f.lowerLimit());
    auto getValue   = [&convMatrix, rows, cols](int r, int c)
    {
        if constexpr (IsScalar<T>::value == true)
        {
            return convMatrix(r, c);
        }
        else
        {
            return convMatrix.block(r * rows, c * cols, rows, cols);
        }
    };
    Cheb<T> returnValue(
        [&](Real) -> T
        {
            if constexpr (IsScalar<T>::value == true)
            {
                return Real(0);
            }
            else
            {
                return T::Zero(k.coefficients()[0].rows(), f.coefficients()[0].cols());
            }
        },
        f.lowerLimit(), f.upperLimit(), N);

    for (int i = 0; i < N; ++i)
    {
        returnValue.coefficients()[i] = getValue(i, 0) * f.coefficients()[0];
        for (int j = 1; j < N; ++j)
        {
            returnValue.coefficients()[i] += getValue(i, j) * f.coefficients()[j];
        }
    }
    return returnValue;
}

///
/// \ingroup IDE
///
/// @brief          Solves integro-differential equations using Chebyshev methods.
///
/// This method solves integro-differential equations of the form
/// \f[
///     \dot y(t) = g y(t) + \int_{t_0}^t ds k(t-s) y(s) + h(t)
/// \f]
/// with initial value \f$y(t_0)=y_0\f$.
/// Here \f$y(t)\f$ can be (real or complex) scalar- or vector-valued (of type _TSolution_). Thus the memory kernel _k_
/// must correspondingly be also a (real or complex) scalar- or matrix-valued function (of type _TKernel_).
///
/// @param g            A constant (matrix or scalar).
/// @param k            The memory kernel.
/// @param h            Inhomogeneous part.
/// @param y0           Initial value.
/// @param t0           Initial time.
/// @param tMax         Time until the solution is computed.
/// @param nMinCheb     Minimum number of Chebyshev coefficients used to represent the solution. By default the
///                     same number of Chebyshev points are used as the number of points to represent the
///                     memory kernel.
///
/// @tparam TSolution   Type of the solution.
/// @tparam TKernel     Type of the memory kernel.
///
template <MatrixOrScalarType TSolution, MatrixOrScalarType TKernel>
Cheb<TSolution> solveIdeCheb(
    const TKernel& g,
    const Cheb<TKernel>& k,
    const Cheb<TSolution>& h,
    const TSolution& y0,
    Real t0,
    Real tMax,
    int nMinCheb = 0)
{
    // Do some consistency checks
    if constexpr (IsScalar<TSolution>::value == false)
    {
        static_assert(TSolution::ColsAtCompileTime == 1, "Solution must be a scalar or vector type.");
    }

    if (k.lowerLimit() != 0)
    {
        throw std::runtime_error("solveIdeCheb: kernel must have lower limit zero.");
    }

    if (std::abs(k.upperLimit() - (tMax - t0)) > 10 * std::numeric_limits<Real>::epsilon())
    {
        throw std::runtime_error("solveIdeCheb: arguments have inconsistent limits");
    }

    if (std::abs(h.lowerLimit() - t0) > 10 * std::numeric_limits<Real>::epsilon())
    {
        throw std::runtime_error("solveIdeCheb: inhomogeneous part has inconsistent lower limit.");
    }

    if (std::abs(h.upperLimit() - tMax) > 10 * std::numeric_limits<Real>::epsilon())
    {
        throw std::runtime_error("solveIdeCheb: inhomogeneous part has inconsistent upper limit.");
    }

    using DynamicMatrix = Eigen::Matrix<typename GetScalarType<TKernel>::type, Eigen::Dynamic, Eigen::Dynamic>;
    using DynamicVector = Eigen::Matrix<typename GetScalarType<TKernel>::type, Eigen::Dynamic, 1>;

    //
    // Copy Chebyshev coefficients of k and h into new vectors with the same size
    // (at least of size nMinCheb). This introduces potentially unnecessary copies (FIXME).
    //
    int nCheb =
        (int)nextPowerOf2(std::max((size_t)nMinCheb, std::max(h.coefficients().size(), k.coefficients().size())));
    TKernel zeroT;
    TSolution zeroY;
    if constexpr (IsScalar<TKernel>::value == true)
    {
        zeroT = Real(0);
        zeroY = Real(0);
    }
    else
    {
        zeroT = TKernel::Zero(k.coefficients()[0].rows(), k.coefficients()[0].cols());
        zeroY = TSolution::Zero(y0.rows());
    }
    std::vector<TKernel> kCoeff(nCheb, zeroT);
    std::vector<TSolution> hCoeff(nCheb, zeroY);
    for (size_t i = 0; i < k.coefficients().size(); ++i)
    {
        kCoeff[i] = k.coefficients()[i];
    }
    for (size_t i = 0; i < h.coefficients().size(); ++i)
    {
        hCoeff[i] = h.coefficients()[i];
    }

    DynamicMatrix systemMatrix = Detail::convolutionMatrix(kCoeff, tMax - t0);
    Detail::process(systemMatrix, g, tMax - t0);

    //
    // Write inhomogeneous part and initial condition vector
    //
    DynamicVector b = DynamicVector::Zero(systemMatrix.rows());
    int y0Rows      = 1;
    if constexpr (IsMatrix<TSolution>::value == true)
    {
        y0Rows = y0.rows();
        assert(y0.cols() == 1);
        assert(hCoeff[0].rows() == y0Rows);
        assert(hCoeff[0].cols() == 1);
    }

    // Inhomogeneous part
    for (int i = 0; i < nCheb - 1; ++i)
    {
        if constexpr (IsScalar<TKernel>::value == true)
        {
            b[i] = -hCoeff[i];
        }
        else
        {
            b.segment(i * y0Rows, y0Rows) = -hCoeff[i];
        }
    }

    // Initial condition
    if constexpr (IsScalar<TKernel>::value == true)
    {
        b[systemMatrix.rows() - 1] = y0;
    }
    else
    {
        b.tail(y0Rows) = y0;
    }

    //
    // Solve for Chebyshev coefficients
    //
    Eigen::PartialPivLU<Eigen::Ref<DynamicMatrix>> lu(systemMatrix);
    DynamicVector coeff = lu.solve(b);

    Cheb<TSolution> returnValue([&](Real) -> TSolution { return zeroY; }, t0, tMax, nCheb);

    for (int i = 0; i < nCheb; ++i)
    {
        if constexpr (IsScalar<TKernel>::value == true)
        {
            returnValue.coefficients()[i] = coeff[i];
        }
        else
        {
            returnValue.coefficients()[i] = coeff.segment(i * y0Rows, y0Rows);
        }
    }

    return returnValue;
}

///
/// \ingroup IDE
///
/// @brief      Convenience version without inhomogeneous part.
///
template <MatrixOrScalarType TSolution, MatrixOrScalarType TKernel>
Cheb<TSolution> solveIdeCheb(
    const TKernel& g,
    const Cheb<TKernel>& k,
    const TSolution& y0,
    Real t0,
    Real tMax,
    int nMinCheb)
{
    Cheb<TSolution> h(
        [&](Real) -> TSolution
        {
            if constexpr (IsScalar<TSolution>::value == true)
            {
                return Real(0);
            }
            else
            {
                return TSolution::Zero(y0.rows());
            }
        },
        t0, tMax, 2);

    return solveIdeCheb(g, k, h, y0, t0, tMax, nMinCheb);
}

///
/// \ingroup IDE
///
/// @brief      Solves integro-differential equations using Chebyshev methods.
///
/// This method solves integro-differential equations of the form
/// \f[
///     \dot y(t) = g y(t) + \int_{t_0}^t ds k(t-s) y(s)
/// \f]
/// with initial value \f$y(t_0)=y_0\f$.
/// Here \f$y(t)\f$ can be (real or complex) scalar- or vector-valued (of type _TSolution_). Thus the memory kernel _k_
/// must correspondingly be also a (real or complex) scalar- or matrix-valued function (of type _TKernel_).
///
/// @param g            A constant (matrix or scalar).
/// @param k            The memory kernel.
/// @param y0           Initial value.
/// @param t0           Initial time.
/// @param tMax         Time until the solution is computed.
/// @param epsAbs       Absolute error goal.
/// @param epsRel       Relative error goal.
/// @param nMinCheb     Minimum number of Chebyshev coefficients used to construct the solution in each interval of _k_.
///                     By default the same number of Chebyshev points are used as the number of points to represent the
///                     memory kernel. After the construction of the solution in each interval the number of used coefficients
///                     are reduced such that the error goal (specified by _epsAbs_ and _epsRel_) is satisfied.
///
/// @tparam TSolution   Type of the solution.
/// @tparam TKernel     Type of the memory kernel.
///
template <MatrixOrScalarType TSolution, MatrixOrScalarType TKernel>
ChebAdaptive<TSolution> solveIdeCheb(
    const TKernel& g,
    const ChebAdaptive<TKernel>& k,
    const TSolution& y0,
    Real t0,
    Real tMax,
    Real epsAbs,
    Real epsRel,
    int nMinCheb)
{
    if constexpr (IsScalar<TSolution>::value == false)
    {
        static_assert(TSolution::ColsAtCompileTime == 1, "Solution must be a scalar or vector type.");
    }

    if (k.lowerLimit() != 0)
    {
        throw std::runtime_error("solveIdeCheb: kernel must have lower limit zero.");
    }

    if (std::abs(k.upperLimit() - (tMax - t0)) > 10 * std::numeric_limits<Real>::epsilon())
    {
        throw std::runtime_error("solveIdeCheb: arguments have inconsistent limits.");
    }

    ChebAdaptive<TSolution> returnValue;

    Real tau  = k.sections()[1];
    int nCheb = nextPowerOf2((size_t)k.numCoefficients()[0]);

    int i     = 0;
    bool exit = false;
    while (exit == false)
    {
        const Real tStart = std::fma((Real)i, tau, t0);
        const Real tEnd   = std::fma((Real)(i+1), tau, t0);
        if (tEnd - tMax > -10 * std::numeric_limits<Real>::epsilon())
        {
            exit = true;
        }

        if (i == 0)
        {
            // kernel always has lower limit 0, so we use k.getCheb(0)
            Cheb sol = solveIdeCheb(g, k.getCheb(0), y0, tStart, tEnd, nMinCheb);
            sol.chopCoefficients(epsAbs, epsRel);
            returnValue.pushBack(std::move(sol));
        }
        else
        {
            Cheb<TSolution> h(
                [&](Real t) -> TSolution
                {
                    // In cases where k has narrow features, these could be missed by the
                    // integration routine if we simply used
                    //      integrateAdaptive(k(s) * returnValue(t-s))
                    // Thus we set up integration sections based on k.
                    RealVector sections  = k.sections();
                    sections.array()    += t - t0 - i * tau;
                    int l                = 1;
                    for (; l < sections.size(); ++l)
                    {
                        if (sections[l] >= t - t0)
                        {
                            sections[l] = t - t0;
                            ++l;
                            break;
                        }
                    }
                    sections.conservativeResize(l);
                    // assert(sections.size() >= 2);
                    // assert(sections[0] == t-t0-i*tau);
                    // assert(sections[sections.size()-1] == t-t0);

                    return integrateAdaptive(
                        [&](Real s) -> TSolution { return k(s) * returnValue(t - s); }, sections, epsAbs, epsRel);
                    //return integrateAdaptive(
                    //    [&](Real s) -> TSolution { return k(t - s) * returnValue(s); }, t0, t0 + i * tau, epsAbs,
                    //    epsRel);
                },
                tStart, (exit == false) ? tEnd : tMax, nCheb);

            if (exit == false)
            {
                // kernel always has lower limit 0, so we use k.getCheb(0)
                Cheb sol = solveIdeCheb(g, k.getCheb(0), h, returnValue(tStart), tStart, tEnd, nMinCheb);
                sol.chopCoefficients(epsAbs, epsRel);
                returnValue.pushBack(std::move(sol));
            }
            else
            {
                if (std::abs(k.getCheb(0).upperLimit() - (tMax - tStart)) < 10*std::numeric_limits<Real>::epsilon())
                {
                    Cheb sol = solveIdeCheb(g, k.getCheb(0), h, returnValue(tStart), tStart, tEnd, nMinCheb);
                    sol.chopCoefficients(epsAbs, epsRel);
                    returnValue.pushBack(std::move(sol));
                }
                else
                {
                    // Build "shortened" memory kernel going from 0...(tMax - tStart)
                    // The original memory kernel is going from 0...upperLimit>(tMax - tStart)
                    assert(k.getCheb(0).upperLimit() > tMax - tStart);
                    Cheb<TKernel> kShort([&](Real s) { return k.getCheb(0)(s); }, 0, tMax - tStart, nextPowerOf2(k.getCheb(0).coefficients().size()));
                    kShort.chopCoefficients(epsAbs, epsRel);
                    Cheb sol = solveIdeCheb(g, kShort, h, returnValue(tStart), tStart, tMax, nMinCheb);
                    sol.chopCoefficients(epsAbs, epsRel);
                    returnValue.pushBack(std::move(sol));
                }
            }
        }

        ++i;
    }

    return returnValue;
}

namespace Detail
{

template <MatrixType T>
Cheb<T> solveIdeChebMatrixValued(
    const T& g,
    const Cheb<T>& k,
    const Cheb<T>& h,
    const T& y0,
    Real t0,
    Real tMax,
    int nMinCheb)
{
    // Other invariants should be
    //      k.upperLimit() == (tMax - t0)
    //      h.lowerLimit() == t0
    //      h.upperLimit() == tMax
    // But due to floating point rounding errors there
    // can practically be violations, also depending on the
    // absolute values of tMax. We don't check these explicitly.
    if (k.lowerLimit() != 0)
    {
        throw std::runtime_error("solveIdeChebMatrixValued: kernel must have lower limit zero.");
    }

    //
    // Copy Chebyshev coefficients of k and h into new vectors with the same size
    // (at least of size nMinCheb). This introduces potentially unnecessary copies (FIXME).
    //
    int nCheb =
        (int)nextPowerOf2(std::max((size_t)nMinCheb, std::max(h.coefficients().size(), k.coefficients().size())));
    T zeroT = T::Zero(y0.rows(), y0.cols());
    std::vector<T> kCoeff(nCheb, zeroT);
    std::vector<T> hCoeff(nCheb, zeroT);
    for (size_t i = 0; i < k.coefficients().size(); ++i)
    {
        kCoeff[i] = k.coefficients()[i];
    }
    for (size_t i = 0; i < h.coefficients().size(); ++i)
    {
        hCoeff[i] = h.coefficients()[i];
    }

    using ElementType          = typename GetScalarType<T>::type;
    using DynamicMatrix        = Eigen::Matrix<ElementType, Eigen::Dynamic, Eigen::Dynamic>;
    DynamicMatrix systemMatrix = Detail::convolutionMatrix(kCoeff, tMax - t0);
    Detail::process(systemMatrix, g, tMax - t0);

    //
    // Write inhomogeneous part and initial condition vector
    //
    DynamicMatrix B = DynamicMatrix::Zero(systemMatrix.rows(), y0.cols());

    // Inhomogeneous part
    for (int i = 0; i < nCheb - 1; ++i)
    {
        B.block(i * y0.rows(), 0, y0.rows(), y0.cols()) = -hCoeff[i];
    }

    // Initial condition
    B.block((nCheb - 1) * y0.rows(), 0, y0.rows(), y0.cols()) = y0;

    //
    // Solve for Chebyshev coefficients
    //
    Eigen::PartialPivLU<Eigen::Ref<DynamicMatrix>> lu(systemMatrix);
    DynamicMatrix coeff = lu.solve(B);

    Cheb<T> returnValue([&](Real) -> T { return zeroT; }, t0, tMax, nCheb);

    for (int i = 0; i < nCheb; ++i)
    {
        returnValue.coefficients()[i] = coeff.block(i * y0.rows(), 0, y0.rows(), y0.cols());
    }

    return returnValue;
}

} // namespace Detail

///
/// \ingroup IDE
///
/// @brief      Computes the propagator for certain integro-differential equations using Chebyshev methods.
///
/// Consider an integro-differential equation of the form
/// \f[
///     \dot y(t) = g y(t) + \int_{t_0}^t ds k(t-s) y(s)
/// \f]
/// with initial value \f$y(t_0)=y_0\f$, where \f$y(t)\f$ is a (real or complex) vector and the memory kernel
/// _k_ is a quadratic matrix. This method computes the propagator \f$P(t-t_0)\f$ which is defined as the
/// matrix such that
/// \f[
///     y(t) = P(t-t_0) y_0
/// \f]
/// for arbitrary initial values \f$y_0\f$.
///
/// @param g            A constant (matrix or scalar).
/// @param k            The memory kernel.
/// @param t0           Initial time.
/// @param tMax         Time until the solution is computed.
/// @param nMinCheb     Minimum number of Chebyshev coefficients used to represent the solution.
///
/// @tparam T           Type of a quadratic matrix.
///
template <MatrixType T>
Cheb<T> computePropagatorIde(const T& g, const Cheb<T>& k, Real t0, Real tMax, int nMinCheb)
{
    int r = g.rows();
    int c = g.cols();
    Cheb<T> zero([r, c](Real) -> T { return T::Zero(r, c); }, t0, tMax, 2);
    T id = T::Identity(r, c);

    return Detail::solveIdeChebMatrixValued(g, k, zero, id, t0, tMax, nMinCheb);
}

///
/// \ingroup IDE
///
/// @brief      Computes the propagator for certain integro-differential equations using Chebyshev methods.
///
/// Consider an integro-differential equation of the form
/// \f[
///     \dot y(t) = g y(t) + \int_{t_0}^t ds k(t-s) y(s)
/// \f]
/// with initial value \f$y(t_0)=y_0\f$, where \f$y(t)\f$ is a (real or complex) vector and the memory kernel
/// _k_ is a quadratic matrix. This method computes the propagator \f$P(t-t_0)\f$ which is defined as the
/// matrix such that
/// \f[
///     y(t) = P(t-t_0) y_0
/// \f]
/// for arbitrary initial values \f$y_0\f$.
///
/// @param g            A constant (matrix or scalar).
/// @param k            The memory kernel.
/// @param t0           Initial time.
/// @param tMax         Time until the solution is computed.
/// @param epsAbs       Absolute error goal.
/// @param epsRel       Relative error goal.
/// @param nMinCheb     Minimum number of Chebyshev coefficients used to construct the solution in each interval of _k_.
///                     By default the same number of Chebyshev points are used as the number of points to represent the
///                     memory kernel. After the construction of the solution in each interval the number of used coefficients
///                     are reduced such that the error goal (specified by _epsAbs_ and _epsRel_) is satisfied.
///
/// @tparam T           Type of a quadratic matrix.
///
template <MatrixType T>
ChebAdaptive<T> computePropagatorIde(
    const T& g,
    const ChebAdaptive<T>& k,
    Real t0,
    Real tMax,
    Real epsAbs,
    Real epsRel,
    int nMinCheb)
{
    static_assert(T::ColsAtCompileTime == T::RowsAtCompileTime, "Matrix must be quadratic.");

    if (k.lowerLimit() != 0)
    {
        throw std::runtime_error("propagatorIdeCheb: kernel must have lower limit zero.");
    }

    if (std::abs(k.upperLimit() - (tMax - t0)) > 10 * std::numeric_limits<Real>::epsilon())
    {
        throw std::runtime_error("propagatorIdeCheb: arguments have inconsistent limits.");
    }

    ChebAdaptive<T> returnValue;

    Real tau  = k.sections()[1];
    int nCheb = nextPowerOf2((size_t)k.numCoefficients()[0]);

    int i     = 0;
    bool exit = false;
    while (exit == false)
    {
        const Real tStart = t0 + i * tau;
        const Real tEnd   = t0 + (i + 1) * tau;
        if (tEnd - tMax > -10 * std::numeric_limits<Real>::epsilon())
        {
            exit = true;
        }

        if (i == 0)
        {
            // kernel always has lower limit 0, so we use k.getCheb(0)
            Cheb sol = computePropagatorIde(g, k.getCheb(0), tStart, tEnd, nMinCheb);
            sol.chopCoefficients(epsAbs, epsRel);
            returnValue.pushBack(std::move(sol));
        }
        else
        {
            Cheb<T> h(
                [&](Real t) -> T
                {
                    // In cases where k has narrow features, these could be missed by the
                    // integration routine if we simply used
                    //      integrateAdaptive(k(s) * returnValue(t-s))
                    // Thus we set up integration sections based on k.
                    RealVector sections  = k.sections();
                    sections.array()    += t - t0 - i * tau;
                    int l                = 1;
                    for (; l < sections.size(); ++l)
                    {
                        if (sections[l] >= t - t0)
                        {
                            sections[l] = t - t0;
                            ++l;
                            break;
                        }
                    }
                    sections.conservativeResize(l);

                    return integrateAdaptive(
                        [&](Real s) -> T { return k(s) * returnValue(t - s); }, sections, epsAbs, epsRel);
                },
                tStart, (exit == false) ? tEnd : tMax, nCheb);

            if (exit == false)
            {
                // kernel always has lower limit 0, so we use k.getCheb(0)
                Cheb sol =
                    Detail::solveIdeChebMatrixValued(g, k.getCheb(0), h, returnValue(tStart), tStart, tEnd, nMinCheb);
                sol.chopCoefficients(epsAbs, epsRel);
                returnValue.pushBack(std::move(sol));
            }
            else
            {
                if (std::abs(k.getCheb(0).upperLimit() - (tMax - tStart)) < 10*std::numeric_limits<Real>::epsilon())
                {
                    Cheb sol =
                        Detail::solveIdeChebMatrixValued(g, k.getCheb(0), h, returnValue(tStart), tStart, tEnd, nMinCheb);
                    sol.chopCoefficients(epsAbs, epsRel);
                    returnValue.pushBack(std::move(sol));
                }
                else
                {
                    // Build "shortened" memory kernel going from 0...(tMax - tStart)
                    // The original memory kernel is going from 0...upperLimit>(tMax - tStart)
                    assert(k.getCheb(0).upperLimit() > tMax - tStart);
                    Cheb kShort([&](Real s) { return k.getCheb(0)(s); }, 0, tMax - tStart, nextPowerOf2(k.getCheb(0).coefficients().size()));
                    kShort.chopCoefficients(epsAbs, epsRel);
                    Cheb sol = Detail::solveIdeChebMatrixValued(g, kShort, h, returnValue(tStart), tStart, tMax, nMinCheb);
                    sol.chopCoefficients(epsAbs, epsRel);
                    returnValue.pushBack(std::move(sol));
                }
            }
        }

        ++i;
    }

    return returnValue;
}

} // namespace SciCore

#endif // SCICORE_IDE_CHEB_H
