//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

///
/// \file   integration.h
///
/// \brief  Quadrature of functions.
///

#ifndef SCICORE_INTEGRATION_H
#define SCICORE_INTEGRATION_H

#include <concepts>
#include <limits>
#include <queue>
#include <stdexcept>
#include <type_traits>

#include "Definitions.h"
#include "extern/FunctionTraits.h"
#include "IntegrationWeights.h"
#include "Utility.h"

namespace SciCore
{

namespace Detail
{
template <typename T>
struct IntegrationErrorType
{
    using type = typename std::conditional<IsScalar<T>::value, Real, RealMatrix>::type;
};

static_assert(std::is_same_v<typename IntegrationErrorType<Real>::type, Real>, "Implementation error.");
static_assert(std::is_same_v<typename IntegrationErrorType<Complex>::type, Real>, "Implementation error.");
static_assert(std::is_same_v<typename IntegrationErrorType<RealMatrix>::type, RealMatrix>, "Implementation error.");
static_assert(std::is_same_v<typename IntegrationErrorType<Matrix>::type, RealMatrix>, "Implementation error.");
} // namespace Detail

///
///  \defgroup Integration Integration
///
///  \brief Contains methods to integrate scalar and matrix-valued functions.
///
///  The module contains methods to integrate scalar and matrix-valued functions.
///
///  \{
///

template <typename T>
T integrateAdaptiveWorker_15(const T* fVals, typename SciCore::Detail::IntegrationErrorType<T>::type* absErr)
{
    using namespace SciCore::Detail;
    constexpr int numWeights = 15;
    static_assert(
        sizeof(integrateAdaptiveWorker_15_nullRule) / sizeof(Real) == numWeights,
        "Size of Gauss and Null rule doesn't match.");

    T value = Detail::gaussLegendre_1D_14[1] * fVals[0];
    T error = integrateAdaptiveWorker_15_nullRule[0] * fVals[0];
    for (int i = 1; i < numWeights; ++i)
    {
        value += Detail::gaussLegendre_1D_14[2 * i + 1] * fVals[i];
        error += integrateAdaptiveWorker_15_nullRule[i] * fVals[i];
    }

    *absErr = cwiseAbs(error);

    return value;
}

//
// \brief      Computes the integral \f$\int_a^b f(x) dx\f$ where \f$f(x)\f$ is a real or complex scalar or a matrix.
//
// \copybrief integrateAdaptive()
// Instead of specifying only the lower and upper bound of the integral, here one can also specify the initial
// subdivision of the interval in specific _sections_. The algorithm might have to further bisect each subinterval in
// order to reach the requested error goal. The final subdivions used are written to _sections_.
//
// \param      f          The function that is integrated.
// \param      sections   Initial subdivision used by the integration routine is _[sections[0],sections[1]]_,
//                        _[sections[1],sections[2]]_, ..., _[sections[n-2],sections[n-1]]_.
//                        Must be strictly increasing and contain at least two elements.
// \param      epsAbs     The absolute error goal.
// \param      epsRel     The relative error goal.
// \param      hMin       Minimum allowed interval size.
//
// \return     Returns an instance of the type \link ResultWithError ResultWithError<T> \endlink where _T_ is the
// return type of _f_.
//
template <typename FunctionT>
    requires std::invocable<FunctionT, Real>
auto integrateAdaptive(FunctionT&& f, RealVector& sections, Real epsAbs, Real epsRel, Real hMin = 1e-4)
{
    using ReturnType = typename std::invoke_result_t<FunctionT, Real>;
    using ErrorType  = typename SciCore::Detail::IntegrationErrorType<ReturnType>::type;

    struct PartialResult
    {
        ReturnType value;
        ErrorType absError;
        bool relativeErrorOk;

        Real a;
        Real b;

        bool operator<(const PartialResult& rhs) const noexcept
        {
            if (relativeErrorOk == rhs.relativeErrorOk)
            {
                return maxNorm(absError) < maxNorm(rhs.absError);
            }
            else if (relativeErrorOk == false && rhs.relativeErrorOk == true)
            {
                return false;
            }
            else // relativeErrorOk == true && rhs.relativeErrorOk == false
            {
                return true;
            }
        }
    };

    if (sections.size() < 2)
    {
        throw std::runtime_error("integrateAdaptive: needs at least two sections");
    }

    constexpr int numWeights = 15;
    std::array<ReturnType, numWeights> fVals;

    //
    // Add initial invervals from sections
    //

    // Initial interval sections[0]...sections[1]
    Real a   = sections[0];
    Real b   = sections[1];
    Real bma = (b - a) / 2;
    Real bpa = (b + a) / 2;
    for (int i = 0; i < numWeights; ++i)
    {
        Real x   = Detail::gaussLegendre_1D_14[2 * i];
        fVals[i] = bma * f(std::fma(x, bma, bpa));
    }

    ErrorType totalErr;
    ReturnType totalValue = integrateAdaptiveWorker_15(fVals.data(), &totalErr);
    Real maxRelErr        = maxNorm(cwiseQuotient(totalErr, ErrorType(cwiseAbs(totalValue)), 0.0));

    std::priority_queue<PartialResult> results;
    results.emplace(totalValue, totalErr, (maxRelErr < epsRel), a, b);

    // All other intervals: sections[1]...sections[2], sections[2]...sections[3], ....
    for (int k = 1; k < sections.size() - 1; ++k)
    {
        a   = sections[k];
        b   = sections[k + 1];
        bma = (b - a) / 2;
        bpa = (b + a) / 2;
        for (int i = 0; i < numWeights; ++i)
        {
            Real x   = Detail::gaussLegendre_1D_14[2 * i];
            fVals[i] = bma * f(std::fma(x, bma, bpa));
        }

        ErrorType err;
        ReturnType value = integrateAdaptiveWorker_15(fVals.data(), &err);
        maxRelErr        = maxNorm(cwiseQuotient(err, ErrorType(cwiseAbs(value)), 0));

        totalValue += value;
        totalErr   += err;

        results.emplace(std::move(value), std::move(err), (maxRelErr < epsRel), a, b);
    }

    // Start bisecting
    ReturnType sum12;
    ErrorType err12, relErr1, relErr2, relErr12, sumErr;
    while (maxNorm(totalErr) > epsAbs && results.top().relativeErrorOk == false &&
           maxNorm(cwiseQuotient(totalErr, ErrorType(cwiseAbs(totalValue)), 0.0)) > epsRel)
    {
        // Compute
        a           = results.top().a;
        b           = results.top().b;
        Real middle = (a + b) / 2;

        if ((middle - a < hMin) || (b - middle < hMin))
        {
            throw std::runtime_error(
                "integrateAdaptive: minimum allowed interval size reached without achieving convergence");
        }

        bma = (middle - a) / 2;
        bpa = (middle + a) / 2;
        for (int i = 0; i < numWeights; ++i)
        {
            Real x   = Detail::gaussLegendre_1D_14[2 * i];
            fVals[i] = bma * f(std::fma(x, bma, bpa));
        }
        ErrorType err1;
        ReturnType val1 = integrateAdaptiveWorker_15(fVals.data(), &err1);

        bma = (b - middle) / 2;
        bpa = (b + middle) / 2;
        for (int i = 0; i < numWeights; ++i)
        {
            Real x   = Detail::gaussLegendre_1D_14[2 * i];
            fVals[i] = bma * f(std::fma(x, bma, bpa));
        }
        ErrorType err2;
        ReturnType val2 = integrateAdaptiveWorker_15(fVals.data(), &err2);

        totalValue += val1 + val2 - results.top().value;

        relErr1 = cwiseQuotient(err1, ErrorType(cwiseAbs(totalValue)), Real(0));
        relErr2 = cwiseQuotient(err2, ErrorType(cwiseAbs(totalValue)), Real(0));

        // Update error estimates
        sum12  = val1 + val2;                           // New estimate for area 1 + 2 together
        err12  = cwiseAbs(results.top().value - sum12); // First error estimate for area 1 + 2 together
        sumErr = err1 + err2;                           // Second error estimate for area 1 + 2 together
        relErr12 =
            cwiseQuotient(err12, ErrorType(cwiseAbs(totalValue)), 0); // Relative error estimate for area 1 + 2 together
        if constexpr (IsScalar<ReturnType>::value == true)
        {
            if (err12 < sumErr)
            {
                err1 = err1 / sumErr * err12;
                err2 = err2 / sumErr * err12;
            }
            if (relErr12 < relErr1)
            {
                relErr1 = relErr12;
            }
            if (relErr12 < relErr2)
            {
                relErr2 = relErr12;
            }
        }
        else
        {
            for (int j = 0; j < err12.cols(); ++j)
            {
                for (int i = 0; i < err12.rows(); ++i)
                {
                    if (err12(i, j) < sumErr(i, j))
                    {
                        err1(i, j) = err1(i, j) / sumErr(i, j) * err12(i, j);
                        err2(i, j) = err2(i, j) / sumErr(i, j) * err12(i, j);
                    }
                    if (relErr12(i, j) < relErr1(i, j))
                    {
                        relErr1(i, j) = relErr12(i, j);
                    }
                    if (relErr12(i, j) < relErr2(i, j))
                    {
                        relErr2(i, j) = relErr12(i, j);
                    }
                }
            }
        }

        totalErr += err1 + err2 - results.top().absError;

        results.pop();

        results.emplace(std::move(val1), std::move(err1), (maxNorm(relErr1) < epsRel), a, middle);
        results.emplace(std::move(val2), std::move(err2), (maxNorm(relErr2) < epsRel), middle, b);
    }

    std::priority_queue<Real, std::vector<Real>, std::greater<Real>> intervals;
    while (results.empty() == false)
    {
        intervals.push(results.top().a);
        results.pop();
    }
    intervals.push(sections[sections.size() - 1]);
    RealVector intervalsVec(intervals.size());
    for (int i = 0; i < intervalsVec.size(); ++i)
    {
        intervalsVec[i] = intervals.top();
        intervals.pop();
    }
    sections = std::move(intervalsVec);

    return ResultWithError<ReturnType>{std::move(totalValue), maxNorm(totalErr)};
}

///
/// \brief      Computes the integral \f$\int_a^b f(x) dx\f$ where \f$f(x)\f$ is a real or complex scalar or a matrix.
///
/// \copybrief integrateAdaptive()
/// The algorithm is globally adaptive using Gauss quadrature. It works best if \f$f\f$ is smooth, but can also be
/// applied to functions having jumps or endpoint singularities (as long as they are integrable). The routine
/// successfully integrated all the functions in the Gander and Gautschi "battery test" with a relative error goal of
/// _epsRel=1e-12_ (Gander, W. and Gautschi, W. 1998. Adaptive quadrature - revisited. Tech. Rep. 306, Department of
/// Computer Science, ETH Zurich, Switzerland.)
///
/// **Example:** We want to compute the oscillatory integral
/// \f[
///     I = \int_0^{10} \frac{\sin(50x)^2}{(50x)^2} dx.
/// \f]
/// \snippet tests/IntegrationTest.cpp integrate_adaptive_example
/// This prints:
/// \code{.sh}
/// Integral is 0.0313959100209013 +- 9.08053695007558e-11 (computed with 2025 function evalations).
/// \endcode
/// The error estimate is very pessimistic in this case, because the true error is actually of the order \f$10^{-17}\f$.
///
/// \param      f          The function that is integrated.
/// \param      a          Lower integration bound.
/// \param      b          Upper integration bound.
/// \param      epsAbs     The absolute error goal.
/// \param      epsRel     The relative error goal.
/// \param      hMin       Minimum allowed interval size.
///
/// \return     Returns an instance of the type \link ResultWithError ResultWithError<T> \endlink where _T_ is the
/// return type of _f_.
///
template <typename FunctionT>
    requires std::invocable<FunctionT, Real>
auto integrateAdaptive(FunctionT&& f, Real a, Real b, Real epsAbs, Real epsRel, Real hMin = 1e-4)
{
    RealVector sections{
        {a, b}
    };
    return integrateAdaptive(std::forward<FunctionT>(f), sections, epsAbs, epsRel, hMin);
}

enum SingularityType
{
    Mild,
    Strong
};

template <typename FunctionT>
    requires std::invocable<FunctionT, Real, Real>
auto integrateDE(FunctionT&& f, Real epsAbs, Real epsRel, SingularityType s)
{
    using ReturnType = typename std::invoke_result_t<FunctionT, Real, Real>;

    int maxOrder = (s == SingularityType::Mild) ? sizeof(Detail::doubleExponential_1D_4) / sizeof(Real*)
                                                : sizeof(Detail::doubleExponential_1D_5) / sizeof(Real*);
    const auto* doubleExponential =
        (s == SingularityType::Mild) ? Detail::doubleExponential_1D_4 : Detail::doubleExponential_1D_5;
    Real prefactor = (s == SingularityType::Mild) ? 4.0 : 5.0;

    ReturnType result = prefactor * f(0.0, 1.0);

    if (SciCore::isFinite(result) == false)
    {
        throw std::runtime_error("integrateDE: result is not finite");
    }

    ReturnType oldResult = result;
    ReturnType sum;

    Real errAbs = std::numeric_limits<Real>::max();
    for (int i = 0; i < maxOrder; ++i)
    {
        int numPoints               = std::pow(2, i);
        Real delta                  = 1.0 / (2.0 * numPoints);
        const Real* nodesAndWeights = doubleExponential[i];

        sum = (prefactor * delta * nodesAndWeights[2]) *
              (f(nodesAndWeights[0], nodesAndWeights[1]) + f(-nodesAndWeights[0], nodesAndWeights[1]));

        for (int j = 1; j < numPoints; ++j)
        {
            Real x              = nodesAndWeights[3 * j + 0];
            Real distanceBorder = nodesAndWeights[3 * j + 1];
            Real w              = nodesAndWeights[3 * j + 2];
            sum                 += (prefactor * delta * w) * (f(x, distanceBorder) + f(-x, distanceBorder));
        }

        result /= 2.0;
        result += sum;

        if (SciCore::isFinite(result) == false)
        {
            throw std::runtime_error("integrateDE: result is not finite");
        }

        errAbs      = maxNorm(oldResult - result);
        Real errRel = relError(oldResult, result);
        if (errAbs < epsAbs || errRel < epsRel)
        {
            return ResultWithError<ReturnType>{std::move(result), errAbs};
        }

        oldResult = result;
    }

    return ResultWithError<ReturnType>{std::move(result), errAbs};
}

///
/// \brief      Computes the integral \f$\int_a^b f(x) dx\f$ using double exponential integration.
///
/// Computes the integral \f$\int_a^b f(x) dx\f$ using double exponential integration.
/// This method is rather insensitive to singular endpoint behavior. Use _SingularityType::Mild_ for mild endpoint
/// singularities like \f$\log{x}\f$ (or functions without endpoint singularities) and _SingularityType::Strong_ for
/// stronger ones like  \f$1/\sqrt{x}\f$. The function \f$f(x)\f$ can return either scalars or matrices (dense or
/// sparse).
///
/// The function \f$f\f$ should be coded as a function of two parameters, \f$f=f(x,\delta)\f$. The parameter
/// \f$\delta=\min\{x-a,b-x\}\geq0\f$ will be given by the distance to the closest endpoint in a numerically stable way.
///
/// **Example:** The integral \f$\int_{-2}^4 \frac{1}{x+2} dx\f$ can be computed as:
/// \snippet tests/IntegrationTest.cpp doubleExponential_simple_integral
///
/// \param      f           The function we want to integrate.
/// \param      a           Lower integration bound.
/// \param      b           Upper integration bound.
/// \param      epsAbs      Absolute error goal.
/// \param      epsRel      Relative error goal.
/// \param      s           Type of singularity.
///
template <typename FunctionT>
auto integrateDE(FunctionT&& f, Real a, Real b, Real epsAbs, Real epsRel, SingularityType s = SingularityType::Mild)
{
    using FuncTraits = StdExt::FunctionTraits<FunctionT>;
    using ReturnType = typename FuncTraits::ReturnType;

    Real b_minus_a = b - a;
    Real b_plus_a  = b + a;

    auto [returnValue, absErr] = integrateDE(
        [&f, b_minus_a, b_plus_a](Real x, Real delta) -> ReturnType
        {
            if constexpr (FuncTraits::ArgCount == 2)
            {
                return (b_minus_a / 2)*f(b_minus_a / 2 * x + b_plus_a / 2, b_minus_a / 2 * delta);
            }
            else if constexpr (FuncTraits::ArgCount == 1)
            {
                return (b_minus_a / 2)*f(b_minus_a / 2 * x + b_plus_a / 2);
            }
            else
            {
                static_assert(
                    FuncTraits::ArgCount == 1 || FuncTraits::ArgCount == 2, "Invalid number of function arguments.");
            }
        },
        epsAbs, epsRel, s);

    return ResultWithError<ReturnType>{std::move(returnValue), absErr};
}

///
/// \brief      Computes the integral \f$\int_{0}^{\infty} f(x) dx\f$ using double exponential integration.
///
/// \copybrief integrate0ToInfDE()
/// If \f$f(x)\f$ is exponentially decaying, then \ref integrate0ToInfDEDecaying() should be used for maximum
/// efficiency. The routine handles (integrable) singularities at the origin well. Slowly decaying, oscillating
/// functions are however problematic.
///
/// \see      integrateDE()
///
template <typename FunctionT>
    requires std::invocable<FunctionT, Real>
auto integrate0ToInfDE(
    FunctionT&& f,
    Real epsAbs,
    Real epsRel,
    SingularityType s    = SingularityType::Mild,
    size_t maxIterations = 20)
{
    using ReturnType     = typename std::invoke_result_t<FunctionT, Real>;
    auto transformedFunc = [&f](Real t) -> ReturnType
    {
        using std::exp;
        using std::sinh;
        using std::cosh;

        Real x = exp(std::numbers::pi_v<Real> * sinh(t));
        if constexpr (IsStdVector_v<ReturnType> == false)
        {
            return f(x) * std::numbers::pi_v<Real> * cosh(t) * x;
        }
        else // ReturnType is std::vector<Type>
        {
            ReturnType f_x = f(x);
            for (size_t i = 0; i < f_x.size(); ++i)
            {
                f_x[i] *= std::numbers::pi_v<Real> * cosh(t) * x;
            }
            return f_x;
        }
    };

    Real hMax = (s == SingularityType::Mild) ? 4.0 : 5.0;
    return integrateTrapezoid(transformedFunc, -hMax, hMax, epsAbs, epsRel, maxIterations, 5);
}

///
/// \brief      Computes the integral \f$\int_{0}^{\infty} f(x) dx\f$ using double exponential integration optimized for
/// exponentially decaying functions.
///
/// \copybrief integrate0ToInfDEDecaying()
/// Use this method if \f$f(x)\f$ is exponentially decaying for maximum efficiency, otherwise \ref integrate0ToInfDE()
/// should be used. The routine handles (integrable) singularities at the origin well.
///
/// \see      integrate0ToInfDE()
///
template <typename FunctionT>
    requires std::invocable<FunctionT, Real>
auto integrate0ToInfDEDecaying(
    FunctionT&& f,
    Real epsAbs,
    Real epsRel,
    SingularityType s    = SingularityType::Mild,
    size_t maxIterations = 20)
{
    using ReturnType     = typename std::invoke_result_t<FunctionT, Real>;
    auto transformedFunc = [&f](Real t) -> ReturnType
    {
        using std::exp;
        using std::sinh;
        using std::cosh;

        Real x    = exp(t - exp(-t));
        Real dxdt = x * (1 + exp(-t));

        if constexpr (IsStdVector_v<ReturnType> == false)
        {
            return f(x) * dxdt;
        }
        else // ReturnType is std::vector<Type>
        {
            ReturnType f_x = f(x);
            for (size_t i = 0; i < f_x.size(); ++i)
            {
                f_x[i] *= dxdt;
            }
            return f_x;
        }
    };

    Real hMax = (s == SingularityType::Mild) ? 4 : 5;
    return integrateTrapezoid(transformedFunc, -hMax, hMax, epsAbs, epsRel, maxIterations, 5);
}

template <typename FunctionT>
    requires std::invocable<FunctionT, Real>
auto integrateToInfDEDecaying(
    FunctionT&& f,
    Real a,
    Real epsAbs,
    Real epsRel,
    SingularityType s    = SingularityType::Mild,
    size_t maxIterations = 20)
{
    auto g = [&f, a](Real x)
    {
        return f(x + a);
    };

    return integrate0ToInfDEDecaying(g, epsAbs, epsRel, s, maxIterations);
}

template <class ResultT>
struct TrapezoidIntegrator
{
    Real a;
    Real b;
    ResultT lastResult;
    int n;

    TrapezoidIntegrator(Real a, Real b) noexcept : a(a), b(b), n(0)
    {
    }

    template <typename FunctionT>
    ResultT next(FunctionT&& func)
    {
        ++n;
        if (n == 1)
        {
            if constexpr (IsStdVector_v<ResultT> == false)
            {
                lastResult = (Real(0.5) * (b - a)) * (func(a) + func(b));
            }
            else // ResultT is std::vector<Type>
            {
                lastResult  = func(a);
                auto func_b = func(b);
                for (size_t i = 0; i < lastResult.size(); ++i)
                {
                    lastResult[i] += func_b[i];
                    lastResult[i] *= Real(0.5) * (b - a);
                }
            }

            return lastResult;
        }
        else
        {
            int it = 1;
            for (int j = 1; j < n - 1; ++j)
            {
                it <<= 1;
            }

            Real tnm   = it;
            Real delta = (b - a) / tnm;
            Real x     = std::fma(Real(0.5), delta, a);

            ResultT sum = func(x);
            x           += delta;

            for (int j = 1; j < it; ++j)
            {
                if constexpr (IsStdVector_v<ResultT> == false)
                {
                    sum += func(x);
                }
                else // ResultT is std::vector<Type>
                {
                    ResultT func_x = func(x);
                    for (size_t i = 0; i < sum.size(); ++i)
                    {
                        sum[i] += func_x[i];
                    }
                }

                x += delta;
            }

            if constexpr (IsStdVector_v<ResultT> == false)
            {
                lastResult = Real(0.5) * (lastResult + (b - a) * sum / tnm);
            }
            else // ResultT is std::vector<Type>
            {
                for (size_t i = 0; i < lastResult.size(); ++i)
                {
                    lastResult[i] = Real(0.5) * (lastResult[i] + (b - a) * sum[i] / tnm);
                }
            }

            return lastResult;
        }
    }
};

///
/// \brief      Integrates \f$f(x)\f$ using the [Trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule).
///
/// \copybrief integrateTrapezoid()
/// This is a rather simple quadrature algorithm, which however works extremely well for rapidly decaying functions.
/// \f$f(x)\f$ can be real or complex, and either scalar or a matrix type.
///
/// \param      f           The function we want to integrate.
/// \param      a           Lower integration bound.
/// \param      b           Upper integration bound.
/// \param      epsAbs      Absolute error goal.
/// \param      epsRel      Relative error goal.
/// \param      maxIterations  The maximum number of intervall bisections.
/// \param      minIterations  The number of intervall bisections to start with.
///
template <typename FunctionT>
    requires std::invocable<FunctionT, Real>
auto integrateTrapezoid(
    FunctionT&& f,
    Real a,
    Real b,
    Real epsAbs          = 1e-6,
    Real epsRel          = 1e-6,
    size_t maxIterations = 20,
    size_t minIterations = 5)
{
    using ResultT = typename std::invoke_result_t<FunctionT, Real>;

    assert(b > a);
    assert(minIterations > 0);
    assert(maxIterations >= minIterations);

    TrapezoidIntegrator<ResultT> trIntegrator(a, b);

    auto computeAbsError = [](const ResultT& x, const ResultT& y) -> Real
    {
        if constexpr (IsStdVector_v<ResultT> == false)
        {
            return maxNorm(x - y);
        }
        else // ResultT is std::vector<Type>
        {
            Real maxErr = 0;
            for (size_t i = 0; i < x.size(); ++i)
            {
                Real error = maxNorm(x[i] - y[i]);
                if (error > maxErr)
                {
                    maxErr = error;
                }
            }
            return maxErr;
        }
    };

    //
    // Perform the minimum number of iterations
    //

    trIntegrator.next(f);
    if (isFinite(trIntegrator.lastResult) == false)
    {
        throw std::domain_error("integrateTrapezoid: result is not finite");
    }

    ResultT oldResult;
    for (size_t i = 0; i < minIterations; ++i)
    {
        oldResult = trIntegrator.lastResult;
        trIntegrator.next(f);
        if (isFinite(trIntegrator.lastResult) == false)
        {
            throw std::domain_error("integrateTrapezoid: result is not finite");
        }
    }

    Real newError = computeAbsError(trIntegrator.lastResult, oldResult);
    if (newError < epsAbs || newError < maxNorm(trIntegrator.lastResult) * epsRel)
    {
        return ResultWithError<ResultT>{std::move(trIntegrator.lastResult), newError};
    }

    //
    // Perform all the iterations between the minimum and the maximum
    //

    for (size_t i = 0; i < maxIterations - minIterations; ++i)
    {
        oldResult = trIntegrator.lastResult;
        trIntegrator.next(f);

        if (isFinite(trIntegrator.lastResult) == false)
        {
            throw std::domain_error("integrateTrapezoid: result is not finite");
        }

        newError = computeAbsError(trIntegrator.lastResult, oldResult);
        if (newError < epsAbs || newError < maxNorm(trIntegrator.lastResult) * epsRel)
        {
            return ResultWithError<ResultT>{std::move(trIntegrator.lastResult), newError};
        }
    }

    return ResultWithError<ResultT>{std::move(trIntegrator.lastResult), newError};
}

// Integrates f(x,y) with x,y \in [-1,1] with a rule exact for all polynomials of degree <= 21.
// Error estimate by comparison of embedded rule exact for all polynomials of degree <= 11.
template <typename FunctionT>
    requires std::invocable<FunctionT, Real, Real>
auto integrateRectangleAdaptiveWorker_21(
    FunctionT&& f,
    typename SciCore::Detail::IntegrationErrorType<typename std::invoke_result_t<FunctionT, Real, Real>>::type* error)
{
    using T = typename std::invoke_result<FunctionT, Real, Real>::type;

    constexpr int numWeights         = sizeof(Detail::rectangle_21_weights) / sizeof(Real);
    constexpr int numWeightsEmbedded = sizeof(Detail::rectangle_21_weights_embedded_11) / sizeof(Real);

    static_assert(numWeights == 21, "Logic error in integration weights.");
    static_assert(numWeights + 20 == numWeightsEmbedded, "Logic error in integration weights.");

    Real x = Detail::rectangle_21_xy[0][0];
    Real y = Detail::rectangle_21_xy[0][1];

    T f0     = f(x, y);
    T value  = Detail::rectangle_21_weights[0] * f0;
    T value2 = Detail::rectangle_21_weights_embedded_11[0] * f0;

    static_assert(Detail::rectangle_21_weights_embedded_11[1] == 0);
    static_assert(Detail::rectangle_21_weights_embedded_11[1 + 20] == 0);
    x     = Detail::rectangle_21_xy[1][0];
    y     = Detail::rectangle_21_xy[1][1];
    f0    = f(x, y);
    T f1  = f(-y, x);
    T f2  = f(y, -x);
    T f3  = f(-x, -y);
    value += Detail::rectangle_21_weights[1] * (f0 + f1 + f2 + f3);

    for (int i = 2; i < numWeights; ++i)
    {
        x      = Detail::rectangle_21_xy[i][0];
        y      = Detail::rectangle_21_xy[i][1];
        f0     = f(x, y);
        f1     = f(-y, x);
        f2     = f(y, -x);
        f3     = f(-x, -y);
        value  += Detail::rectangle_21_weights[i] * (f0 + f1 + f2 + f3);
        value2 += Detail::rectangle_21_weights_embedded_11[i] * (f0 + f3) +
                  Detail::rectangle_21_weights_embedded_11[i + 20] * (f1 + f2);
    }

    *error = cwiseAbs(value - value2);

    return value;
}

template <typename FunctionT>
    requires std::invocable<FunctionT, Real, Real>
auto integrateRectangleAdaptiveWorker_21(
    FunctionT&& f,
    StaticRealVector<2> lower,
    StaticRealVector<2> upper,
    typename SciCore::Detail::IntegrationErrorType<typename std::invoke_result_t<FunctionT, Real, Real>>::type* error)
{
    using T = typename std::invoke_result<FunctionT, Real, Real>::type;

    Real detJ = std::abs((upper[0] - lower[0]) * (upper[1] - lower[1])) / 4;
    Real a0   = (upper[0] + lower[0]) / 2;
    Real b0   = (upper[0] - lower[0]) / 2;
    Real a1   = (upper[1] + lower[1]) / 2;
    Real b1   = (upper[1] - lower[1]) / 2;

    return integrateRectangleAdaptiveWorker_21(
        [&f, detJ, a0, b0, a1, b1](Real x, Real y) -> T { return detJ * f(a0 + b0 * x, a1 + b1 * y); }, error);
}

///
/// \brief      Computes the integral \f$ \int_{\text{lower}_0}^{\text{upper}_0} dx
///             \int_{\text{lower}_1}^{\text{upper}_1} dy f(x,y) \f$ using adaptive bisection.
///
/// \copybrief integrate2DAdaptive()
/// The used algorithm is optimized for smooth integrands. This means that the error estimates for non-smooth
/// functions are likely to be unreliable. Non-analycities should always be moved to the boundary of the
/// integration region.
///
/// **Example:** We want to compute the Gaussian integral
/// \f[
///     I = \int_{-3}^{3} dx \int_{-3}^{3} dy \, e^{-25 \left(x - \frac{1}{2} \right)^2 - 25 \left(y - \frac{1}{2}
///     \right)^2}.
/// \f]
/// \snippet tests/IntegrationTest.cpp integrate2DAdaptive_example
/// This prints:
/// \code{.sh}
/// Integral = 0.125663706143589 +- 2.91341992800507e-13
/// Computed with 9153 function evalations.
/// \endcode
/// The true absolute error is rougly \f$2.5 \cdot 10^{-15}\f$ (and the true relative error is \f$2 \cdot 10^{-14}\f$).
/// Note that the number of evaluations may seem high, but compares favourably to many commonly used quadrature algorithms,
/// such as the standard one used in _SciPy_: Computing the same integral using _dblquad_ from _SciPy_ we get
/// \code{.py}
/// import numpy as np
/// from scipy import integrate
///
/// count = 0
/// def f(x, y):
///     global count
///     count += 1
///     return np.exp(-25*pow(x - 0.5, 2) - 25*pow(y - 0.5, 2))
///
/// I = integrate.dblquad(f, -3, 3, -3, 3, epsrel = 1e-12, epsabs = 0)
/// print('Integral = ', I[0], ' +- ', I[1])
/// print('Computed with ', count, 'function evalations.')
/// \endcode
/// which prints
/// \code{.sh}
/// Integral =  0.12566370614359174  +-  4.553843017609735e-15
/// Computed with  159201 function evalations.
/// \endcode
///
/// \param      f               The function that is integrated.
/// \param      lower           Lower integral boundaries.
/// \param      upper           Upper integral boundaries.
/// \param      epsAbs          The absolute error goal.
/// \param      epsRel          The relative error goal.
///
/// \return     Returns an instance of the type \link ResultWithError ResultWithError<T> \endlink where _T_ is the return type of _f_.
///
template <typename FunctionT>
    requires std::invocable<FunctionT, Real, Real>
auto integrate2DAdaptive(FunctionT&& f, StaticRealVector<2> lower, StaticRealVector<2> upper, Real epsAbs, Real epsRel)
{
    using ReturnType = typename std::invoke_result_t<FunctionT, Real, Real>;
    using ErrorType  = typename SciCore::Detail::IntegrationErrorType<ReturnType>::type;

    struct PartialResult
    {
        ReturnType value;
        ErrorType absError;
        bool relativeErrorOk;

        StaticRealVector<2> lower;
        StaticRealVector<2> upper;

        bool operator<(const PartialResult& rhs) const noexcept
        {
            if (relativeErrorOk == rhs.relativeErrorOk)
            {
                return maxNorm(absError) < maxNorm(rhs.absError);
            }
            else if (relativeErrorOk == false && rhs.relativeErrorOk == true)
            {
                return false;
            }
            else // relativeErrorOk == true && rhs.relativeErrorOk == false
            {
                return true;
            }
        }
    };

    std::priority_queue<PartialResult> results;

    //
    // Compute initial approximation
    //
    ErrorType totalErr;
    ReturnType totalValue = integrateRectangleAdaptiveWorker_21(std::forward<FunctionT>(f), lower, upper, &totalErr);
    Real maxRelErr        = maxNorm(cwiseQuotient(totalErr, ErrorType(cwiseAbs(totalValue)), 0.0));
    results.emplace(totalValue, totalErr, (maxRelErr < epsRel), lower, upper);

    //
    // Cut rectangle
    //
    ReturnType sum12;
    ErrorType err12, relErr1, relErr2, relErr12, sumErr;
    while (maxNorm(totalErr) > epsAbs && results.top().relativeErrorOk == false &&
           maxNorm(cwiseQuotient(totalErr, ErrorType(cwiseAbs(totalValue)), 0.0)) > epsRel)
    {
        lower = results.top().lower;
        upper = results.top().upper;

        StaticRealVector<2> lower1, upper1, lower2, upper2;
        if ((upper[1] - lower[1]) > (upper[0] - lower[0])) // split horizontally
        {
            lower1 = lower;
            upper1 = StaticRealVector<2>{upper[0], (upper[1] + lower[1]) / 2};
            lower2 = StaticRealVector<2>{lower[0], (upper[1] + lower[1]) / 2};
            upper2 = upper;
        }
        else // split vertically
        {
            lower1 = lower;
            upper1 = StaticRealVector<2>{(upper[0] + lower[0]) / 2, upper[1]};
            lower2 = StaticRealVector<2>{(upper[0] + lower[0]) / 2, lower[1]};
            upper2 = upper;
        }

        ErrorType err1, err2;
        ReturnType val1 = integrateRectangleAdaptiveWorker_21(std::forward<FunctionT>(f), lower1, upper1, &err1);
        ReturnType val2 = integrateRectangleAdaptiveWorker_21(std::forward<FunctionT>(f), lower2, upper2, &err2);

        totalValue += val1 + val2 - results.top().value;

        relErr1 = cwiseQuotient(err1, ErrorType(cwiseAbs(totalValue)), 0.0);
        relErr2 = cwiseQuotient(err2, ErrorType(cwiseAbs(totalValue)), 0.0);

        // Update error estimates
        sum12    = val1 + val2;                           // New estimate for area 1 + 2 together
        err12    = cwiseAbs(results.top().value - sum12); // First error estimate for area 1 + 2 together
        sumErr   = err1 + err2;                           // Second error estimate for area 1 + 2 together
        relErr12 = cwiseQuotient(
            err12, ErrorType(cwiseAbs(totalValue)), 0.0); // Relative error estimate for area 1 + 2 together
        if constexpr (IsScalar<ReturnType>::value == true)
        {
            if (err12 < sumErr)
            {
                err1 = err1 / sumErr * err12;
                err2 = err2 / sumErr * err12;
            }
            if (relErr12 < relErr1)
            {
                relErr1 = relErr12;
            }
            if (relErr12 < relErr2)
            {
                relErr2 = relErr12;
            }
        }
        else
        {
            for (int j = 0; j < err12.cols(); ++j)
            {
                for (int i = 0; i < err12.rows(); ++i)
                {
                    if (err12(i, j) < sumErr(i, j))
                    {
                        err1(i, j) = err1(i, j) / sumErr(i, j) * err12(i, j);
                        err2(i, j) = err2(i, j) / sumErr(i, j) * err12(i, j);
                    }
                    if (relErr12(i, j) < relErr1(i, j))
                    {
                        relErr1(i, j) = relErr12(i, j);
                    }
                    if (relErr12(i, j) < relErr2(i, j))
                    {
                        relErr2(i, j) = relErr12(i, j);
                    }
                }
            }
        }

        totalErr += err1 + err2 - results.top().absError;

        results.pop();

        results.emplace(std::move(val1), std::move(err1), (maxNorm(relErr1) < epsRel), lower1, upper1);
        results.emplace(std::move(val2), std::move(err2), (maxNorm(relErr2) < epsRel), lower2, upper2);
    }

    return ResultWithError<ReturnType>{std::move(totalValue), maxNorm(totalErr)};
}

/// \} // end of Integration

} // namespace SciCore

#endif // SCICORE_INTEGRATION_H
