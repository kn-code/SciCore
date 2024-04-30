//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

///
/// \file   ChebAdaptive.h
///
/// \brief  Adaptive Chebyshev interpolation.
///

#ifndef SCICORE_CHEB_ADAPTIVE_H
#define SCICORE_CHEB_ADAPTIVE_H

#include <map>

#include "Serialization.h"
#include <cereal/types/map.hpp>

#include "Cheb.h"

namespace SciCore
{

///
/// \ingroup Interpolation
///
/// \brief    Constructs approximations to functions adaptively by bisecting the interval until each subinterval is
/// well-represented by a Chebyshev approximation.
///
/// This class approximates functions by bisecting an interval adaptively until the function is well-represented by a
/// Chebyshev approximant on each subinterval seprarately. To illustrate this consider the function \f[
///     f(x) = \cos\big(\cos(x) + 3 \sin(x) + 2 \cos(2 x) + 3 \cos(3 x)\big) + \frac{1}{\cosh(400 (x-0.4))}
/// \f]
///
/// \image html doc-adaptive-cheb-f.png width=60%
///
/// \f$f\f$ has a narrow peak at \f$x=0.4\f$ ontop of a otherwise smooth structure, which would be difficult to capture
/// for a single Cheb approximant. Using ChebAdaptive this function is however easily approximated:
///
/// \snippet tests/ChebAdaptiveTest.cpp Basic adaptive cheb usage
///
/// The interval is automatically decomposed into the following sections:
/// \code{.cpp}
/// std::cout << "Sections:" << chebf.sections().transpose() << "\n"
///           << "Number coefficients: " << chebf.numCoefficients().transpose() << "\n";
/// \endcode
/// \code{.sh}
/// Sections 1:       0  0.21875 0.328125 0.382812 0.410156   0.4375    0.875     1.75      3.5
/// Number coefficients 1:  16  14  26 113  27  47  32  53
/// \endcode
///
/// This means that in the interval \f$[0, 0.21875]\f$ an approximant with 16 coefficients is used, in the interval
/// \f$[0.21875, 0.328125]\f$ an approximant with 14 coefficients is used and so on. With the approximation in hand one
/// can now easily compute the derivative or integral using diff() or integrate() respectively.
///
/// \headerfile ChebAdaptive.h <SciCore/ChebAdaptive.h>
///
template <MatrixOrScalarType T>
class ChebAdaptive
{
  public:
    ///
    /// \brief The return type of the interpolated function (Real, Complex, Matrix, ...)
    ///
    using ElementType = T;

    ///
    /// \brief Constructs an empty object.
    ///
    ChebAdaptive() noexcept
    {
    }

    ChebAdaptive(const ChebAdaptive& other)     = default;
    ChebAdaptive(ChebAdaptive&& other) noexcept = default;

    ChebAdaptive& operator=(const ChebAdaptive& other)     = default;
    ChebAdaptive& operator=(ChebAdaptive&& other) noexcept = default;

    ///
    /// \brief      Creates a piecewise Chebyshev interpolation of the function \f$f(x)\f$ with \f$x\in[a,b]\f$.
    ///
    /// Creates a piecewise Chebyshev interpolation of the function \f$f(x)\f$ with \f$x\in[a,b]\f$.
    /// The routine first constructs a Chebyshev approximation over the entire interval
    /// \f$[a,b]\f$ with increasing orders \f$ n_\text{Start},3\cdot n_\text{Start},...,3^{n_\text{Trip}}\cdot
    /// n_\text{Start}\f$. If the error goal is not achieved then the interval is bisected in the middle and the
    /// procedure is repeated recursively. The bisection is stopped if the error goal is met or the
    /// interval length is less than _hMin_.
    ///
    /// \param      f               The function for which the piecewise Chebyshev approximation is computed.
    /// \param      a               Lower interval point.
    /// \param      b               Upper interval point.
    /// \param      epsAbs          Absolute error goal.
    /// \param      epsRel          Relative error goal.
    /// \param      hMin            The minimum allowed interval length.
    /// \param      ok              Set to _true_ if error goal was achieved, otherwise set to _false_.
    /// \param      nStart          Order of the initial Chebyshev approximation.
    /// \param      nTrip           The maximum number of triplings.
    /// \param      buffer          A buffer containing \f$n_{\text{Start}} 3^{n_{\text{Trip}}}\f$ elements.
    ///                             If set to _nullptr_ then the buffer is allocated automatically.
    ///
    template <typename FunctionT>
        requires std::invocable<FunctionT, Real>
    ChebAdaptive(
        FunctionT&& f,
        Real a,
        Real b,
        Real epsAbs,
        Real epsRel,
        Real hMin,
        bool* ok   = nullptr,
        int nStart = 16,
        int nTrip  = 2,
        T* buffer  = nullptr)
    {
        bool errorOk;
        if (buffer == nullptr)
        {
            std::vector<T> newBuffer((size_t)(nStart * std::pow(3, nTrip)));
            errorOk =
                _addInterval(std::forward<FunctionT>(f), a, b, epsAbs, epsRel, hMin, nStart, nTrip, newBuffer.data());
        }
        else
        {
            errorOk = _addInterval(std::forward<FunctionT>(f), a, b, epsAbs, epsRel, hMin, nStart, nTrip, buffer);
        }

        if (ok != nullptr)
        {
            *ok = errorOk;
        }
    }

    template <typename FunctionT>
        requires std::invocable<FunctionT, Real>
    ChebAdaptive(
        FunctionT&& f,
        Real a,
        Real b,
        Real epsAbs,
        Real epsRel,
        Real hMin,
        tf::Executor& executor,
        bool* ok   = nullptr,
        int nStart = 16,
        int nTrip  = 2,
        T* buffer  = nullptr)
    {
        bool errorOk;
        if (buffer == nullptr)
        {
            std::vector<T> newBuffer((size_t)(nStart * std::pow(3, nTrip)));
            errorOk = _addInterval(
                std::forward<FunctionT>(f), a, b, epsAbs, epsRel, hMin, executor, nStart, nTrip, newBuffer.data());
        }
        else
        {
            errorOk =
                _addInterval(std::forward<FunctionT>(f), a, b, epsAbs, epsRel, hMin, executor, nStart, nTrip, buffer);
        }

        if (ok != nullptr)
        {
            *ok = errorOk;
        }
    }

    ///
    /// \brief      Creates a piecewise Chebyshev interpolation of the function \f$f(x)\f$ such that at least one
    /// Chebyshev interpolant is used in each of the given \a sections.
    ///
    /// Creates a piecewise Chebyshev interpolation of the function \f$f(x)\f$ such that at least one
    /// Chebyshev interpolant is used in each of the given \a sections.
    /// The routine first constructs a Chebyshev approximation in each section of \a
    /// sections with increasing orders \f$ n_\text{Start},3\cdot n_\text{Start},...,3^{n_\text{Trip}}\cdot
    /// n_\text{Start}\f$. If the error estimate is too high then the interval is bisected in the middle and the
    /// procedure is repeated recursively. The bisection is stopped if the error estimate becomes small enough or the
    /// interval length is less than _hMin_.
    ///
    /// \param      f               The function for which the piecewise Chebyshev approximation is computed.
    /// \param      sections        Example: if \a sections \f$=\{a,b,c\}\f$ then at least one Chebyshev interpolant is
    ///                             used in each of \f$[a,b]\f$ and \f$[b,c]\f$.
    /// \param      epsAbs          Absolute error goal.
    /// \param      epsRel          Relative error goal.
    /// \param      hMin            The minimum allowed interval length.
    /// \param      ok              Set to _true_ if error goal was achieved, otherwise set to _false_.
    /// \param      nStart          Order of the initial Chebyshev approximation.
    /// \param      nTrip           The maximum number of triplings. Can currently at most be equal to 2.
    /// \param      buffer          A buffer containing \f$n_{\text{Start}} 3^{n_{\text{Trip}}}\f$ elements.
    ///                             If set to _nullptr_ then the buffer is allocated automatically.
    ///
    template <typename FunctionT>
        requires std::invocable<FunctionT, Real>
    ChebAdaptive(
        FunctionT&& f,
        const RealVector& sections,
        Real epsAbs,
        Real epsRel,
        Real hMin,
        bool* ok   = nullptr,
        int nStart = 16,
        int nTrip  = 2,
        T* buffer  = nullptr)
    {
        assert(sections.size() >= 2);

        bool errorOk = true;
        if (buffer == nullptr)
        {
            std::vector<T> newBuffer((size_t)(nStart * pow(3, nTrip)));
            for (int i = 0; i < sections.size() - 1; ++i)
            {
                if (_addInterval(
                        std::forward<FunctionT>(f), sections[i], sections[i + 1], epsAbs, epsRel, hMin, nStart, nTrip,
                        newBuffer.data()) == false)
                {
                    errorOk = false;
                }
            }
        }
        else
        {
            for (int i = 0; i < sections.size() - 1; ++i)
            {
                if (_addInterval(
                        std::forward<FunctionT>(f), sections[i], sections[i + 1], epsAbs, epsRel, hMin, nStart, nTrip,
                        buffer) == false)
                {
                    errorOk = false;
                }
            }
        }

        if (ok != nullptr)
        {
            *ok = errorOk;
        }
    }

    template <typename FunctionT>
        requires std::invocable<FunctionT, Real>
    ChebAdaptive(
        FunctionT&& f,
        const RealVector& sections,
        Real epsAbs,
        Real epsRel,
        Real hMin,
        tf::Executor& executor,
        bool* ok   = nullptr,
        int nStart = 16,
        int nTrip  = 2,
        T* buffer  = nullptr)
    {
        assert(sections.size() >= 2);

        bool errorOk = true;
        if (buffer == nullptr)
        {
            std::vector<T> newBuffer((size_t)(nStart * pow(3, nTrip)));
            for (int i = 0; i < sections.size() - 1; ++i)
            {
                if (_addInterval(
                        std::forward<FunctionT>(f), sections[i], sections[i + 1], epsAbs, epsRel, hMin, executor,
                        nStart, nTrip, newBuffer.data()) == false)
                {
                    errorOk = false;
                }
            }
        }
        else
        {
            for (int i = 0; i < sections.size() - 1; ++i)
            {
                if (_addInterval(
                        std::forward<FunctionT>(f), sections[i], sections[i + 1], epsAbs, epsRel, hMin, executor,
                        nStart, nTrip, buffer) == false)
                {
                    errorOk = false;
                }
            }
        }

        if (ok != nullptr)
        {
            *ok = errorOk;
        }
    }

    ChebAdaptive(Cheb<T>* chebs, int n)
    {
        assert(n > 0);
        Real oldUpper;
        for (int i = 0; i < n; ++i)
        {
            Real lower = chebs[i].lowerLimit();
            if (i > 0)
            {
                if (lower != oldUpper)
                {
                    throw std::runtime_error("ChebAdaptive(Cheb<T>* , int ): Sections must be connected.");
                }
            }
            oldUpper = chebs[i].upperLimit();
            _chebs.emplace(lower, std::move(chebs[i]));
        }
    }

    const Cheb<T>& getCheb(Real x) const
    {
        if (_chebs.empty() == true)
        {
            throw std::runtime_error("ChebAdaptive::get(): Can't be used on empty object.");
        }

        auto it = _chebs.upper_bound(x);
        if (it == _chebs.end())
        {
            return _chebs.rbegin()->second;
        }
        else if (it == _chebs.begin())
        {
            return it->second;
        }
        else
        {
            --it;
            return it->second;
        }
    }

    ///
    /// \brief Evaluates the approximation at \f$x\f$ and returns the result.
    ///
    T operator()(Real x) const
    {
        const Cheb<T>& c = getCheb(x);
        return c(x);
    }

    ///
    /// \brief Evaluates the interpolation at at all points in \a x. Only works if \a ElementType is \a Real.
    ///
    RealVector operator()(const RealVector& x)
        requires std::is_same_v<std::decay_t<T>, Real>
    {
        RealVector returnValue(x.size());
        for (int i = 0; i < x.size(); ++i)
        {
            returnValue[i] = this->operator()(x[i]);
        }
        return returnValue;
    }

    ///
    /// \brief Evaluates the interpolation at at all points in \a x. Only works if \a ElementType is \a Complex.
    ///
    Vector operator()(const RealVector& x)
        requires std::is_same_v<std::decay_t<T>, Complex>
    {
        Vector returnValue(x.size());
        for (int i = 0; i < x.size(); ++i)
        {
            returnValue[i] = this->operator()(x[i]);
        }
        return returnValue;
    }

    ///
    /// \brief      Equality operator.
    ///
    bool operator==(const ChebAdaptive& other) const noexcept
    {
        return (_chebs == other._chebs);
    }

    ///
    /// \brief      Unequality operator.
    ///
    bool operator!=(const ChebAdaptive& other) const noexcept
    {
        return !operator==(other);
    }

    ///
    /// \brief  Computes the derivative of the represented function and returns it as a new object.
    ///
    /// \copybrief diff
    ///
    /// \warning Computing the derivative is numerically not very stable,
    ///          see also <a href="https://www.sciencedirect.com/science/article/pii/0021999192902743">On the errors
    ///          incurred calculating derivatives using Chebyshev polynomials.</a>
    ///
    ChebAdaptive diff() const
    {
        if (_chebs.empty() == true)
        {
            throw std::runtime_error("ChebAdaptive::diff(): Can't be used on empty object.");
        }

        ChebAdaptive returnValue;

        for (const auto& pair : _chebs)
        {
            returnValue._chebs.emplace(pair.first, pair.second.diff());
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
    /// is returned, where \a a corresponds to lowerLimit().
    ///
    ChebAdaptive integrate() const
    {
        if (_chebs.empty() == true)
        {
            throw std::runtime_error("ChebAdaptive::integrate(): Can't be used on empty object.");
        }

        auto it = _chebs.begin();
        ChebAdaptive returnValue;
        returnValue._chebs[it->first] = it->second.integrate();

        auto rBeg  = returnValue._chebs.rbegin();
        T constant = rBeg->second(rBeg->second.upperLimit());

        ++it;
        for (; it != _chebs.end(); ++it)
        {
            auto result                    = it->second.integrate();
            result                        += constant;
            returnValue._chebs[it->first]  = result;
            rBeg                           = returnValue._chebs.rbegin();
            constant                       = rBeg->second(rBeg->second.upperLimit());
        }

        return returnValue;
    }

    ///
    /// \brief      Returns a list of the used intervals, where each interval is represented by one _Cheb_ approximant.
    ///
    RealVector sections() const
    {
        if (_chebs.empty() == true)
        {
            throw std::runtime_error("ChebAdaptive::sections(): Can't be used on empty object.");
        }

        int i = 0;
        RealVector returnValue(_chebs.size() + 1);
        for (const auto& pair : _chebs)
        {
            returnValue[i] = pair.first;
            ++i;
        }
        returnValue[_chebs.size()] = _chebs.rbegin()->second.upperLimit();
        return returnValue;
    }

    ///
    /// \brief      Returns the number of coefficients used in each interval.
    ///
    IntVector numCoefficients() const
    {
        int i = 0;
        IntVector returnValue(_chebs.size());
        for (const auto& pair : _chebs)
        {
            returnValue[i] = (int)pair.second.coefficients().size();
            ++i;
        }
        return returnValue;
    }

    ///
    /// \brief  Returns the lower limit of the approximation.
    ///
    Real lowerLimit() const
    {
        if (_chebs.empty() == true)
        {
            throw std::runtime_error("ChebAdaptive::lowerLimit(): Can't be used on empty object.");
        }

        RealVector sec = this->sections();
        return this->sections()[0];
    }

    ///
    /// \brief  Returns the upper limit of the approximation.
    ///
    Real upperLimit() const
    {
        if (_chebs.empty() == true)
        {
            throw std::runtime_error("ChebAdaptive::upperLimit(): Can't be used on empty object.");
        }

        RealVector sec = this->sections();
        return sec[sec.size() - 1];
    }

    void pushBack(Cheb<T>&& c)
    {
        Real lowerLimit = c.lowerLimit();
        if (_chebs.empty() == true)
        {
            _chebs[lowerLimit] = std::move(c);
            return;
        }

        auto last = _chebs.rbegin();
        if (last->second.upperLimit() != lowerLimit)
        {
            if (std::abs(last->second.upperLimit() - lowerLimit) < 10 * std::numeric_limits<Real>::epsilon())
            {
                lowerLimit = last->second.upperLimit();
            }
            else
            {
                throw std::runtime_error("ChebAdaptive::pushBack(): Inconsistent limits.");
            }
        }

        if (_chebs.find(lowerLimit) != _chebs.end())
        {
            throw std::runtime_error("ChebAdaptive::pushBack(): Logic error.");
        }

        _chebs[lowerLimit] = std::move(c);
    }

    void pushBack(const Cheb<T>& c)
    {
        pushBack(Cheb<T>(c));
    }

    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(CEREAL_NVP(_chebs));
    }

  private:
    template <typename FunctionT>
    bool _addInterval(
        FunctionT&& f,
        Real a,
        Real b,
        Real epsAbs,
        Real epsRel,
        Real hMin,
        int nStart,
        int nTrip,
        T* buffer)
    {
        bool ok = false;
        Cheb<T> c(std::forward<FunctionT>(f), a, b, epsAbs, epsRel, nStart, nTrip, &ok, buffer);

        if (ok == true)
        {
            _chebs[a] = std::move(c);
        }
        else if (b - a < hMin)
        {
            _chebs[a] = std::move(c);
        }
        else
        {
            bool okLeft =
                _addInterval(std::forward<FunctionT>(f), a, (a + b) / 2, epsAbs, epsRel, hMin, nStart, nTrip, buffer);
            bool okRight =
                _addInterval(std::forward<FunctionT>(f), (a + b) / 2, b, epsAbs, epsRel, hMin, nStart, nTrip, buffer);
            ok = (okLeft == true) && (okRight == true);
        }
        return ok;
    }

    template <typename FunctionT>
    bool _addInterval(
        FunctionT&& f,
        Real a,
        Real b,
        Real epsAbs,
        Real epsRel,
        Real hMin,
        tf::Executor& executor,
        int nStart,
        int nTrip,
        T* buffer)
    {
        bool ok = false;
        Cheb<T> c(std::forward<FunctionT>(f), a, b, epsAbs, epsRel, executor, nStart, nTrip, &ok, buffer);

        if (ok == true)
        {
            _chebs[a] = std::move(c);
        }
        else if (b - a < hMin)
        {
            _chebs[a] = std::move(c);
        }
        else
        {
            bool okLeft = _addInterval(
                std::forward<FunctionT>(f), a, (a + b) / 2, epsAbs, epsRel, hMin, executor, nStart, nTrip, buffer);
            bool okRight = _addInterval(
                std::forward<FunctionT>(f), (a + b) / 2, b, epsAbs, epsRel, hMin, executor, nStart, nTrip, buffer);
            ok = (okLeft == true) && (okRight == true);
        }
        return ok;
    }

    std::map<Real, Cheb<T>> _chebs; // Key: lower limit, value: Chebyshev approximation starting at the lower limit
};

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, Real a, Real b, Real epsAbs, Real epsRel, Real hMin)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, Real a, Real b, Real epsAbs, Real epsRel, Real hMin, tf::Executor&)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, Real a, Real b, Real epsAbs, Real epsRel, Real hMin, bool* ok)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, Real a, Real b, Real epsAbs, Real epsRel, Real hMin, tf::Executor&, bool* ok)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, Real a, Real b, Real epsAbs, Real epsRel, Real hMin, bool* ok, int nStart)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, Real a, Real b, Real epsAbs, Real epsRel, Real hMin, tf::Executor&, bool* ok, int nStart)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, Real a, Real b, Real epsAbs, Real epsRel, Real hMin, bool* ok, int nStart, int nTrip)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(
    FunctionT&& f,
    Real a,
    Real b,
    Real epsAbs,
    Real epsRel,
    Real hMin,
    tf::Executor&,
    bool* ok,
    int nStart,
    int nTrip) -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, const RealVector& sections, Real epsAbs, Real epsRel, Real hMin)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, const RealVector& sections, Real epsAbs, Real epsRel, Real hMin, tf::Executor&)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, const RealVector& sections, Real epsAbs, Real epsRel, Real hMin, bool* ok)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, const RealVector& sections, Real epsAbs, Real epsRel, Real hMin, tf::Executor&, bool* ok)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, const RealVector& sections, Real epsRel, Real hMin, bool* ok, int nStart)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, const RealVector& sections, Real epsRel, Real hMin, tf::Executor&, bool* ok, int nStart)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(FunctionT&& f, const RealVector& sections, Real epsRel, Real hMin, bool* ok, int nStart, int nTrip)
    -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

template <typename FunctionT>
ChebAdaptive(
    FunctionT&& f,
    const RealVector& sections,
    Real epsRel,
    Real hMin,
    tf::Executor&,
    bool* ok,
    int nStart,
    int nTrip) -> ChebAdaptive<std::invoke_result_t<FunctionT, Real>>;

#ifndef SCICORE_DONT_PRECOMPILE_TEMPLATES
extern template class ChebAdaptive<Real>;
extern template class ChebAdaptive<Complex>;
extern template class ChebAdaptive<RealVector>;
extern template class ChebAdaptive<Vector>;
extern template class ChebAdaptive<RealMatrix>;
extern template class ChebAdaptive<Matrix>;
#endif // SCICORE_DONT_PRECOMPILE_TEMPLATES

} // namespace SciCore

#endif // SCICORE_CHEB_ADAPTIVE_H
