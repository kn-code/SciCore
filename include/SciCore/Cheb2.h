///
/// \file   Cheb2.h
///
/// \brief  Chebyshev interpolation in two dimensions.
///

#ifndef SCICORE_CHEB2_H
#define SCICORE_CHEB2_H

#include <concepts>

#include "Definitions.h"
#include "Cheb.h"
#include "DCT.h"
#include "DCT2.h"
#include "Parallel.h"
#include "Utility.h"
#include "SciCore_export.h"
#include "Serialization.h"

namespace SciCore {

///
/// \ingroup ChebyshevApproximations
///
/// \brief Computes Chebyshev interpolations for scalar- and matrix-valued functions of two real parameters.
/// 
/// This class computes Chebyshev interpolations for scalar- and matrix-valued functions of two real parameters.
/// 
/// **Example:**
/// \snippet tests/Cheb2Test.cpp Basic cheb2 usage
/// 
/// \headerfile Cheb2.h <SciCore/Cheb2.h>
///
template <MatrixOrScalarType T>
class SCICORE_EXPORT Cheb2
{
public:

    ///
    /// \brief Scalar or matrix type.
    ///
    using ElementType = T;

    ///
    /// \brief Constructs an empty object.
    ///
    Cheb2() noexcept
        : _lower{0,0}, _upper{0,0}, _n{0,0}
    {}

    ///
    /// \brief Copy constructor.
    ///
    Cheb2(const Cheb2& other) = default;

    ///
    /// \brief Move constructor.
    ///
    Cheb2(Cheb2&& other) noexcept = default;

    ///
    /// \brief Copy assignment operator.
    ///
    Cheb2& operator=(const Cheb2& other) = default;

    ///
    /// \brief Move assignment operator.
    ///
    Cheb2& operator=(Cheb2&& other) noexcept = default;

    ///
    /// \brief      Equality operator.
    ///
    bool operator==(const Cheb2& other) const noexcept
    {
        return (_c == other._c) &&
               (_lower == other._lower) &&
               (_upper == other._upper) &&
               (_n == other._n) &&
               (_evaluatedRows == other._evaluatedRows) &&
               (_evaluatedCols == other._evaluatedCols);
    }

    ///
    /// \brief      Unequality operator.
    ///
    bool operator!=(const Cheb2& other) const noexcept
    {
        return !operator==(other);
    }

    ///
    /// \brief  Constructs the Chebyshev approximation of the function \f$f(x,y)\f$.
    /// 
    /// Constructs the Chebyshev approximation of the function \f$f(x,y)\f$.
    ///
    /// \param  f          The function that is approximated.
    /// \param  lower      The lower limit of both function parameters.
    /// \param  upper      The upper limit of both function parameters.
    /// \param  n          The number of Chebyshev coefficients per dimension.
    ///
    /// \tparam FunctionT  A function type that takes two Real arguments.
    ///
    template <typename FunctionT>
        requires std::invocable<FunctionT, Real, Real>
    Cheb2(FunctionT&& f, StaticRealVector<2> lower, StaticRealVector<2> upper, StaticIntVector<2> n)
        : _lower(lower), _upper(upper), _n(n)
    {
        _compute(std::forward<FunctionT>(f));
    }

    ///
    /// \brief  Constructs the Chebyshev approximation of the function \f$f(x,y)\f$ in parallel.
    /// 
    /// Constructs the Chebyshev approximation of the function \f$f(x,y)\f$.
    ///
    /// \param  f          The function that is approximated.
    /// \param  lower      The lower limit of both function parameters.
    /// \param  upper      The upper limit of both function parameters.
    /// \param  n          The number of Chebyshev coefficients per dimension.
    /// \param  executor   Taskflow executor for parallel execution.
    ///
    /// \tparam FunctionT  A function type that takes two Real arguments.
    ///
    template <typename FunctionT>
        requires std::invocable<FunctionT, Real, Real>
    Cheb2(FunctionT&& f, StaticRealVector<2> lower, StaticRealVector<2> upper, StaticIntVector<2> n, tf::Executor& executor)
        : _lower(lower), _upper(upper), _n(n)
    {
        _compute(std::forward<FunctionT>(f), executor);
    }

    ///
    /// \brief  Evaluates the approximation at arguments _x_ and _y_.
    ///
    T operator()(Real x, Real y) const
    {
        const Real z1 = (2.0*x - _lower[0] - _upper[0])/(_upper[0] - _lower[0]);
        const Real z2 = (2.0*y - _lower[1] - _upper[1])/(_upper[1] - _lower[1]);

        std::vector<T> tmp(_evaluatedRows);
        for (int i = 0; i < _evaluatedRows; ++i)
        {
            // tmp[i] = chebyshevClenshawRecurrence(&_c[i*_n[1]], _evaluatedCols[i], z2);
            if (z2 < -0.95 || z2 > 0.95)
            {
                tmp[i] = Detail::modifiedChebyshevClenshawRecurrence(&_c[i*_n[1]], _evaluatedCols[i], z2);
            }
            else
            {
                tmp[i] = Detail::chebyshevClenshawRecurrence(&_c[i*_n[1]], _evaluatedCols[i], z2);
            }
        }

        // return chebyshevClenshawRecurrence(&tmp[0], _evaluatedRows, z1);
        if (z1 < -0.95 || z1 > 0.95)
        {
            return Detail::modifiedChebyshevClenshawRecurrence(&tmp[0], _evaluatedRows, z1);
        }
        else
        {
            return Detail::chebyshevClenshawRecurrence(&tmp[0], _evaluatedRows, z1);
        }
    }

    ///
    /// \brief      Removes all coefficients (starting from the highest index) with absolute value less then _norm_ until the first coefficient greater than _norm_ appears.
    /// 
    /// \copybrief setMaxNormCoefficients
    /// This is useful to speed up the evaluation of the Chebyshev polynomial later.
    ///
    void setMaxNormCoefficients(Real norm)
    {
        if (_n[0] < 2 || _n[1] < 2)
        {
            throw std::runtime_error("Cheb2::setMaxNormCoefficients(): Can't use on empty object.");
        }

        //
        // Compute _evaluatedRows
        //
        int i = _n[0] - 1;
        for (; i > 1; --i)
        {
            if (std::max(SciCore::maxNorm(_c[i*_n[1]]), SciCore::maxNorm(_c[i*_n[1]+1])) > norm)
            {
                break ;
            }
        }
        _evaluatedRows = i+1;

        //
        // For each row compute _evaluatedCols
        //
        for (i = 0; i < _evaluatedRows; ++i)
        {
            int j = _n[1] - 1;
            for (; j > 1; --j)
            {
                if (SciCore::maxNorm(_c[i*_n[1]+j]) > norm)
                {
                    break ;
                }
            }
            _evaluatedCols[i] = j+1;
        }
    }

    ///
    /// \brief  Computes the derivative of the represented function in _x_ direction and returns it as a new Cheb object.
    /// 
    /// \copybrief diffX
    /// 
    /// \warning Computing the derivative is numerically not very stable,
    ///          see also <a href="https://www.sciencedirect.com/science/article/pii/0021999192902743">On the errors incurred calculating derivatives using Chebyshev polynomials.</a>
    ///
    Cheb2 diffX() const
    {
        if (_n[0] < 2 || _n[1] < 2)
        {
            throw std::runtime_error("Cheb2::diffX(): Can't use on empty object.");
        }

        Cheb2 returnValue;
        returnValue._c.resize(_c.size());
        returnValue._lower = _lower;
        returnValue._upper = _upper;
        returnValue._n = _n;
        returnValue._evaluatedRows = _n[0];
        returnValue._evaluatedCols = std::vector<int>(_n[0], _n[1]);

        const Real con = 2.0 / (_upper[0] - _lower[0]);

        for (int j = 0; j < _n[1]; ++j)
        {
            if constexpr (IsMatrix<T>::value == true)
            {
                returnValue._c[_n[1]*(_n[0]-1) + j] = T::Zero(_c[0].rows(), _c[0].cols());
            }
            else // is Real or Complex
            {
                returnValue._c[_n[1]*(_n[0]-1) + j] = 0.0;
            }

            if(_n[0] > 1)
            {
                returnValue._c[_n[1]*(_n[0]-2) + j] = (2.0 *(_n[0]-1.0)) * _c[_n[1]*(_n[0]-1) + j];

                for(int i = _n[0]; i >= 3; --i) 
                {
                    returnValue._c[_n[1]*(i-3) + j] = returnValue._c[_n[1]*(i-1) + j] + (2.0 *(i-2.0)) * _c[_n[1]*(i-2) + j];
                }

                for(int i = 0  ; i < _n[0] ; ++i) 
                {
                    returnValue._c[_n[1]*i + j] *= con;
                }
            }
        }

        return returnValue;
    }

    ///
    /// \brief  Computes the derivative of the represented function in _y_ direction and returns it as a new Cheb object.
    /// 
    /// \copybrief diffY
    /// 
    /// \warning Computing the derivative is numerically not very stable,
    ///          see also <a href="https://www.sciencedirect.com/science/article/pii/0021999192902743">On the errors incurred calculating derivatives using Chebyshev polynomials.</a>
    ///
    Cheb2 diffY() const
    {
        if (_n[0] < 2 || _n[1] < 2)
        {
            throw std::runtime_error("Cheb2::diffY(): Can't use on empty object.");
        }

        Cheb2 returnValue;
        returnValue._c.resize(_c.size());
        returnValue._lower = _lower;
        returnValue._upper = _upper;
        returnValue._n = _n;
        returnValue._evaluatedRows = _n[0];
        returnValue._evaluatedCols = std::vector<int>(_n[0], _n[1]);

        const Real con = 2.0 / (_upper[1] - _lower[1]);

        for (int i = 0; i < _n[0]; ++i)
        {
            if constexpr (IsMatrix<T>::value == true)
            {
                returnValue._c[i*_n[1] + _n[1]-1] = T::Zero(_c[0].rows(), _c[0].cols());
            }
            else // is Real or Complex
            {
                returnValue._c[i*_n[1] + _n[1]-1] = 0.0;
            }

            if(_n[1] > 1)
            {
                returnValue._c[i*_n[1] + _n[1]-2] = (2.0 *(_n[1]-1.0)) * _c[i*_n[1] + _n[1]-1];

                for(int j = _n[1]; j >= 3; --j)
                {
                    returnValue._c[i*_n[1] + j-3] = returnValue._c[i*_n[1] + j-1] + (2.0 *(j-2.0)) * _c[i*_n[1] + j-2];
                }

                for(int j = 0  ; j < _n[1] ; ++j)
                {
                    returnValue._c[i*_n[1] + j] *= con;
                }
            }
        }

        return returnValue;
    }

    ///
    /// \brief  Returns the _x_ integral of the represented function.
    /// 
    /// \copybrief integrateX
    /// If the object represents the function \f$f(x,y)\f$, then the antiderivative defined by
    /// \f[
    ///     \int_{a}^x dx' f(x', y),
    /// \f]
    /// is returned where \a a corresponds to \a lowerLimit()[0].
    /// 
    Cheb2 integrateX() const
    {
        if (_n[0] < 2 || _n[1] < 2)
        {
            throw std::runtime_error("Cheb2::integrateX(): Can't use on empty object.");
        }

        Cheb2 returnValue;
        returnValue._c.resize(_c.size());
        returnValue._lower = _lower;
        returnValue._upper = _upper;
        returnValue._n = _n;
        returnValue._evaluatedRows = _n[0];
        returnValue._evaluatedCols = std::vector<int>(_n[0], _n[1]);

        const Real con = 0.25 * (_upper[0] - _lower[0]);

        if(_n[0] == 2)
        {
            for (int j = 0; j < _n[1]; ++j)
            {
                returnValue._c[_n[1]+j] = con * _c[j];
                returnValue._c[j] = 2.0 * returnValue._c[_n[1]+j];
            }
        }
        else
        {
            for (int j = 0; j < _n[1]; ++j)
            {
                T sum;
                if constexpr (IsMatrix<T>::value == true)
                {
                    sum = T::Zero(_c[0].rows(), _c[0].cols());
                }
                else // is Real or Complex
                {
                    sum = 0.0;
                }

                Real fac = 1.0;
                for(int i = 1; i <= _n[0]-2; ++i)
                {
                    returnValue._c[_n[1]*i + j] = (con/((Real)i)) * (_c[_n[1]*(i-1) + j] - _c[_n[1]*(i+1) + j]);
                    sum += fac * returnValue._c[_n[1]*i + j];
                    fac = -fac;
                }
                returnValue._c[_n[1]*(_n[0]-1) + j] = (con/(_n[0]-1.0)) * _c[_n[1]*(_n[0]-2) + j];
                sum += fac * returnValue._c[_n[1]*(_n[0]-1) + j];
                returnValue._c[j] = 2.0 * sum;
            }
        }

        return returnValue;
    }

    ///
    /// \brief  Returns the _y_ integral of the represented function.
    /// 
    /// \copybrief integrateY
    /// If the object represents the function \f$f(x,y)\f$, then the antiderivative defined by
    /// \f[
    ///     \int_{a}^y dy' f(x, y'),
    /// \f]
    /// is returned where \a a corresponds to \a lowerLimit()[1].
    /// 
    Cheb2 integrateY() const
    {
        if (_n[0] < 2 || _n[1] < 2)
        {
            throw std::runtime_error("Cheb2::integrateX(): Can't use on empty object.");
        }

        Cheb2 returnValue;
        returnValue._c.resize(_c.size());
        returnValue._lower = _lower;
        returnValue._upper = _upper;
        returnValue._n = _n;
        returnValue._evaluatedRows = _n[0];
        returnValue._evaluatedCols = std::vector<int>(_n[0], _n[1]);

        const Real con = 0.25 * (_upper[1] - _lower[1]);

        if(_n[1] == 2)
        {
            for (int i = 0; i < _n[0]; ++i)
            {
                returnValue._c[i*_n[1] + 1] = con * _c[i*_n[1]];
                returnValue._c[i*_n[1] + 0] = 2.0 * returnValue._c[i*_n[1] + 1];
            }
        }
        else
        {
            for (int i = 0; i < _n[0]; ++i)
            {
                T sum;
                if constexpr (IsMatrix<T>::value == true)
                {
                    sum = T::Zero(_c[i*_n[1] + 0].rows(), _c[i*_n[1] + 0].cols());
                }
                else // is Real or Complex
                {
                    sum = 0.0;
                }

                Real fac = 1.0;
                for(int j = 1; j <= _n[1]-2; ++j)
                {
                    returnValue._c[i*_n[1] + j] = (con/((Real)j)) * (_c[i*_n[1] + j-1] - _c[i*_n[1] + j+1]);
                    sum += fac * returnValue._c[i*_n[1] + j];
                    fac = -fac;
                }
                returnValue._c[i*_n[1] + _n[1]-1] = (con/(_n[1]-1.0)) * _c[i*_n[1] + _n[1]-2];
                sum += fac * returnValue._c[i*_n[1] + _n[1]-1];
                returnValue._c[i*_n[1] + 0] = 2.0 * sum;
            }
        }

        return returnValue;
    }

    ///
    /// \brief  Returns the lower limits of the approximation.
    ///
    StaticRealVector<2> lowerLimit() const noexcept
    {
        return _lower;
    }

    ///
    /// \brief  Returns the upper limits of the approximation.
    ///
    StaticRealVector<2> upperLimit() const noexcept
    {
        return _upper;
    }

    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(CEREAL_NVP(_c), CEREAL_NVP(_lower), CEREAL_NVP(_upper), CEREAL_NVP(_n), CEREAL_NVP(_evaluatedRows), CEREAL_NVP(_evaluatedCols));
    }

private:

    std::vector<T> _c; // Coefficients saved in row major storage as, c_{i j} = c[i*_n[1] + j]
    StaticRealVector<2> _lower;
    StaticRealVector<2> _upper;
    StaticIntVector<2> _n; // _n[0] is number of rows in _c, _n[1] is number of columns in _c

    // We allow to evalute less rows (columns per row) than we computed.
    // The idea is that we first compute _n[0] rows and _n[1] columns, check which
    // matrix elements are smaller than some error bound, and then ignore all these
    // when we evaluate the approximation.
    int _evaluatedRows;
    std::vector<int> _evaluatedCols; // in row i evaluate _evaluatedCols [i] elements

    template <typename FunctionT>
    void _compute(FunctionT&& f)
    {
        _c.resize(_n[0]*_n[1]);

        const Real pi = std::numbers::pi_v<Real>;
        const Real bma1 = (_upper[0] - _lower[0])/2;
        const Real bma2 = (_upper[1] - _lower[1])/2;
        const Real bpa1 = (_upper[0] + _lower[0])/2;
        const Real bpa2 = (_upper[1] + _lower[1])/2;

        // Sample the function at the Chebyshev points
        for (int i = 0; i < _n[0]; ++i)
        {
            for (int j = 0; j < _n[1]; ++j)
            {
                _c[i*_n[1] + j] = f(cos(pi*(i+0.5)/_n[0])*bma1 + bpa1, cos(pi*(j+0.5)/_n[1])*bma2 + bpa2);
            }
        }

        std::vector<T> buffer(std::max(_n[0], _n[1]));
        dct2<T, Eigen::StorageOptions::RowMajor>(_c.data(), _n[0], _n[1], buffer.data());
        for (auto& coeff : _c)
        {
            coeff *= 4.0/(_n[0]*_n[1]);
        }

        _evaluatedRows = _n[0];
        _evaluatedCols = std::vector<int>(_n[0], _n[1]);
    }

    template <typename FunctionT>
    void _compute(FunctionT&& f, tf::Executor& executor)
    {
        _c.resize(_n[0]*_n[1]);

        const Real pi = std::numbers::pi_v<Real>;
        const Real bma1 = (_upper[0] - _lower[0])/2;
        const Real bma2 = (_upper[1] - _lower[1])/2;
        const Real bpa1 = (_upper[0] + _lower[0])/2;
        const Real bpa2 = (_upper[1] + _lower[1])/2;

        // Sample the function at the Chebyshev points
        parallelFor(
            [this, &f, pi, bma1, bma2, bpa1, bpa2](int i)
            {
                for (int j = 0; j < _n[1]; ++j)
                {
                    _c[i*_n[1] + j] = f(cos(pi*(i+0.5)/_n[0])*bma1 + bpa1, cos(pi*(j+0.5)/_n[1])*bma2 + bpa2);
                }
            }, 0, _n[0], executor.num_workers(), executor);

        dct2Parallel<T, Eigen::StorageOptions::RowMajor>(_c.data(), _n[0], _n[1], executor);

        for (auto& coeff : _c)
        {
            coeff *= 4.0/(_n[0]*_n[1]);
        }

        _evaluatedRows = _n[0];
        _evaluatedCols = std::vector<int>(_n[0], _n[1]);
    }
};

template <typename FunctionT>
Cheb2(FunctionT&& f, StaticRealVector<2> lower, StaticRealVector<2> upper, StaticIntVector<2> n)
    -> Cheb2<std::invoke_result_t<FunctionT, Real, Real>>;

template <typename FunctionT>
Cheb2(FunctionT&& f, StaticRealVector<2> lower, StaticRealVector<2> upper, StaticIntVector<2> n, tf::Executor& executor)
    -> Cheb2<std::invoke_result_t<FunctionT, Real, Real>>;

} // namespace SciCore

#endif // SCICORE_CHEB2_H
