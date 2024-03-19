//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include "SciCore/Utility.h"

namespace SciCore
{

// clang-format off
#ifndef SCICORE_DONT_PRECOMPILE_TEMPLATES
template bool isFinite<RealVector>(const RealVector& );
template bool isFinite<Vector>(const Vector& );
template bool isFinite<RealMatrix>(const RealMatrix& );
template bool isFinite<Matrix>(const Matrix& );
#endif // SCICORE_DONT_PRECOMPILE_TEMPLATES
// clang-format on

Real relError(Real testValue, Real trueValue) noexcept
{
    // includes the case trueValue==testValue==0
    if (testValue == trueValue)
    {
        return 0;
    }
    else if (trueValue != 0)
    {
        return std::abs(testValue - trueValue) / std::abs(trueValue);
    }
    else
    {
        return std::numeric_limits<Real>::max();
    }
}

Real relError(Complex testValue, Complex trueValue) noexcept
{
    // includes the case trueValue==testValue==0
    if (testValue == trueValue)
    {
        return 0;
    }
    else if (trueValue != Real(0))
    {
        return std::abs(testValue - trueValue) / std::abs(trueValue);
    }
    else
    {
        return std::numeric_limits<Real>::max();
    }
}

size_t nextPowerOf2(size_t n)
{
    size_t maxAllowedN = std::numeric_limits<size_t>::max() / 2;

    if (n >= maxAllowedN)
    {
        throw std::domain_error("nextPowerOf2: argument loo large");
    }

    return (n == 0) ? 1 : std::pow(2, std::ceil(log(n) / log(2)));
}

} // namespace SciCore
