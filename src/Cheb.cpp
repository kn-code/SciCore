//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include "SciCore/Cheb.h"

namespace SciCore
{

void chebNodes(Real a, Real b, int n, Real* nodes) noexcept
{
    assert(n > 0);
    for (int k = 0; k < n; ++k)
    {
        nodes[k] = (a + b) / 2 - (b - a) / 2 * std::cos(Real(2 * k + 1) / Real(2 * n) * std::numbers::pi_v<Real>);
    }
}

RealVector chebNodes(Real a, Real b, int n)
{
    assert(n > 0);
    RealVector returnValue(n);

    chebNodes(a, b, n, returnValue.data());

    return returnValue;
}

#ifndef SCICORE_DONT_PRECOMPILE_TEMPLATES
template class Cheb<Real>;
template class Cheb<Complex>;
template class Cheb<RealVector>;
template class Cheb<Vector>;
template class Cheb<RealMatrix>;
template class Cheb<Matrix>;
#endif // SCICORE_DONT_PRECOMPILE_TEMPLATES

} // namespace SciCore
