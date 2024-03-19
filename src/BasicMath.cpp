//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include "SciCore/BasicMath.h"

namespace SciCore
{

Real wrapAngle(Real angle)
{
    angle = std::fmod(angle, 2 * std::numbers::pi_v<Real>);
    if (angle < 0)
    {
        angle += 2 * std::numbers::pi_v<Real>;
    }
    return angle;
}

Real cosm1(Real x)
{
    x = wrapAngle(x);
    if (std::abs(x) < Real(0.05))
    {
        Real x2 = x * x;
        return -x2 / (2 * (1 + x2 * (Real(1.0 / 12.0) + x2 * (Real(1.0 / 240.0) + x2 / 6048))));
    }
    else
    {
        return std::cos(x) - 1;
    }
}

Complex expm1(Complex z)
{
    Real x = z.real();
    Real y = z.imag();

    if ((std::abs(x) >= 1.0) || (std::abs(y) >= 1.0))
    {
        return std::exp(z) - 1.0;
    }

    return Complex(std::expm1(x) * std::cos(y) + SciCore::cosm1(y), std::exp(x) * std::sin(y));
}

} // namespace SciCore