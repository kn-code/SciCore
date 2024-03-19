//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

///
/// \file   BasicMath.h
///
/// \brief  Basic math functions.
///

#ifndef SCICORE_BASIC_MATH_H
#define SCICORE_BASIC_MATH_H

#include "Definitions.h"

namespace SciCore
{

///
/// \brief      Wraps the angle _angle_ into the range \f$[0, 2\pi)\f$.
///
Real wrapAngle(Real angle);

///
/// \brief      Computes \f$\cos(x)-1\f$ in a numerically stable way for small \f$x\f$.
///
Real cosm1(Real x);

///
/// \brief      Computes \f$\exp(z)-1\f$ in a numerically stable way for small complex \f$z\f$.
///
Complex expm1(Complex z);

} // namespace SciCore

#endif // SCICORE_BASIC_MATH_H
