//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include "SciCore/DCT.h"

namespace SciCore
{

// clang-format off
#ifndef SCICORE_DONT_PRECOMPILE_TEMPLATES

template void dct<Real>(Real *, int, Real *);
template void dct<Complex>(Complex *, int, Complex *);
template void dct<RealVector>(RealVector *, int, RealVector *);
template void dct<Vector>(Vector *, int, Vector *);
template void dct<RealMatrix>(RealMatrix *, int, RealMatrix *);
template void dct<Matrix>(Matrix *, int, Matrix *);

template void dct<Real>(Real *, int);
template void dct<Complex>(Complex *, int);
template void dct<RealVector>(RealVector *, int);
template void dct<Vector>(Vector *, int);
template void dct<RealMatrix>(RealMatrix *, int);
template void dct<Matrix>(Matrix *, int);

#endif // SCICORE_DONT_PRECOMPILE_TEMPLATES
// clang-format on

} // namespace SciCore
