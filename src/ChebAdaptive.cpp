//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include "SciCore/ChebAdaptive.h"

namespace SciCore
{

#ifndef SCICORE_DONT_PRECOMPILE_TEMPLATES
template class ChebAdaptive<Real>;
template class ChebAdaptive<Complex>;
template class ChebAdaptive<RealVector>;
template class ChebAdaptive<Vector>;
template class ChebAdaptive<RealMatrix>;
template class ChebAdaptive<Matrix>;
#endif // SCICORE_DONT_PRECOMPILE_TEMPLATES

} // namespace SciCore

