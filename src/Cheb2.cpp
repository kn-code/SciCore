#include "SciCore/Cheb2.h"

namespace SciCore {

#ifndef SCICORE_DONT_PRECOMPILE_TEMPLATES
template class Cheb2<Real>;
template class Cheb2<Complex>;
template class Cheb2<RealVector>;
template class Cheb2<Vector>;
template class Cheb2<RealMatrix>;
template class Cheb2<Matrix>;
#endif // SCICORE_DONT_PRECOMPILE_TEMPLATES

} // namespace SciCore
