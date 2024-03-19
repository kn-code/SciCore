//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include "SciCore/DCT2.h"

namespace SciCore
{

// clang-format off
#ifndef SCICORE_DONT_PRECOMPILE_TEMPLATES

template void dct2<Real, Eigen::StorageOptions::ColMajor>(Real *, int, int);
template void dct2<Complex, Eigen::StorageOptions::ColMajor>(Complex *, int, int);
template void dct2<RealVector, Eigen::StorageOptions::ColMajor>(RealVector *, int, int);
template void dct2<Vector, Eigen::StorageOptions::ColMajor>(Vector *, int, int);
template void dct2<RealMatrix, Eigen::StorageOptions::ColMajor>(RealMatrix *, int, int);
template void dct2<Matrix, Eigen::StorageOptions::ColMajor>(Matrix *, int, int);

template void dct2<Real, Eigen::StorageOptions::RowMajor>(Real *, int, int);
template void dct2<Complex, Eigen::StorageOptions::RowMajor>(Complex *, int, int);
template void dct2<RealVector, Eigen::StorageOptions::RowMajor>(RealVector *, int, int);
template void dct2<Vector, Eigen::StorageOptions::RowMajor>(Vector *, int, int);
template void dct2<RealMatrix, Eigen::StorageOptions::RowMajor>(RealMatrix *, int, int);
template void dct2<Matrix, Eigen::StorageOptions::RowMajor>(Matrix *, int, int);

template void dct2Parallel<Real, Eigen::StorageOptions::ColMajor>(Real*, int, int, tf::Executor&);
template void dct2Parallel<Complex, Eigen::StorageOptions::ColMajor>(Complex*, int, int, tf::Executor&);
template void dct2Parallel<RealVector, Eigen::StorageOptions::ColMajor>(RealVector*, int, int, tf::Executor&);
template void dct2Parallel<Vector, Eigen::StorageOptions::ColMajor>(Vector*, int, int, tf::Executor&);
template void dct2Parallel<RealMatrix, Eigen::StorageOptions::ColMajor>(RealMatrix*, int, int, tf::Executor&);
template void dct2Parallel<Matrix, Eigen::StorageOptions::ColMajor>(Matrix*, int, int, tf::Executor&);

template void dct2Parallel<Real, Eigen::StorageOptions::RowMajor>(Real*, int, int, tf::Executor&);
template void dct2Parallel<Complex, Eigen::StorageOptions::RowMajor>(Complex*, int, int, tf::Executor&);
template void dct2Parallel<RealVector, Eigen::StorageOptions::RowMajor>(RealVector*, int, int, tf::Executor&);
template void dct2Parallel<Vector, Eigen::StorageOptions::RowMajor>(Vector*, int, int, tf::Executor&);
template void dct2Parallel<RealMatrix, Eigen::StorageOptions::RowMajor>(RealMatrix*, int, int, tf::Executor&);
template void dct2Parallel<Matrix, Eigen::StorageOptions::RowMajor>(Matrix*, int, int, tf::Executor&);

#endif // SCICORE_DONT_PRECOMPILE_TEMPLATES
// clang-format on

} // namespace SciCore
