//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

///
/// \file   DCT2.h
///
/// \brief  Fast discrete two-dimensional cosine transform.
///

#ifndef SCICORE_DCT2_H
#define SCICORE_DCT2_H

#include "DCT.h"
#include "Parallel.h"

namespace SciCore
{

///
/// \ingroup DCT
///
/// \brief      Computes the two dimensional discrete cosine transform of \a x.
///
/// \copybrief dct2()
/// This is defined as
///
/// \f[
///     X_{ij} = \sum_{n=0}^{N-1} \sum_{m=0}^{M-1} x_{nm} \cos \left[ \frac{\pi}{N} \left(n + \frac{1}{2} \right) i
///     \right] \cos \left[ \frac{\pi}{M} \left(m + \frac{1}{2} \right) j \right].
/// \f]
///
/// Here \a x points to a matrix with \a N rows and \a M columns stored with a given \a StorageOrder.
///
/// \param      x               The data that is transformed.
/// \param      N               Number of rows.
/// \param      M               Number of columns.
/// \param      buffer          A buffer of size \a max(N, M).
///
/// \tparam     T               Data type.
/// \tparam     StorageOrder    Either \a Eigen::StorageOptions::ColMajor or \a Eigen::StorageOptions::RowMajor.
///
template <MatrixOrScalarType T, int StorageOrder = Eigen::StorageOptions::ColMajor>
void dct2(T* x, int N, int M, T* buffer)
{
    if constexpr (StorageOrder == Eigen::StorageOptions::RowMajor)
    {
        // Apply DCT to each row
        T* start = x;
        for (int i = 0; i < N; ++i)
        {
            dct(start, M, buffer, 1, 1);
            start += M;
        }

        // Apply DCT to each column
        start = x;
        for (int j = 0; j < M; ++j)
        {
            dct(start, N, buffer, M, 1);
            start += 1;
        }
    }
    else if constexpr (StorageOrder == Eigen::StorageOptions::ColMajor)
    {
        // Apply DCT to each row
        T* start = x;
        for (int i = 0; i < N; ++i)
        {
            dct(start, M, buffer, N, 1);
            start += 1;
        }

        // Apply DCT to each column
        start = x;
        for (int j = 0; j < M; ++j)
        {
            dct(start, N, buffer, 1, 1);
            start += N;
        }
    }
    else
    {
        static_assert(sizeof(T) + 1 == 0, "Invalid storage order");
    }
}

///
/// \brief      Convenience version that automatically allocates the required buffer.
///
template <MatrixOrScalarType T, int StorageOrder = Eigen::StorageOptions::ColMajor>
void dct2(T* x, int N, int M)
{
    std::vector<T> buffer(std::max(N, M));
    dct2<T, StorageOrder>(x, N, M, buffer.data());
}

// clang-format off
#ifndef SCICORE_DONT_PRECOMPILE_TEMPLATES
extern template void dct2<Real, Eigen::StorageOptions::ColMajor>(Real*, int, int);
extern template void dct2<Complex, Eigen::StorageOptions::ColMajor>(Complex*, int, int);
extern template void dct2<RealVector, Eigen::StorageOptions::ColMajor>(RealVector*, int, int);
extern template void dct2<Vector, Eigen::StorageOptions::ColMajor>(Vector*, int, int);
extern template void dct2<RealMatrix, Eigen::StorageOptions::ColMajor>(RealMatrix*, int, int);
extern template void dct2<Matrix, Eigen::StorageOptions::ColMajor>(Matrix*, int, int);

extern template void dct2<Real, Eigen::StorageOptions::RowMajor>(Real*, int, int);
extern template void dct2<Complex, Eigen::StorageOptions::RowMajor>(Complex*, int, int);
extern template void dct2<RealVector, Eigen::StorageOptions::RowMajor>(RealVector*, int, int);
extern template void dct2<Vector, Eigen::StorageOptions::RowMajor>(Vector*, int, int);
extern template void dct2<RealMatrix, Eigen::StorageOptions::RowMajor>(RealMatrix*, int, int);
extern template void dct2<Matrix, Eigen::StorageOptions::RowMajor>(Matrix*, int, int);
#endif // SCICORE_DONT_PRECOMPILE_TEMPLATES
// clang-format on

///
/// \ingroup DCT
///
/// \brief      Computes the two dimensional discrete cosine transform of \a x in parallel.
///
/// Here \a x points to a matrix with \a N rows and \a M columns stored with a given \a StorageOrder.
///
/// \param      x               The data that is transformed.
/// \param      N               Number of rows.
/// \param      M               Number of columns.
/// \param      executor        Taskflow executor.
///
/// \tparam     T               Data type.
/// \tparam     StorageOrder    Either \a Eigen::StorageOptions::ColMajor or \a Eigen::StorageOptions::RowMajor.
///
template <MatrixOrScalarType T, int StorageOrder = Eigen::StorageOptions::ColMajor>
void dct2Parallel(T* x, int N, int M, tf::Executor& executor)
{
    if constexpr (StorageOrder == Eigen::StorageOptions::RowMajor)
    {
        parallelFor(
            [x, M](int i)
            {
                T* start  = x;
                start    += i * M;
                std::vector<T> buffer(M);
                dct(start, M, buffer.data(), 1, 1);
            },
            0, N, executor.num_workers(), executor);

        parallelFor(
            [x, M, N](int i)
            {
                T* start  = x;
                start    += i;
                std::vector<T> buffer(N);
                dct(start, N, buffer.data(), M, 1);
            },
            0, M, executor.num_workers(), executor);
    }
    else if constexpr (StorageOrder == Eigen::StorageOptions::ColMajor)
    {
        parallelFor(
            [x, M, N](int i)
            {
                T* start  = x;
                start    += i;
                std::vector<T> buffer(M);
                dct(start, M, buffer.data(), N, 1);
            },
            0, N, executor.num_workers(), executor);

        parallelFor(
            [x, N](int i)
            {
                T* start  = x;
                start    += i * N;
                std::vector<T> buffer(N);
                dct(start, N, buffer.data(), 1, 1);
            },
            0, M, executor.num_workers(), executor);
    }
    else
    {
        static_assert(sizeof(T) + 1 == 0, "Invalid storage order");
    }
}

// clang-format off
#ifndef SCICORE_DONT_PRECOMPILE_TEMPLATES
extern template void dct2Parallel<Real, Eigen::StorageOptions::ColMajor>(Real*, int, int, tf::Executor&);
extern template void dct2Parallel<Complex, Eigen::StorageOptions::ColMajor>(Complex*, int, int, tf::Executor&);
extern template void dct2Parallel<RealVector, Eigen::StorageOptions::ColMajor>(RealVector*, int, int, tf::Executor&);
extern template void dct2Parallel<Vector, Eigen::StorageOptions::ColMajor>(Vector*, int, int, tf::Executor&);
extern template void dct2Parallel<RealMatrix, Eigen::StorageOptions::ColMajor>(RealMatrix*, int, int, tf::Executor&);
extern template void dct2Parallel<Matrix, Eigen::StorageOptions::ColMajor>(Matrix*, int, int, tf::Executor&);

extern template void dct2Parallel<Real, Eigen::StorageOptions::RowMajor>(Real*, int, int, tf::Executor&);
extern template void dct2Parallel<Complex, Eigen::StorageOptions::RowMajor>(Complex*, int, int, tf::Executor&);
extern template void dct2Parallel<RealVector, Eigen::StorageOptions::RowMajor>(RealVector*, int, int, tf::Executor&);
extern template void dct2Parallel<Vector, Eigen::StorageOptions::RowMajor>(Vector*, int, int, tf::Executor&);
extern template void dct2Parallel<RealMatrix, Eigen::StorageOptions::RowMajor>(RealMatrix*, int, int, tf::Executor&);
extern template void dct2Parallel<Matrix, Eigen::StorageOptions::RowMajor>(Matrix*, int, int, tf::Executor&);
#endif // SCICORE_DONT_PRECOMPILE_TEMPLATES
// clang-format on

/// \} // end of DCT

} // namespace SciCore

#endif // SCICORE_DCT2_H
