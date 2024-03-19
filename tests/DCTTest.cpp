//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include <iostream>
#include <vector>

#include <gtest/gtest.h>

#include <SciCore/DCT.h>
#include <SciCore/DCT2.h>
#include <SciCore/Random.h>
#include <SciCore/Utility.h>

using namespace SciCore;

template <typename VectorT, typename ReturnT = VectorT>
ReturnT dctNaive(const VectorT& x)
{
    int N = std::ssize(x);
    ReturnT X(N);

    for (int k = 0; k < N; ++k)
    {
        X[k] = x[0] * std::cos(std::numbers::pi_v<Real> / Real(2 * N) * Real(k));
        for (int n = 1; n < N; ++n)
        {
            X[k] += x[n] * std::cos(std::numbers::pi_v<Real> / Real(N) * (Real(n + 0.5)) * Real(k));
        }
    }

    return X;
}

template <typename MatrixT>
MatrixT dct2Naive(const MatrixT& A)
{
    using std::cos;

    int N = A.rows();
    int M = A.cols();

    MatrixT returnValue(N, M);

    Real pi = std::numbers::pi_v<Real>;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            auto dct_ij  = A(0, 0);
            dct_ij      *= std::cos(i * Real(0.5) * pi / N) * std::cos(j * Real(0.5) * pi / M);
            for (int alpha = 0; alpha < N; ++alpha)
            {
                for (int beta = (alpha == 0) ? 1 : 0; beta < M; ++beta)
                {
                    dct_ij += (std::cos(pi * i / N * Real(alpha + 0.5)) * std::cos(pi * j / M * Real(beta + 0.5))) *
                              A(alpha, beta);
                }
            }
            returnValue(i, j) = dct_ij;
        }
    }

    return returnValue;
}

TEST(dct, CompareWithNaiveDCT)
{
    uint64_t seed = 1234;
    Xoshiro256 rng(seed);
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    std::vector<int> sizes{2, 3, 6, 8, 9, 16, 24, 128, 256};

    for (size_t i = 0; i < sizes.size(); ++i)
    {
        int size = sizes[i];

        RealVector data  = randomVector<RealVector>(size, dis, rng);
        RealVector data2 = data;

        dct(data.data(), data.size());
        data2 = dctNaive(data2);

        for (int j = 0; j < std::ssize(data); ++j)
        {
            EXPECT_LT(std::abs(data[j] - data2[j]), 9e-13) << "size = " << size;
        }
    }
}

TEST(dct, SparseMatrix)
{
    std::vector<Complex> scalarInput{Complex(2, 1), Complex(3, 2), Complex(5, -4),  Complex(7, -1),
                                      11.0,          -13.0,         Complex(17, -2), Complex(19, 1)};

    std::vector<Matrix> denseInput1{
        Matrix::Constant(2, 2, scalarInput[0]), Matrix::Constant(2, 2, scalarInput[1]),
        Matrix::Constant(2, 2, scalarInput[2]), Matrix::Constant(2, 2, scalarInput[3]),
        Matrix::Constant(2, 2, scalarInput[4]), Matrix::Constant(2, 2, scalarInput[5]),
        Matrix::Constant(2, 2, scalarInput[6]), Matrix::Constant(2, 2, scalarInput[7]),
    };

    std::vector<SparseMatrix> sparseInputCopy;
    for (int i = 0; i < std::ssize(denseInput1); ++i)
    {
        sparseInputCopy.push_back(denseInput1[i].sparseView());
    }

    std::vector<SparseMatrix> sparseInputCopy1;

    for (int i = 0; i < std::ssize(denseInput1); ++i)
    {
        sparseInputCopy1.push_back(denseInput1[i].sparseView());
    }

    std::vector<SparseMatrix> sparseInput = sparseInputCopy1;
    std::vector<SparseMatrix> buffer(sparseInputCopy1.size());

    dct(sparseInputCopy1.data(), sparseInputCopy1.size(), buffer.data());
    sparseInput = dctNaive(sparseInput);

    for (int i = 0; i < std::ssize(sparseInputCopy); ++i)
    {
        SparseMatrix diff = sparseInputCopy1[i] - sparseInput[i];
        EXPECT_LT(diff.coeffs().abs().maxCoeff(), 1e-13);
    }
}

TEST(dct, StridedInput)
{
    std::vector<Complex> scalarInput1{Complex(2, 1), Complex(3, 2), Complex(5, -4),  Complex(7, -1),
                                      11.0,          -13.0,         Complex(17, -2), Complex(19, 1)};
    Vector scalarResult = dctNaive<std::vector<Complex>, Vector>(scalarInput1);

    Vector scalarInputStride2{
        {Complex(2, 1), 1.0, Complex(3, 2), 2.0, Complex(5, -4), 3.0, Complex(7, -1), 4.0, 11.0, 5.0, -13.0, 6.0,
         Complex(17, -2), 7.0, Complex(19, 1), 8.0}
    };

    Vector scalarResultStride2{
        {scalarResult[0], 1.0, scalarResult[1], 2.0, scalarResult[2], 3.0, scalarResult[3], 4.0, scalarResult[4], 5.0,
         scalarResult[5], 6.0, scalarResult[6], 7.0, scalarResult[7], 8.0}
    };

    Vector scalarInputStride3{
        {Complex(2, 1), 1.0, 1.0, Complex(3, 2), 2.0, 2.0, Complex(5, -4), 3.0, 3.0, Complex(7, -1), 4.0, 4.0,
         11.0, 5.0, 5.0, -13.0, 6.0, 6.0, Complex(17, -2), 7.0, 7.0, Complex(19, 1), 8.0, 8.0}
    };

    Vector scalarResultStride3{
        {scalarResult[0], 1.0, 1.0, scalarResult[1], 2.0, 2.0, scalarResult[2], 3.0, 3.0, scalarResult[3], 4.0, 4.0,
         scalarResult[4], 5.0, 5.0, scalarResult[5], 6.0, 6.0, scalarResult[6], 7.0, 7.0, scalarResult[7], 8.0, 8.0}
    };

    Vector result2 = scalarInputStride2;
    Vector buffer2(result2.size() / 2);
    dct(result2.data(), result2.size() / 2, buffer2.data(), 2, 1);
    EXPECT_LT(maxNorm(result2 - scalarResultStride2), 1e-13);

    Vector result3 = scalarInputStride3;
    Vector buffer3(result3.size() / 3);
    dct(result3.data(), result3.size() / 3, buffer3.data(), 3, 1);
    EXPECT_LT(maxNorm(result3 - scalarResultStride3), 1e-13);
}

TEST(dct2, CompareColMajorWithNaiveDCT2)
{
    uint64_t seed = 1234;
    Xoshiro256 rng(seed);
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    using RealColMajor    = Eigen::Matrix<Real, -1, -1, Eigen::ColMajor>;
    using ComplexColMajor = Eigen::Matrix<Complex, -1, -1, Eigen::ColMajor>;

    std::vector sizes{
        StaticIntVector<2>{2, 2},   StaticIntVector<2>{2, 4},   StaticIntVector<2>{4, 4},
        StaticIntVector<2>{3, 3},   StaticIntVector<2>{3, 4},   StaticIntVector<2>{6, 4},
        StaticIntVector<2>{8, 8},   StaticIntVector<2>{9, 12},  StaticIntVector<2>{16, 16},
        StaticIntVector<2>{32, 16}, StaticIntVector<2>{16, 32}, StaticIntVector<2>{64, 64}};

    for (size_t i = 0; i < sizes.size(); ++i)
    {
        StaticIntVector<2> size = sizes[i];

        const auto input1 = randomMatrix<RealColMajor>(size[0], size[1], dis, rng);
        const auto input2 = randomMatrix<ComplexColMajor>(size[0], size[1], dis, rng);

        RealColMajor naiveResult1    = dct2Naive(input1);
        ComplexColMajor naiveResult2 = dct2Naive(input2);

        RealColMajor result1    = input1;
        ComplexColMajor result2 = input2;
        dct2<Real, Eigen::StorageOptions::ColMajor>(result1.data(), size[0], size[1]);
        dct2<Complex, Eigen::StorageOptions::ColMajor>(result2.data(), size[0], size[1]);

        EXPECT_LT(maxNorm(result1 - naiveResult1), 1e-13 * std::max(size[0], size[1])) << "size = " << size[0] << ", " << size[1];
        EXPECT_LT(maxNorm(result2 - naiveResult2), 1e-13 * std::max(size[0], size[1])) << "size = " << size[0] << ", " << size[1];
    }
}

TEST(dct2, CompareRowMajorWithNaiveDCT2)
{
    uint64_t seed = 1234;
    Xoshiro256 rng(seed);
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    using RealRowMajor    = Eigen::Matrix<Real, -1, -1, Eigen::RowMajor>;
    using ComplexRowMajor = Eigen::Matrix<Complex, -1, -1, Eigen::RowMajor>;

    std::vector sizes{
        StaticIntVector<2>{2, 2},   StaticIntVector<2>{2, 4},   StaticIntVector<2>{4, 4},
        StaticIntVector<2>{3, 3},   StaticIntVector<2>{3, 4},   StaticIntVector<2>{6, 4},
        StaticIntVector<2>{8, 8},   StaticIntVector<2>{9, 12},  StaticIntVector<2>{16, 16},
        StaticIntVector<2>{32, 16}, StaticIntVector<2>{16, 32}, StaticIntVector<2>{64, 64}};

    for (size_t i = 0; i < sizes.size(); ++i)
    {
        StaticIntVector<2> size = sizes[i];

        const auto input1 = randomMatrix<RealRowMajor>(size[0], size[1], dis, rng);
        const auto input2 = randomMatrix<ComplexRowMajor>(size[0], size[1], dis, rng);

        RealRowMajor naiveResult1    = dct2Naive(input1);
        ComplexRowMajor naiveResult2 = dct2Naive(input2);

        RealRowMajor result1    = input1;
        ComplexRowMajor result2 = input2;
        dct2<Real, Eigen::StorageOptions::RowMajor>(result1.data(), size[0], size[1]);
        dct2<Complex, Eigen::StorageOptions::RowMajor>(result2.data(), size[0], size[1]);

        EXPECT_LT(maxNorm(result1 - naiveResult1), 1e-13 * std::max(size[0], size[1])) << "size = " << size[0] << ", " << size[1];
        EXPECT_LT(maxNorm(result2 - naiveResult2), 1e-13 * std::max(size[0], size[1])) << "size = " << size[0] << ", " << size[1];
    }
}

TEST(dct2Parallel, CompareRowMajorWithNaiveDCT2)
{
    uint64_t seed = 4321;
    Xoshiro256 rng(seed);
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    using RealRowMajor    = Eigen::Matrix<Real, -1, -1, Eigen::RowMajor>;
    using ComplexRowMajor = Eigen::Matrix<Complex, -1, -1, Eigen::RowMajor>;

    tf::Executor executor(8);

    std::vector sizes{
        StaticIntVector<2>{2, 2},   StaticIntVector<2>{2, 4},   StaticIntVector<2>{4, 4},
        StaticIntVector<2>{3, 3},   StaticIntVector<2>{3, 4},   StaticIntVector<2>{6, 4},
        StaticIntVector<2>{8, 8},   StaticIntVector<2>{9, 12},  StaticIntVector<2>{16, 16},
        StaticIntVector<2>{32, 16}, StaticIntVector<2>{16, 32}, StaticIntVector<2>{64, 64}};

    for (size_t i = 0; i < sizes.size(); ++i)
    {
        StaticIntVector<2> size = sizes[i];

        const auto input1 = randomMatrix<RealRowMajor>(size[0], size[1], dis, rng);
        const auto input2 = randomMatrix<ComplexRowMajor>(size[0], size[1], dis, rng);

        RealRowMajor naiveResult1    = dct2Naive(input1);
        ComplexRowMajor naiveResult2 = dct2Naive(input2);

        RealRowMajor resultParallel1    = input1;
        ComplexRowMajor resultParallel2 = input2;
        RealRowMajor resultSerial1      = input1;
        ComplexRowMajor resultSerial2   = input2;
        dct2Parallel<Real, Eigen::StorageOptions::RowMajor>(resultParallel1.data(), size[0], size[1], executor);
        dct2Parallel<Complex, Eigen::StorageOptions::RowMajor>(resultParallel2.data(), size[0], size[1], executor);
        dct2<Real, Eigen::StorageOptions::RowMajor>(resultSerial1.data(), size[0], size[1]);
        dct2<Complex, Eigen::StorageOptions::RowMajor>(resultSerial2.data(), size[0], size[1]);

        EXPECT_LT(maxNorm(resultParallel1 - naiveResult1), 1e-13 * std::max(size[0], size[1]));
        EXPECT_LT(maxNorm(resultParallel2 - naiveResult2), 1e-13 * std::max(size[0], size[1]));

        bool sameAsSerial1 = (resultSerial1 == resultParallel1);
        bool sameAsSerial2 = (resultSerial2 == resultParallel2);
        EXPECT_EQ(sameAsSerial1, true);
        EXPECT_EQ(sameAsSerial2, true);
    }
}
