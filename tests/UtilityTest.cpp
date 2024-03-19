//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include <gtest/gtest.h>

#include <SciCore/Utility.h>

using namespace SciCore;

TEST(maxNorm, RealVector)
{
    RealVector x{
        {1, 4, -7, 6}
    };

    EXPECT_EQ(maxNorm(x), 7);
}

TEST(maxNorm, StdVector)
{
    std::vector<Real> x{1, 4, -7, 6};

    EXPECT_EQ(maxNorm(x), 7);
}

TEST(maxNorm, VectorOfMatrices)
{
    Matrix A{
        {-1.,  5.,  3.},
        { 3., 12., -9.}
    };
    Matrix B{
        {-1., 0.,  1.},
        { 1., 1., 99.}
    };
    Matrix C{
        {-10.,  5., 23.},
        {102., 12., -9.}
    };
    std::vector<Matrix> x{A, B, C};

    EXPECT_EQ(maxNorm(x), 102);
}
