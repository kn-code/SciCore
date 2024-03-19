//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include <gtest/gtest.h>

#include <SciCore/BasicMath.h>

using namespace SciCore;

TEST(wrapAngle, Test)
{
    EXPECT_LE(abs(wrapAngle(-4 * std::numbers::pi_v<Real> + 0.5) - 0.5), 1e-17);
    EXPECT_LE(abs(wrapAngle(-2 * std::numbers::pi_v<Real> + 0.5) - 0.5), 1e-17);
    EXPECT_LE(abs(wrapAngle(0.5) - 0.5), 1e-17);
    EXPECT_LE(abs(wrapAngle(2 * std::numbers::pi_v<Real> + 0.5) - 0.5), 1e-17);
    EXPECT_LE(abs(wrapAngle(4 * std::numbers::pi_v<Real> + 0.5) - 0.5), 1e-17);

    EXPECT_LE(abs(wrapAngle(-4 * std::numbers::pi_v<Real>)), 1e-17);
    EXPECT_LE(abs(wrapAngle(-2 * std::numbers::pi_v<Real>)), 1e-17);
    EXPECT_LE(abs(wrapAngle(0.0)), 1e-17);
    EXPECT_LE(abs(wrapAngle(2 * std::numbers::pi_v<Real>)), 1e-17);
    EXPECT_LE(abs(wrapAngle(4 * std::numbers::pi_v<Real>)), 1e-17);
}

TEST(cosm1, Test)
{
    // Vector of Argument/Result
    std::vector<StaticRealVector<2>> testData{
        StaticRealVector<2>{1e-2, -0.49999583334722219742},
        StaticRealVector<2>{1e-4, -0.49999999958333333347},
        StaticRealVector<2>{1e-6, -0.49999999999995833333},
        StaticRealVector<2>{1e-8, -0.499999999999999995833333333333},
        StaticRealVector<2>{1e-10, -0.499999999999999999999583333333},
        StaticRealVector<2>{1e-12, -0.499999999999999999999999958333}};

    for (size_t i = 0; i < testData.size(); ++i)
    {
        Real x = testData[i][0];
        Real result = testData[i][1];
        EXPECT_LE(abs(cosm1(x) / (x * x) - result), 2e-16) << "x=" << x;
    }
}

TEST(expm1, Test)
{
    // Vector of Argument/Result
    std::vector<StaticVector<2>> testData{
        StaticVector<2>{Complex(1e-2, 1e-2), Complex(1.00499991633277777977, 0.00503341666610952183)},
        StaticVector<2>{Complex(1e-4, 1e-4), Complex(1.00004999999991666333, 0.00005000333341666667)},
        StaticVector<2>{Complex(1e-6, 1e-6), Complex(1.00000049999999999992, 5.0000033333342e-7)},
        StaticVector<2>{Complex(1e-8, 1e-8), Complex(1.0000000049999999999999999, 5.0000000333333334e-9)},
        StaticVector<2>{Complex(1e-10, 1e-10), Complex(1.00000000005, 5.0000000003333333333e-11)},
        StaticVector<2>{Complex(1e-12, 1e-12), Complex(1.0000000000005, 5.00000000000333333e-13)}};

    for (size_t i = 0; i < testData.size(); ++i)
    {
        Complex z = testData[i][0];
        Complex result = testData[i][1];
        EXPECT_LE(abs(expm1(z) / z - result), 4e-16) << "z=" << z;
    }

}
