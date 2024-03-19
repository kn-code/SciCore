//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include <gtest/gtest.h>

#include <fstream>

#include <SciCore/ChebAdaptive.h>

using namespace SciCore;

using std::cos;
using std::cosh;
using std::exp;
using std::sin;

TEST(ChebAdaptive, ConstructScalarFunctionRelErr)
{
// clang-format off
//! [Basic adaptive cheb usage]
auto f = [](Real x) -> Real
{
    return cos(cos(x) + 3 * sin(x) + 2 * cos(2 * x) + 3 * cos(3 * x)) + 1 / (cosh(400 * (x - 0.4)));
};

Real a = 0;
Real b = 3.5;
Real epsAbs = 0;
Real epsRel = 1e-12;
ChebAdaptive chebf(f, a, b, epsAbs, epsRel, 0);
//! [Basic adaptive cheb usage]
// clang-format on

    std::cout << "Sections:" << chebf.sections().transpose() << "\n"
              << "Number coefficients: " << chebf.numCoefficients().transpose() << "\n";

    RealVector xValues = RealVector::LinSpaced(1000, a, b);
    for (Real x : xValues)
    {
        if (std::abs(f(x)) > 0.1)
        {
            EXPECT_LT(relError(chebf(x), f(x)), 3e-12)
                << "x =" << x << "f(x) =" << f(x) << "chebf.sections() =" << chebf.sections().transpose();
        }
        else
        {
            EXPECT_LT(relError(chebf(x), f(x)), 1e-11)
                << "x =" << x << "f(x) =" << f(x) << "chebf.sections() =" << chebf.sections().transpose();
        }
    }
}

TEST(ChebAdaptive, ConstructScalarFunctionAbsErr)
{
    auto f = [](Real x) -> Real
    {
        return cos(cos(x) + 3 * sin(x) + 2 * cos(2 * x) + 3 * cos(3 * x)) + 1 / (cosh(400 * (x - 0.4)));
    };

    Real a      = 0;
    Real b      = 3.5;
    Real epsAbs = 1e-12;
    Real epsRel = 0;
    ChebAdaptive chebf(f, a, b, epsAbs, epsRel, 0);

    RealVector xValues = RealVector::LinSpaced(1000, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(std::abs(chebf(x) - f(x)), 1e-12)
            << "x =" << x << "f(x) =" << f(x) << "chebf.sections() =" << chebf.sections().transpose();
    }
}

TEST(ChebAdaptive, ConstructWithSections)
{
    auto f = [](Real x) -> Real
    {
        return cos(cos(x) + 3 * sin(x) + 2 * cos(2 * x) + 3 * cos(3 * x)) + 1 / (cosh(400 * (x - 0.4)));
    };

    bool ok = false;
    RealVector sections{
        {0, 0.35, 0.45, 3.5}
    };
    Real epsAbs  = 0;
    Real epsRel  = 1e-12;
    Real hMin    = 0.01;
    ChebAdaptive chebf(f, sections, epsAbs, epsRel, hMin, &ok);

    EXPECT_EQ(ok, true);

    RealVector xValues = RealVector::LinSpaced(1000, sections[0], sections[sections.size() - 1]);
    for (Real x : xValues)
    {
        EXPECT_LT(relError(chebf(x), f(x)), 5e-12)
            << "x =" << x << "f(x) =" << f(x) << "chebf.sections() =" << chebf.sections().transpose()
            << "chebf.numCoefficients() =" << chebf.numCoefficients().transpose();
    }
}

TEST(ChebAdaptive, ConstructMatrixFunction)
{
    auto f = [](Real x) -> Matrix
    {
        return Matrix{
            {sin(x) / (x * x + 0.1), cos(x) / (x * x + 0.1)}
        };
    };

    Real a       = 0;
    Real b       = 1;
    Real epsAbs  = 0;
    Real epsRel  = 1e-12;
    Real hMin    = 0.01;
    ChebAdaptive chebf(f, a, b, epsAbs, epsRel, hMin);

    RealVector xValues = RealVector::LinSpaced(1000, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(maxNorm(f(x) - chebf(x)), 3e-12);
    }

    epsAbs = 1e-12;
    epsRel = 0;

    chebf = ChebAdaptive(f, a, b, epsAbs, epsRel, hMin);
    for (Real x : xValues)
    {
        EXPECT_LT(maxNorm(f(x) - chebf(x)), 1e-12);
    }
}

TEST(ChebAdaptive, ConstructFromPiecewiseChebs)
{
    auto f1 = [](Real x)
    {
        return sin(x);
    };

    auto f2 = [](Real x)
    {
        return x * x / (x + 1);
    };

    auto f3 = [](Real x)
    {
        return sin(x) * exp(-x);
    };

    RealVector sections{
        {0, 4, 5, 32}
    };
    std::vector<Cheb<Real>> singleChebs;
    singleChebs.push_back(Cheb(f1, sections[0], sections[1], 32));
    singleChebs.push_back(Cheb(f2, sections[1], sections[2], 32));
    singleChebs.push_back(Cheb(f3, sections[2], sections[3], 32));

    RealVector xValues1 = RealVector::LinSpaced(1000, sections[0], sections[1] - 1e-15);
    RealVector xValues2 = RealVector::LinSpaced(1000, sections[1], sections[2] - 1e-15);
    RealVector xValues3 = RealVector::LinSpaced(1000, sections[2], sections[3]);
    Real maxErr1        = 0;
    Real maxErr2        = 0;
    Real maxErr3        = 0;
    for (Real x : xValues1)
    {
        Real err = maxNorm(f1(x) - singleChebs[0](x));
        if (err > maxErr1)
        {
            maxErr1 = err;
        }
    }
    for (Real x : xValues2)
    {
        Real err = maxNorm(f2(x) - singleChebs[1](x));
        if (err > maxErr2)
        {
            maxErr2 = err;
        }
    }
    for (Real x : xValues3)
    {
        Real err = maxNorm(f3(x) - singleChebs[2](x));
        if (err > maxErr3)
        {
            maxErr3 = err;
        }
    }

    ChebAdaptive<Real> chebAll(singleChebs.data(), singleChebs.size());
    for (Real x : xValues1)
    {
        EXPECT_LE(maxNorm(f1(x) - chebAll(x)), maxErr1) << "x=" << x;
    }
    for (Real x : xValues2)
    {
        EXPECT_LE(maxNorm(f2(x) - chebAll(x)), maxErr2) << "x=" << x;
    }
    for (Real x : xValues3)
    {
        EXPECT_LE(maxNorm(f3(x) - chebAll(x)), maxErr3) << "x=" << x;
    }
}

TEST(ChebAdaptive, Serialization)
{
    auto f = [](Real x) -> Matrix
    {
        return Matrix{
            {sin(x) / (x * x + 0.1), cos(x) / (x * x + 0.1)}
        };
    };

    Real a       = 0;
    Real b       = 1;
    Real epsAbs  = 0;
    Real epsRel  = 1e-12;
    Real hMin    = 0.01;
    ChebAdaptive interpolation(f, a, b, epsAbs, epsRel, hMin);

    std::string archiveFilename = "chebAdaptive_test_out.cereal";
    std::remove(archiveFilename.c_str());

    // Serialize to file
    {
        std::ofstream os(archiveFilename, std::ios::binary);
        cereal::BinaryOutputArchive archive(os);
        archive(interpolation);
    };

    // Deserialize from file
    {
        ChebAdaptive<Matrix> fromFile;
        {
            std::ifstream is(archiveFilename, std::ios::binary);
            cereal::BinaryInputArchive archive(is);
            archive(fromFile);
        }

        EXPECT_EQ(fromFile, interpolation);
    };
}
