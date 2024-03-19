//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include <gtest/gtest.h>

#include <SciCore/Cheb2.h>

using namespace SciCore;

TEST(Cheb2, ConstructScalarRealFunction)
{
    using std::sin;
    using std::cos;

//! [Basic cheb2 usage]
// Some function of 2 arguments
auto K = [](Real x, Real y)
{
    return sin(y)*x-cos(x)*y;
};

// Parameters: -1 <= x <= 2
//             -2 <= y <= 1
const StaticRealVector<2> lower{-1, -2};
const StaticRealVector<2> upper{ 2,  1};
const StaticIntVector<2>  n{16, 16};

// Compute approximation
Cheb2 c(K, lower, upper, n);
//! [Basic cheb2 usage]

    RealVector gridX = RealVector::LinSpaced(30, lower[0], upper[0]);
    RealVector gridY = RealVector::LinSpaced(30, lower[1], upper[1]);

    for (Real x : gridX)
    {
        for (Real y : gridY)
        {
            EXPECT_LT(std::abs(K(x,y) - c(x,y)), 1e-14);
        }
    }
}

TEST(Cheb2, ConstructScalarComplexFunction)
{
    auto K = [](Real x, Real y)
    {
        return Complex(std::sin(y)*x, y-std::cos(x));
    };

    const StaticRealVector<2> lower{-0.4, -0.5};
    const StaticRealVector<2> upper{ 2,  2.1};
    const StaticIntVector<2>  n{16, 16};

    Cheb2 c(K, lower, upper, n);

    RealVector gridX = RealVector::LinSpaced(30, lower[0], upper[0]);
    RealVector gridY = RealVector::LinSpaced(30, lower[1], upper[1]);

    for (Real x : gridX)
    {
        for (Real y : gridY)
        {
            EXPECT_LT(maxNorm(K(x,y) - c(x,y)), 1e-14);
        }
    }
}

TEST(Cheb2, ConstructScalarMatrixFunction)
{
    auto K = [](Real x, Real y)
    {
        return Matrix
        {
            {Complex(std::sin(y)*x, y-std::cos(x)), Complex(x,y), 0},
            {x*x*x+y, Complex(x,y), Complex(0,1)*std::cos(Complex(0,1)*x)*std::exp(x+y)},
            {Complex(std::sin(y)*std::sin(x*y), y-std::erf(x)), Complex(x,y), Complex(3*x+y*x, 2)}
        };
    };

    const StaticRealVector<2> lower{-0.5, -0.5};
    const StaticRealVector<2> upper{   2,  2};
    const StaticIntVector<2>  n{32, 32};

    Cheb2 c(K, lower, upper, n);

    RealVector gridX = RealVector::LinSpaced(30, lower[0], upper[0]);
    RealVector gridY = RealVector::LinSpaced(30, lower[1], upper[1]);

    for (Real x : gridX)
    {
        for (Real y : gridY)
        {
            EXPECT_LT(maxNorm(K(x,y) - c(x,y)), 1e-13) << "x=" << x << ", y=" << y;
        }
    }
}

TEST(Cheb2, DifferentiateScalarFunction)
{
    auto f = [](Real x, Real y)
    {
        return Complex(std::sin(y)*x, y-std::cos(x));
    };

    auto dfdx = [](Real x, Real y)
    {
        return Complex(std::sin(y), std::sin(x));
    };

    auto dfdy = [](Real x, Real y)
    {
        return Complex(std::cos(y)*x, 1.0);
    };

    const StaticRealVector<2> lower{-0.4, -0.5};
    const StaticRealVector<2> upper{   2,  1.8};
    const StaticIntVector<2>  n{16, 16};

    Cheb2 c(f, lower, upper, n);
    Cheb2 dcdx = c.diffX();
    Cheb2 dcdy = c.diffY();

    RealVector gridX = RealVector::LinSpaced(30, lower[0], upper[0]);
    RealVector gridY = RealVector::LinSpaced(30, lower[1], upper[1]);

    for (Real x : gridX)
    {
        for (Real y : gridY)
        {
            EXPECT_LT(maxNorm(dfdx(x,y) - dcdx(x,y)), 1e-13) << "x=" << x << ", y=" << y;
            EXPECT_LT(maxNorm(dfdy(x,y) - dcdy(x,y)), 1e-13) << "x=" << x << ", y=" << y;
        }
    }
}

TEST(Cheb2, IntegrateScalarFunction)
{
    auto f = [](Real x, Real y)
    {
        return Complex(std::sin(y)*x, y-std::cos(x));
    };

    auto Fx = [](Real x, Real y)
    {
        return Complex((-0.5 + x*x/2)*sin(y), y + x*y - sin(1) - sin(x));
    };

    auto Fy = [](Real x, Real y)
    {
        return Complex(x*(cos(2) - cos(y)), -2 + y*y/2 - (2 + y)*cos(x));
    };

    const StaticRealVector<2> lower{ -1,-2};
    const StaticRealVector<2> upper{1.5, 1};
    const StaticIntVector<2>  n{16, 16};

    Cheb2 c(f, lower, upper, n);
    Cheb2 cX = c.integrateX();
    Cheb2 cY = c.integrateY();

    RealVector gridX = RealVector::LinSpaced(30, lower[0], upper[0]);
    RealVector gridY = RealVector::LinSpaced(30, lower[1], upper[1]);

    for (Real x : gridX)
    {
        for (Real y : gridY)
        {
            EXPECT_LT(maxNorm(Fx(x,y) - cX(x,y)), 1e-14) << "x=" << x << ", y=" << y;
            EXPECT_LT(maxNorm(Fy(x,y) - cY(x,y)), 1e-14) << "x=" << x << ", y=" << y;
        }
    }
}

TEST(Cheb2, SerializeMatrixFunction)
{
    auto K = [](Real x, Real y)
    {
        return Matrix
        {
            {Complex(std::sin(y)*x, y-std::cos(x)), Complex(x,y), 0},
            {x*x*x+y, Complex(x,y), Complex(0,1)*std::cos(Complex(0,1)*x)*std::exp(x+y)},
            {Complex(std::sin(y)*std::sin(x*y), y-std::erf(x)), Complex(x,y), Complex(3*x+y*x, 2)}
        };
    };

    const StaticRealVector<2> lower{-0.5, -0.5};
    const StaticRealVector<2> upper{   2,  2};
    const StaticIntVector<2>  n{128, 128};

    const Cheb2 c(K, lower, upper, n);

    const std::string archiveFilename = "cheb2_test_out.cereal";
    std::remove(archiveFilename.c_str());

    // Serialize to file
    {
        std::ofstream os(archiveFilename, std::ios::binary);
        cereal::BinaryOutputArchive archive(os);
        archive(c);
    }

    // Deserialize from file
    {
        Cheb2<Matrix> fromFile;
        {
            std::ifstream is(archiveFilename, std::ios::binary);
            cereal::BinaryInputArchive archive(is);
            archive(fromFile);
        }

        EXPECT_EQ(fromFile, c);
    };
}
