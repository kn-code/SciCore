//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>

#include <SciCore/Cheb.h>
#include <SciCore/Random.h>

using namespace SciCore;

// Test cases from J. L. Aurentz and L. N. Trefethen, Chopping a Chebyshev series, http://arxiv.org/abs/1512.01803, p. 18
TEST(chopChebSeriesRelative, ExamplesFromPaper)
{
    Real epsilon = std::pow(2.0, -52);

    RealVector coeffs(50);
    RealVector coeffs2(50);

    for (int i = 0; i < coeffs.size(); ++i)
    {
        coeffs[i]  = std::pow(10.0, -(i + 1));
        coeffs2[i] = std::cos(std::pow(i + 1.0, 2));
    }

    EXPECT_EQ(Detail::chopChebSeriesRelative(coeffs.data(), coeffs.size(), epsilon), 18);

    RealVector test2 = coeffs + 1e-16 * coeffs2;
    RealVector test3 = coeffs + 1e-13 * coeffs2;
    RealVector test4 = coeffs + 1e-10 * coeffs2;

    EXPECT_EQ(Detail::chopChebSeriesRelative(test2.data(), test2.size(), epsilon), 15);
    EXPECT_EQ(Detail::chopChebSeriesRelative(test3.data(), test3.size(), epsilon), 13);
    EXPECT_EQ(Detail::chopChebSeriesRelative(test4.data(), test4.size(), epsilon), 50);
    EXPECT_EQ(Detail::chopChebSeriesRelative(test4.data(), test4.size(), 1e-10), 10);
}

// Check that 3 exp(-1/(x + 1)) - (x + 1) has same number of coefficients as in J. L. Aurentz and L. N. Trefethen, Chopping a Chebyshev series, http://arxiv.org/abs/1512.01803, p. 2, p. 7
TEST(chopChebSeriesRelative, ExamplesFromPaper2)
{
    Real epsilon = std::pow(2.0, -52);

    auto f = [](Real x) -> Real
    {
        return 3 * std::exp(-1 / (x + 1)) - (x + 1);
    };

    auto fLarge = [&](Real x) -> Real
    {
        return std::pow(2, 20) * f(x);
    };

    auto fSmall = [&](Real x) -> Real
    {
        return std::pow(2, -20) * f(x);
    };

    Cheb cheb(f, -1, 1, 512);
    Cheb chebLarge(fLarge, -1, 1, 512);
    Cheb chebSmall(fSmall, -1, 1, 512);

    EXPECT_EQ(cheb.chopCoefficients(0, epsilon), true);
    EXPECT_EQ((int)cheb.coefficients().size(), 166);

    EXPECT_EQ(chebLarge.chopCoefficients(0, epsilon), true);
    EXPECT_EQ((int)chebLarge.coefficients().size(), 166);

    EXPECT_EQ(chebSmall.chopCoefficients(0, epsilon), true);
    EXPECT_EQ((int)chebSmall.coefficients().size(), 166);

    EXPECT_EQ(cheb.chopCoefficients(0, 1e-6), true);
    EXPECT_EQ((int)cheb.coefficients().size(), 51);

    EXPECT_EQ(chebLarge.chopCoefficients(0, 1e-6), true);
    EXPECT_EQ((int)chebLarge.coefficients().size(), 51);

    EXPECT_EQ(chebSmall.chopCoefficients(0, 1e-6), true);
    EXPECT_EQ((int)chebSmall.coefficients().size(), 51);
}

TEST(chopChebSeriesRelative, RealMatrix)
{
    Real epsilon = std::pow(2.0, -52);

    auto f1 = [](Real x) -> RealMatrix
    {
        // 0, 1 + x + x * x, pow(x, 7) are captured exactly by low order Chebyshev approximation
        // Therefore number of coefficients should be same as for scalar function 3 * std::exp(-1 / (x + 1)) - (x + 1).
        return RealMatrix{
            {0, 3 * std::exp(-1 / (x + 1)) - (x + 1)},
            {1 + x + x * x, 1 + pow(x, 7)}
        };
    };

    Cheb cheb1(f1, -1, 1, 512);

    EXPECT_EQ(cheb1.chopCoefficients(0, epsilon), true);
    EXPECT_EQ((int)cheb1.coefficients().size(), 166);

    cheb1.chopCoefficients(0, 1e-13);
    RealVector xValues = RealVector::LinSpaced(500, -0.999, 1);
    for (Real x : xValues)
    {
        EXPECT_LT(relError(cheb1(x), f1(x)), 1e-12) << "x =" << x;
    }
}

TEST(Cheb, ConstructScalarSmoothFunction)
{
    Real a = 1;
    Real b = 10;
    int n  = 64;

    bool testMonotonicity = true;
    Real xOld             = a;
    auto f                = [&](Real x) -> Real
    {
        if (testMonotonicity == true)
        {
            EXPECT_GT(x, xOld); // Require that f is sampled at monotonically increasing arguments x
        }
        return 100 + sin(3.0 * x) * exp(-x);
    };

    Cheb interp(f, a, b, n);
    testMonotonicity = false;

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    RealVector fValues = interp(xValues);

    int i = 0;
    for (Real x : xValues)
    {
        EXPECT_LT(std::abs(f(x) - interp(x)), 1e-13);
        EXPECT_EQ(interp(x), fValues[i++]);
    }

    interp.chopCoefficients(0, 1e-11);
    for (Real x : xValues)
    {
        EXPECT_LT(relError(interp(x), f(x)), 1e-11);
    }
}

TEST(Cheb, ConstructScalarNoisyFunction)
{
    Real a = 1;
    Real b = 10;
    int n  = 128;

    Xoshiro256 rng(1234);
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    auto f = [](Real x) -> Real
    {
        return log(0.5 * x) + sin(3.0 * x) * exp(-x);
    };

    Real noise  = 1e-5;
    auto fNoisy = [&](Real x) -> Real
    {
        return log(0.5 * x) + sin(3.0 * x) * exp(-x) + noise * dis(rng);
    };

    // Constructing a "clean" function with epsAbs=10^-6 should require more coefficients
    // than constructing a "noisy" function with epsAbs=10^-6 but noise=1e-5, because
    // the Chebyshev interpolation should detect a plateau that sets in at the noise level
    {
        Real epsAbs = 1e-6;
        Cheb interp(f, a, b, n);
        interp.chopCoefficients(epsAbs, 0);

        noise = 1e-5;
        Cheb interpNoisy(fNoisy, a, b, n);
        interpNoisy.chopCoefficients(epsAbs, 0);

        EXPECT_GT(interp.coefficients().size(), interpNoisy.coefficients().size());
    }

    // Constructing a "clean" function with epsAbs=10^-6 should require same number of coefficients
    // than constructing a "noisy" function with epsAbs=10^-6 and noise=1e-6.
    {
        Real epsAbs = 1e-6;
        Cheb interp(f, a, b, n);
        interp.chopCoefficients(epsAbs, 0);

        noise = 1e-6;
        Cheb interpNoisy(fNoisy, a, b, n);
        interpNoisy.chopCoefficients(epsAbs, 0);

        EXPECT_EQ(interp.coefficients().size(), interpNoisy.coefficients().size());
    }
}

TEST(Cheb, ConstructLog)
{
    Real a = 0.01;
    Real b = 0.8;
    int n  = 256;

    auto f = [](Real x) -> Real
    {
        return std::log(x);
    };

    Cheb interp(f, a, b, n);

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    EXPECT_EQ(interp.chopCoefficients(0, 1e-11), true);
    for (Real x : xValues)
    {
        EXPECT_LT(relError(interp(x), f(x)), 1e-11) << "x =" << x << "f(x) =" << f(x) << "interp(x) = " << interp(x);
    }
}

// Scalar construction test for pow(I * t, 3.2)
TEST(Cheb, ConstructComplexPow)
{
    Real a = 0.1;
    Real b = 5;
    int n  = 128;

    auto f = [](Real x) -> Complex
    {
        Complex I(0, 1);
        return std::pow(I * x, 3.2);
    };

    Cheb interp(f, a, b, n);
    EXPECT_EQ(interp.chopCoefficients(0, 1e-9), true);

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    Vector fValues     = interp(xValues);

    int i = 0;
    for (Real x : xValues)
    {
        EXPECT_EQ(interp(x), fValues[i++]);
        EXPECT_LT(maxNorm(interp(x) - f(x)), 2e-8)
            << "x =" << x << "f(x) =" << f(x) << "interp(x) =" << interp(x) << " err =" << maxNorm(interp(x) - f(x));
    }
}

TEST(Cheb, ConstructFromArray)
{
    Real a = 1;
    Real b = 10;
    int n  = 64;

    auto f = [](Real x) -> Real
    {
        return sin(3.0 * x) * exp(-x);
    };

    RealVector cNodes = chebNodes(a, b, n);
    RealVector fVals  = cNodes.unaryExpr(f);
    Cheb interp(fVals.data(), a, b, n);

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(std::abs(f(x) - interp(x)), 1e-14);
    }

    interp.chopCoefficients(0, 1e-9);
    for (Real x : xValues)
    {
        // test that we are still accurate
        EXPECT_LT(maxNorm(interp(x) - f(x)), 1e-9) << "x =" << x << "f(x) =" << f(x);
    }
}

TEST(Cheb, ConstructMatrixFunction)
{
    using std::cos;
    using std::exp;
    using std::log;
    using std::pow;
    using std::sin;

    // clang-format off
//! [Basic cheb usage]
// Some complicated and hard to compute function
auto f = [](Real t) -> Matrix
{
    Complex I(0, 1);
    return Matrix{
        {cos(t) * I + sqrt(t + 1), exp(I * t + 1.0), 3 * t - I * t * t},
        {log(t), pow(I * t, 3.2), I * t * t * t - 1.0}
    };
};

Real a = 0.1;   // Lower limit of approdimation
Real b = 5;     // Upper limit of approdimation
int n  = 128;   // Number of Chebyshev coefficients to use

// Compute approximation
Cheb interp(f, a, b, n);
//! [Basic cheb usage]
    // clang-format on

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(maxNorm(f(x) - interp(x)), 2e-13) << "x =" << x;
    }

    EXPECT_EQ(interp.chopCoefficients(0, 1e-9), true);
    for (Real x : xValues)
    {
        EXPECT_LT(relError(interp(x), f(x)), 1e-8) << "x =" << x;
    }
}

TEST(Cheb, AdaptiveConstructionRealVectorFunction)
{
    auto f = [](Real x) -> RealVector
    {
        return RealVector{
            {1.0 / (x + 21.0), 1.0 / (x * x + 21.0)}
        };
    };

    Real a      = -20;
    Real b      = 20;
    Real epsAbs = 0;
    Real epsRel = 1e-8;
    bool ok     = false;
    Cheb c(f, a, b, epsAbs, epsRel, 16, 2, &ok);

    EXPECT_EQ(ok, true);

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(relError(c(x), f(x)), 1e-8) << "x =" << x;
    }

    // Construct again with parallel algorithm, check that they are exactly equal
    tf::Executor executor(4);
    Cheb c2(f, a, b, epsAbs, epsRel, executor, 16, 2, &ok);

    EXPECT_EQ(c, c2);
}

TEST(Cheb, AddRealMatrices)
{
    auto f = [](Real x) -> RealMatrix
    {
        return RealMatrix{
            {sin(3.0 * x), exp(-x)}
        };
    };

    auto g = [](Real x) -> RealMatrix
    {
        return RealMatrix{
            {cos(3.0 * x), x * x}
        };
    };

    Real a = 1;
    Real b = 10;
    int n  = 64;
    Cheb f_interpolation(f, a, b, n);
    Cheb g_interpolation(g, a, b, n);

    Cheb<RealMatrix> empty;
    empty           += f_interpolation;
    f_interpolation += g_interpolation;

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(maxNorm(empty(x) - f(x)), 1e-14);
        EXPECT_LT(maxNorm(f_interpolation(x) - f(x) - g(x)), 1e-13);
    }
}

TEST(Cheb, AddConstant)
{
    auto f = [](Real x) -> RealMatrix
    {
        return RealMatrix{
            {sin(3.0 * x), exp(-x)}
        };
    };

    RealMatrix A{
        {-5, 4}
    };

    auto fPlusA = [&](Real x) -> RealMatrix
    {
        return f(x) + A;
    };

    Real a = 1;
    Real b = 10;
    int n  = 64;
    Cheb f_interpolation(f, a, b, n);
    f_interpolation += A;

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(maxNorm(f_interpolation(x) - fPlusA(x)), 1e-3);
    }
}

TEST(Cheb, MultiplyScalarConstant)
{
    auto f = [](Real x) -> RealMatrix
    {
        return RealMatrix{
            {sin(3.0 * x), exp(-x)}
        };
    };

    Real alpha = 0.4;

    Real a = -1;
    Real b = 3;
    int n  = 64;
    Cheb f_interpolation(f, a, b, n);
    f_interpolation *= alpha;

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(maxNorm(f_interpolation(x) - alpha * f(x)), 1e-15);
    }
}

TEST(Cheb, IntegrateScalarSmoothFunction)
{
    auto f = [](Real x) -> Real
    {
        return sin(3.0 * x) * exp(x);
    };

    auto F = [](Real x)
    {
        return (3.0 * exp(1.0) * cos(3.0) - 3.0 * exp(x) * cos(3.0 * x) - exp(1.0) * sin(3.0) + exp(x) * sin(3.0 * x)) /
               10.0;
    };

    Real a = 1;
    Real b = 10;
    int n  = 64;
    Cheb f_interpolation(f, a, b, n);
    Cheb F_interpolation = f_interpolation.integrate();
    ;

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(std::abs(f(x) - f_interpolation(x)), 1e-10);
        EXPECT_LT(maxNorm(F(x) - F_interpolation(x)), 1e-10);
    }
}

TEST(Cheb, IntegrateMatrixFunction)
{
    auto f = [](Real x) -> RealMatrix
    {
        return RealMatrix{{sin(3.0 * x) * exp(x)}};
    };

    auto F = [](Real x) -> RealMatrix
    {
        return RealMatrix{
            {(3.0 * exp(1.0) * cos(3.0) - 3.0 * exp(x) * cos(3.0 * x) - exp(1.0) * sin(3.0) + exp(x) * sin(3.0 * x)) /
             10.0}};
    };

    Real a = 1;
    Real b = 10;
    int n  = 64;
    Cheb f_interpolation(f, a, b, n);
    Cheb F_interpolation = f_interpolation.integrate();

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(maxNorm(F(x) - F_interpolation(x)), 1e-10);
    }
}

TEST(Cheb, IntegrateSparseMatrixFunction)
{
    auto f = [](Real x) -> SparseRealMatrix
    {
        SparseRealMatrix A(5, 5);
        A.insert(3, 4) = sin(3.0 * x) * exp(x);
        return A;
    };

    auto F = [](Real x) -> SparseRealMatrix
    {
        SparseRealMatrix A(5, 5);
        A.insert(3, 4) =
            (3.0 * exp(1.0) * cos(3.0) - 3.0 * exp(x) * cos(3.0 * x) - exp(1.0) * sin(3.0) + exp(x) * sin(3.0 * x)) /
            10.0;
        return A;
    };

    Real a = 1;
    Real b = 10;
    int n  = 64;
    Cheb f_interpolation(f, a, b, n);
    Cheb F_interpolation = f_interpolation.integrate();

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(maxNorm(F(x) - F_interpolation(x)), 1e-10);
    }
}

TEST(Cheb, DifferentiateScalarFunction)
{
    auto f = [](Real x) -> Real
    {
        return sin(3.0 * x) * exp(x);
    };

    auto df = [](Real x) -> Real
    {
        return (3.0 * cos(3.0 * x) + sin(3.0 * x)) * exp(x);
    };

    Real a = 1;
    Real b = 10;
    int n  = 64;
    Cheb f_interpolation(f, a, b, n);
    Cheb f_diff = f_interpolation.diff();

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(maxNorm(df(x) - f_diff(x)), 1e-9);
    }
}

TEST(Cheb, DifferentiateMatrixFunction)
{
    auto f = [](Real x) -> RealMatrix
    {
        return RealMatrix{{sin(3.0 * x) * exp(x)}};
    };

    auto df = [](Real x) -> RealMatrix
    {
        return RealMatrix{{(3.0 * cos(3.0 * x) + sin(3.0 * x)) * exp(x)}};
    };

    Real a = 1;
    Real b = 10;
    int n  = 64;
    Cheb f_interpolation(f, a, b, n);
    Cheb f_diff = f_interpolation.diff();

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(maxNorm(df(x) - f_diff(x)), 2e-9);
    }
}

TEST(Cheb, DifferentiateSparseMatrixFunction)
{
    auto f = [](Real x) -> SparseRealMatrix
    {
        SparseRealMatrix A(5, 5);
        A.insert(3, 4) = sin(3.0 * x) * exp(x);
        return A;
    };

    auto df = [](Real x) -> SparseRealMatrix
    {
        SparseRealMatrix A(5, 5);
        A.insert(3, 4) = (3.0 * cos(3.0 * x) + sin(3.0 * x)) * exp(x);
        return A;
    };

    Real a = 1;
    Real b = 10;
    int n  = 64;
    Cheb f_interpolation(f, a, b, n);
    Cheb f_diff = f_interpolation.diff();

    RealVector xValues = RealVector::LinSpaced(100, a, b);
    for (Real x : xValues)
    {
        EXPECT_LT(maxNorm(df(x) - f_diff(x)), 2e-9);
    }
}

TEST(Cheb, Serialize)
{
    auto f = [](Real t) -> Matrix
    {
        Complex I(0, 1);
        return Matrix{
            {std::cos(t) * I + std::sqrt(t + 1), std::exp(I * t + Real(1)), 3 * t - I * t * t},
            {std::log(3 + t), std::pow(I * t, Real(3.2)), I * t * t * t - Real(1)}
        };
    };

    Real a = 0.1;
    Real b = 5;
    int n  = 128;

    Cheb interpolation(f, a, b, n);

    // Use binary file
    std::string binaryFilename = "cheb_test_serialization.binary";
    std::remove(binaryFilename.c_str());
    {
        std::ofstream os(binaryFilename, std::ios::binary);
        cereal::BinaryOutputArchive archive(os);
        archive(interpolation);
    }

    Cheb<Matrix> fromFile;
    {
        std::ifstream is(binaryFilename, std::ios::binary);
        cereal::BinaryInputArchive archive(is);
        archive(fromFile);
    }

    EXPECT_EQ(fromFile, interpolation);

    // Use JSON file
    std::string jsonFilename = "cheb_test_serialization.json";
    std::remove(jsonFilename.c_str());
    {
        std::ofstream os(jsonFilename);
        cereal::JSONOutputArchive archive(os);
        archive(cereal::make_nvp("cheb", interpolation));
        archive(interpolation);
    }

    {
        std::ifstream is(jsonFilename);
        cereal::JSONInputArchive archive(is);
        archive(cereal::make_nvp("cheb", fromFile));
    }

    EXPECT_EQ(fromFile, interpolation);

    // Serialize empty Cheb object
    jsonFilename   = "empty_cheb_test_serialization.json";
    binaryFilename = "empty_cheb_test_serialization.binary";

    Cheb<Matrix> empty;
    {
        std::ofstream os(jsonFilename);
        cereal::JSONOutputArchive archive(os);
        archive(cereal::make_nvp("cheb", empty));
    }

    Cheb<Matrix> fromJSONFile = interpolation;
    {
        std::ifstream is(jsonFilename);
        cereal::JSONInputArchive archive(is);
        archive(cereal::make_nvp("cheb", fromJSONFile));
    }

    EXPECT_EQ(fromJSONFile, empty);

    {
        std::ofstream os(binaryFilename, std::ios::binary);
        cereal::BinaryOutputArchive archive(os);
        archive(cereal::make_nvp("cheb", empty));
    }

    Cheb<Matrix> fromBinaryFile = interpolation;
    {
        std::ifstream is(binaryFilename, std::ios::binary);
        cereal::BinaryInputArchive archive(is);
        archive(fromBinaryFile);
    }

    EXPECT_EQ(fromBinaryFile, empty);
}