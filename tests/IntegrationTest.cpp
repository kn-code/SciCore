//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include <gtest/gtest.h>

#include <iomanip>
#include <iostream>
#include <vector>

#include <SciCore/Integration.h>

using std::abs;
using std::exp;
using std::sin;
using std::cos;
using std::sqrt;

using namespace SciCore;

// p. 21 of https://arxiv.org/pdf/1006.3962.pdf
TEST(integrateAdaptive, BatteryTest)
{
    int counter = 0;

    Real a1     = 0;
    Real b1     = 1;
    Real exact1 = 1.7182818284590452354;
    auto f1     = [&](Real x) -> Real
    {
        ++counter;
        return exp(x);
    };

    Real a2     = 0;
    Real b2     = 1;
    Real exact2 = 0.7;
    auto f2     = [&](Real x) -> Real
    {
        ++counter;
        return (x > 0.3) ? 1.0 : 0.0;
    };

    Real a3     = 0;
    Real b3     = 1;
    Real exact3 = 2. / 3.;
    auto f3     = [&](Real x) -> Real
    {
        ++counter;
        return sqrt(x);
    };

    Real a4     = -1;
    Real b4     = 1;
    Real exact4 = 0.47942822668880166736;
    auto f4     = [&](Real x) -> Real
    {
        ++counter;
        return 23. / 25. * cosh(x) - cos(x);
    };

    Real a5     = -1;
    Real b5     = 1;
    Real exact5 = 1.5822329637296729331;
    auto f5     = [&](Real x) -> Real
    {
        ++counter;
        return 1.0 / (x * x * x * x + x * x + 0.9);
    };

    Real a6     = 0;
    Real b6     = 1;
    Real exact6 = 0.4;
    auto f6     = [&](Real x) -> Real
    {
        ++counter;
        return pow(x, 1.5);
    };

    Real a7     = 0;
    Real b7     = 1;
    Real exact7 = 2;
    auto f7     = [&](Real x) -> Real
    {
        ++counter;
        return 1.0 / sqrt(x);
    };

    Real a8     = 0;
    Real b8     = 1;
    Real exact8 = 0.86697298733991103757;
    auto f8     = [&](Real x) -> Real
    {
        ++counter;
        return 1.0 / (1 + x * x * x * x);
    };

    Real a9     = 0;
    Real b9     = 1;
    Real exact9 = 1.1547005383792515290;
    auto f9     = [&](Real x) -> Real
    {
        ++counter;
        return 2.0 / (2 + sin(10 * std::numbers::pi_v<Real> * x));
    };

    Real a10     = 0;
    Real b10     = 1;
    Real exact10 = 0.69314718055994530942;
    auto f10     = [&](Real x) -> Real
    {
        ++counter;
        return 1.0 / (1 + x);
    };

    Real a11     = 0;
    Real b11     = 1;
    Real exact11 = 0.37988549304172247537;
    auto f11     = [&](Real x) -> Real
    {
        ++counter;
        return 1.0 / (1 + exp(x));
    };

    Real a12     = 0;
    Real b12     = 1;
    Real exact12 = 0.77750463411224827642;
    auto f12     = [&](Real x) -> Real
    {
        ++counter;
        return x / (exp(x) - 1);
    };

    Real a13     = 0;
    Real b13     = 1;
    Real exact13 = 0.49898680869304550250;
    auto f13     = [&](Real x) -> Real
    {
        ++counter;
        return sin(100 * std::numbers::pi_v<Real> * x) / (std::numbers::pi_v<Real> * x);
    };

    Real a14     = 0;
    Real b14     = 10;
    Real exact14 = 0.5;
    auto f14     = [&](Real x) -> Real
    {
        ++counter;
        return sqrt(50) * exp(-50 * std::numbers::pi_v<Real> * x * x);
    };

    Real a15     = 0;
    Real b15     = 10;
    Real exact15 = 1.0;
    auto f15     = [&](Real x) -> Real
    {
        ++counter;
        return 25 * exp(-25 * x);
    };

    Real a16     = 0;
    Real b16     = 10;
    Real exact16 = 0.49936338107645674464;
    auto f16     = [&](Real x) -> Real
    {
        ++counter;
        return 50 / (std::numbers::pi_v<Real> * (2500 * x * x + 1));
    };

    Real a17     = 0;
    Real b17     = 1;
    Real exact17 = 0.49898680869304550250;
    auto f17     = [&](Real x) -> Real
    {
        ++counter;
        return 50 * sin(50 * std::numbers::pi_v<Real> * x) * sin(50 * std::numbers::pi_v<Real> * x) /
               (50 * 50 * std::numbers::pi_v<Real> * std::numbers::pi_v<Real> * x * x);
    };

    Real a18     = 0;
    Real b18     = std::numbers::pi_v<Real>;
    Real exact18 = 0.29101878286005269852;
    auto f18     = [&](Real x) -> Real
    {
        ++counter;
        return cos(cos(x) + 3 * sin(x) + 2 * cos(2 * x) + 3 * cos(3 * x));
    };

    Real a19     = 0;
    Real b19     = 1;
    Real exact19 = -1;
    auto f19     = [&](Real x) -> Real
    {
        ++counter;
        return std::log(x);
    };

    Real a20     = -1;
    Real b20     = 1;
    Real exact20 = 1.5643964440690497731;
    auto f20     = [&](Real x) -> Real
    {
        ++counter;
        return 1.0 / (1.005 + x * x);
    };

    Real a21     = 0;
    Real b21     = 1;
    Real exact21 = 0.16349494301863722618;
    auto f21     = [&](Real x) -> Real
    {
        ++counter;
        return 1.0 / (cosh(20 * (x - 2.0 / 10.0))) + 1.0 / (cosh(400 * (x - 4.0 / 10.0))) +
               1.0 / (cosh(8000 * (x - 6.0 / 10.0)));
    };

    Real a22     = 0;
    Real b22     = 1;
    Real exact22 = -0.63466518254339257343;
    auto f22     = [&](Real x) -> Real
    {
        ++counter;
        return 4 * std::numbers::pi_v<Real> * std::numbers::pi_v<Real> * x *
               sin(20 * std::numbers::pi_v<Real> * x) * cos(2 * std::numbers::pi_v<Real> * x);
    };

    Real a23     = 0;
    Real b23     = 1;
    Real exact23 = 0.013492485649467772692;
    auto f23     = [&](Real x) -> Real
    {
        ++counter;
        return 1.0 / (1 + (230 * x - 30) * (230 * x - 30));
    };

    int counterFun       = 0;
    Real epsAbs          = 0;
    Real epsRel          = 1e-9;
    auto testIntegration = [&](auto f, Real a, Real b, Real exact)
    {
        counter = 0;
        ++counterFun;
        auto [result, absErr] = SciCore::integrateAdaptive(f, a, b, epsAbs, epsRel, 0);
        // std::cout << counterFun << ": est rel error = " << absErr / result
        //          << ", real rel error = " << abs(exact - result) / abs(exact)
        //          << ", evaluations = " << counter << "\n";

        EXPECT_LE(std::abs((result - exact) / exact), 3 * epsRel)
            << " epsRel=" << epsRel << " fun=" << counterFun << " result=" << result << " exact=" << exact;
    };

    // Test with epsRel=1e-9
    counterFun = 0;
    epsRel     = 1e-9;
    testIntegration(f1, a1, b1, exact1);
    testIntegration(f2, a2, b2, exact2);
    testIntegration(f3, a3, b3, exact3);
    testIntegration(f4, a4, b4, exact4);
    testIntegration(f5, a5, b5, exact5);
    testIntegration(f6, a6, b6, exact6);
    testIntegration(f7, a7, b7, exact7);
    testIntegration(f8, a8, b8, exact8);
    testIntegration(f9, a9, b9, exact9);
    testIntegration(f10, a10, b10, exact10);
    testIntegration(f11, a11, b11, exact11);
    testIntegration(f12, a12, b12, exact12);
    testIntegration(f13, a13, b13, exact13);
    testIntegration(f14, a14, b14, exact14);
    testIntegration(f15, a15, b15, exact15);
    testIntegration(f16, a16, b16, exact16);
    testIntegration(f17, a17, b17, exact17);
    testIntegration(f18, a18, b18, exact18);
    testIntegration(f19, a19, b19, exact19);
    testIntegration(f20, a20, b20, exact20);
    testIntegration(f21, a21, b21, exact21);
    testIntegration(f22, a22, b22, exact22);
    testIntegration(f23, a23, b23, exact23);

    // Test with epsRel=1e-12
    std::cout << "\n";
    counterFun = 0;
    epsRel     = 1e-12;
    testIntegration(f1, a1, b1, exact1);
    testIntegration(f2, a2, b2, exact2);
    testIntegration(f3, a3, b3, exact3);
    testIntegration(f4, a4, b4, exact4);
    testIntegration(f5, a5, b5, exact5);
    testIntegration(f6, a6, b6, exact6);
    testIntegration(f7, a7, b7, exact7);
    testIntegration(f8, a8, b8, exact8);
    testIntegration(f9, a9, b9, exact9);
    testIntegration(f10, a10, b10, exact10);
    testIntegration(f11, a11, b11, exact11);
    testIntegration(f12, a12, b12, exact12);
    testIntegration(f13, a13, b13, exact13);
    testIntegration(f14, a14, b14, exact14);
    testIntegration(f15, a15, b15, exact15);
    testIntegration(f16, a16, b16, exact16);
    testIntegration(f17, a17, b17, exact17);
    testIntegration(f18, a18, b18, exact18);
    testIntegration(f19, a19, b19, exact19);
    testIntegration(f20, a20, b20, exact20);
    testIntegration(f21, a21, b21, exact21);
    testIntegration(f22, a22, b22, exact22);
    testIntegration(f23, a23, b23, exact23);
}

TEST(integrateAdaptive, OscillatingSlowlyDecaying)
{
    Real exact = 0.031395910020901238423;

// clang-format off
//! [integrate_adaptive_example]
int counter = 0;
auto f = [&](Real x) -> Real
{
    ++counter;
    return pow(sin(50 * x), 2) / pow(50 * x, 2);
};

Real a      = 0;
Real b      = 10;
Real epsAbs = 1e-10;
Real epsRel = 0;

auto [result, absErr] = integrateAdaptive(f, a, b, epsAbs, epsRel);
std::cout << std::setprecision(15) << "Integral is " << result << " +- " << absErr
          << " (computed with " << counter << " function evalations).\n";
//! [integrate_adaptive_example]
// clang-format on

    EXPECT_LE(abs(result - exact), absErr);
    EXPECT_LE(absErr, epsAbs);
}

TEST(integrateAdaptive, ScalarComplex)
{
    Complex exact(-0.38141464882934316933, 0.36783404124167983674);

    auto f = [](Real x) -> Complex
    {
        return Complex(sin(3.0 * x), exp(-x));
    };

    Real a                = 1;
    Real b                = 10;
    Real epsAbs           = 0;
    Real epsRel           = 1e-10;
    auto [result, absErr] = integrateAdaptive(f, a, b, epsAbs, epsRel);

    EXPECT_LE(abs((result - exact) / result), epsRel);
    EXPECT_LE(abs(result - exact), absErr);
}

TEST(integrateAdaptive, InitialSubsections)
{
    Real a      = 0;
    Real b      = 10;
    Real epsAbs = 0;
    Real epsRel = 1e-12;
    Real exact  = 0.49936338107645674464;

    int counter = 0;
    auto f      = [&](Real x) -> Real
    {
        ++counter;
        return 50 / (std::numbers::pi_v<Real> * (2500 * x * x + 1));
    };

    RealVector subsections{
        {a, b}
    };
    auto [result, absErr] = integrateAdaptive(f, subsections, epsAbs, epsRel);
    // std::cout << "Test subsection 1: est rel error = " << absErr/result << ", real rel error = " <<
    // abs(exact-result)/abs(exact) << ", evaluations = " << counter << "\n"; std::cout << "Subsections = " <<
    // subsections.transpose() << "\n\n";

    EXPECT_LE(abs((result - exact) / result), 3 * epsRel)
        << " epsRel=" << epsRel << " result=" << result << " exact=" << exact;
    EXPECT_LE(absErr / result, epsRel) << " epsRel=" << epsRel << " result=" << result << " exact=" << exact;

    counter             = 0;
    tie(result, absErr) = integrateAdaptive(f, subsections, epsAbs, epsRel);
    // std::cout << "Test subsection 2: est rel error = " << absErr/result << ", real rel error = " <<
    // abs(exact-result)/abs(exact) << ", evaluations = " << counter << "\n"; std::cout << "Subsections = " <<
    // subsections.transpose() << "\n\n";

    EXPECT_LE(abs((result - exact) / result), 3 * epsRel)
        << " epsRel=" << epsRel << " result=" << result << " exact=" << exact;
    EXPECT_LE(absErr / result, epsRel) << " epsRel=" << epsRel << " result=" << result << " exact=" << exact;

    counter             = 0;
    tie(result, absErr) = integrateAdaptive(f, subsections, epsAbs, epsRel);
    // std::cout << "Test subsection 2: est rel error = " << absErr/result << ", real rel error = " <<
    // abs(exact-result)/abs(exact) << ", evaluations = " << counter << "\n"; std::cout << "Subsections = " <<
    // subsections.transpose() << "\n\n";

    EXPECT_LE(abs((result - exact) / result), 3 * epsRel)
        << " epsRel=" << epsRel << " result=" << result << " exact=" << exact;
    EXPECT_LE(absErr / result, epsRel) << " epsRel=" << epsRel << " result=" << result << " exact=" << exact;
}

TEST(integrateTrapezoid, ScalarRealFunction)
{
    auto [result, errEstimate] =
        integrateTrapezoid([](Real x) -> Real { return exp(-x * x); }, 0, 3, 1e-13, 1e-13, 20, 1);
    Real exact = 0.88620734825952123389L;

    EXPECT_LE(abs(result - exact), 1e-13);
    EXPECT_LE(abs(result - exact), errEstimate);
}

TEST(integrateTrapezoid, StdVectorRealFunction)
{
    auto [result, errEstimate] = integrateTrapezoid(
        [](Real x) -> std::vector<Real> {
            return std::vector<Real>{exp(-x * x), exp(-2 * x * x)};
        },
        0, 3, 1e-13, 1e-13, 20, 1);
    std::vector<Real> exact{0.88620734825952123389L, 0.62665706742124588238L};

    EXPECT_LE(abs(result[0] - exact[0]), 1e-13);
    EXPECT_LE(abs(result[1] - exact[1]), errEstimate);
}

TEST(integrateTrapezoid, ScalarComplexFunction)
{
    auto [result, errEstimate] =
        integrateTrapezoid([](Real x) -> Complex { return exp(Complex(-x, x)); }, 0, 3, 1e-13, 1e-13, 25, 20);
    Complex exact(0.52815738780063440509L, 0.52113143631128428504L);

    EXPECT_LE(maxNorm(result - exact), 1e-13);
    EXPECT_LE(maxNorm(result - exact), errEstimate);
}

TEST(integrateTrapezoid, MatrixFunction)
{
    auto [result, errEstimate] = integrateTrapezoid(
        [](Real x) -> Matrix
        {
            return Matrix{
                {exp(Complex(-x, x)), exp(-x * x)}
            };
        },
        0, 3, 1e-13, 1e-13, 25, 20);
    Matrix exact{
        {Complex(0.52815738780063440509L, 0.52113143631128428504L), 0.88620734825952123389L}
    };

    EXPECT_LE(maxNorm(result - exact), 1e-13);
    EXPECT_LE(maxNorm(result - exact), 2 * errEstimate);
}

TEST(integrateTrapezoid, StdVectorRealMatrixFunction)
{
    auto [result, errEstimate] = integrateTrapezoid(
        [](Real x) -> std::vector<RealMatrix>
        {
            return std::vector<RealMatrix>{
                RealMatrix{{exp(-x * x), exp(-2 * x * x)}},
                RealMatrix{{exp(-3 * x * x)}, {exp(-4 * x * x)}}
            };
        },
        0, 3, 1e-13, 1e-13, 20, 1);
    std::vector<RealMatrix> exact{
        RealMatrix{{0.88620734825952123389L, 0.62665706742124588238L}},
        RealMatrix{{0.51166335397314166105L}, {0.44311346272637899729L}}
    };

    EXPECT_LE(maxNorm(result[0] - exact[0]), 1e-13);
    EXPECT_LE(maxNorm(result[1] - exact[1]), errEstimate);
}

TEST(integrateDE, ScalarRealFunction1)
{
    using std::exp;
    using std::sin;
    using std::cos;

    // We code the function twice, with one or two parameters
    // to check that both versions work.
    auto f1 = [](Real x) -> Real
    {
        return sin(3.0 * x) * exp(x);
    };

    auto f2 = [](Real x, Real) -> Real
    {
        return sin(3.0 * x) * exp(x);
    };

    auto F = [](Real x) -> Real
    {
        return (3.0 * exp(1.0) * cos(3.0) - 3.0 * exp(x) * cos(3.0 * x) - exp(1.0) * sin(3.0) +
                exp(x) * sin(3.0 * x)) /
               10.0;
    };

    Real eps = 1e-10;

    auto [result1, error1] = integrateDE(f1, 1.0, 10.0, eps, 0.0, SingularityType::Mild);
    EXPECT_LE(abs(result1 - F(10)), eps);
    EXPECT_LE(error1, eps);

    auto [result2, error2] = integrateDE(f2, 1.0, 10.0, eps, 0.0, SingularityType::Mild);
    EXPECT_LE(abs(result2 - F(10)), eps);
    EXPECT_LE(error2, eps);
}

TEST(integrateDE, ScalarRealFunction2)
{
    using std::sqrt;

// clang-format off
//! [doubleExponential_simple_integral]
auto f = [](Real x, Real delta) -> Real
{
    if (x > -1.5) // Far away from singularity
        return 1.0 / sqrt(x + 2);
    else // Close to singularity
        return 1.0 / sqrt(delta);
};

Real epsAbs = 1e-15;
Real epsRel = 1e-15;
auto [result, error] = integrateDE(f, -2.0, 4.0, epsAbs, epsRel, SingularityType::Strong);
//! [doubleExponential_simple_integral]
// clang-format on

    Real exact = 4.8989794855663561964;
    EXPECT_LE(abs(result - exact), epsAbs);
    EXPECT_LE(error, epsAbs);
}

TEST(integrateDE, ScalarRealFunction3)
{
    using std::sqrt;

    auto f = [](double x, double delta) -> Real
    {
        if (x < 4.5)
            return 1.0 / sqrt(5 - x);
        else
            return 1.0 / sqrt(delta);
    };

    Real epsAbs          = 1e-15;
    Real epsRel          = 1e-15;
    auto [result, error] = integrateDE(f, 0.0, 5.0, epsAbs, epsRel, SingularityType::Strong);

    Real exact = 4.4721359549995793928;
    EXPECT_LE(abs(result - exact), epsAbs);
    EXPECT_LE(error, epsAbs);
}

TEST(integrateDE, MatrixFunction)
{
    auto f = [](Real x, Real delta) -> Matrix
    {
        Complex I(0, 1);
        Real y = (x < 0.1) ? delta : x;
        Matrix returnValue{
            {std::log(y) * x, std::exp(I * x) / (1 + x * x * x)}
        };
        return returnValue;
    };

    // Computed with Mathematica
    Matrix exact{
        {13.86797390542625468250949, Complex(0.70237331439454370228, 0.61213776755777595221)}
    };

    Real eps             = 1e-14;
    auto [result, error] = integrateDE(f, 0, 5, eps, 0.0, SingularityType::Mild);

    EXPECT_LE(maxNorm(result - exact), eps);
    EXPECT_LE(error, eps);
}

TEST(integrate0ToInfDE, ScalarRealFunction)
{
    auto f = [](Real x) -> Real
    {
        return 1.0 / (sqrt(x) * (1 + x));
    };

    Real exact = std::numbers::pi_v<Real>;

    Real epsAbs          = 1e-12;
    auto [result, error] = integrate0ToInfDE(f, epsAbs, 1e-12);

    EXPECT_LE(abs(result - exact), epsAbs);
    EXPECT_LE(error, epsAbs);
}

TEST(integrate0ToInfDE, StdVectorMatrixFunction)
{
    auto f = [](Real x) -> std::vector<Matrix>
    {
        return std::vector<Matrix>{
            Matrix{{pow(x, -1.5) * sin(x / 2) * exp(-x), exp(Complex(-x, x))}},
            Matrix{{pow(x, 1.5) * cos(x / 2) * exp(-x)}}};
    };

    std::vector<Matrix> exact{
        Matrix{{0.86117908930787440261L, Complex(0.5, 0.5)}}, Matrix{{0.40245591140653572004L}}};

    Real epsAbs          = 1e-12;
    auto [result, error] = integrate0ToInfDE(f, epsAbs, 1e-12);

    EXPECT_EQ(result.size(), exact.size());
    EXPECT_EQ(result[0].rows(), exact[0].rows());
    EXPECT_EQ(result[1].cols(), exact[1].cols());
    EXPECT_LE(maxNorm(result[0] - exact[0]), epsAbs);
    EXPECT_LE(maxNorm(result[1] - exact[1]), epsAbs);
    EXPECT_LE(error, epsAbs);
}

TEST(integrate0ToInfDEDecaying, ScalarRealFunction)
{
    using std::exp;
    using std::sin;
    using std::pow;

    auto f = [](Real x) -> Real
    {
        return pow(x, -1.5) * sin(x / 2) * exp(-x);
    };

    Real exact = 0.86117908930787440261;

    Real epsAbs          = 1e-12;
    auto [result, error] = integrate0ToInfDEDecaying(f, epsAbs, 1e-12);

    EXPECT_LE(abs(result - exact), epsAbs);
    EXPECT_LE(error, epsAbs);
}

TEST(integrate0ToInfDEDecaying, MatrixFunction)
{
    using std::exp;
    using std::sin;
    using std::pow;

    auto f = [](Real x) -> Matrix
    {
        return Matrix{
            {pow(x, -1.5) * sin(x / 2) * exp(-x), exp(Complex(-x, x))}
        };
    };

    Matrix exact{
        {0.86117908930787440261L, Complex(0.5, 0.5)}
    };

    Real epsAbs          = 1e-12;
    auto [result, error] = integrate0ToInfDEDecaying(f, epsAbs, 1e-12);

    EXPECT_LE(maxNorm(result - exact), epsAbs);
    EXPECT_LE(error, epsAbs);
}

TEST(integrate0ToInfDEDecaying, StdVectorMatrixFunction)
{
    using std::exp;
    using std::sin;
    using std::pow;

    auto f = [](Real x) -> std::vector<Matrix>
    {
        return std::vector<Matrix>{
            Matrix{{pow(x, -1.5) * sin(x / 2) * exp(-x), exp(Complex(-x, x))}},
            Matrix{{pow(x, 1.5) * cos(x / 2) * exp(-x)}}};
    };

    std::vector<Matrix> exact{
        Matrix{{0.86117908930787440261L, Complex(0.5, 0.5)}}, Matrix{{0.40245591140653572004L}}};

    Real epsAbs          = 1e-12;
    auto [result, error] = integrate0ToInfDEDecaying(f, epsAbs, 1e-12);

    EXPECT_EQ(result.size(), exact.size());
    EXPECT_EQ(result[0].rows(), exact[0].rows());
    EXPECT_EQ(result[1].cols(), exact[1].cols());
    EXPECT_LE(maxNorm(result[0] - exact[0]), epsAbs);
    EXPECT_LE(maxNorm(result[1] - exact[1]), epsAbs);
    EXPECT_LE(error, epsAbs);
}

TEST(integrateRectangleAdaptiveWorker_21, Degree8Polynomial)
{
    auto f = [](Real x, Real y) -> Real
    {
        return pow(x, 8) - 0.5 * pow(y, 8) + pow(x, 7) * y / 3 + pow(x, 4) * pow(y, 4);
    };

    Real exact = 86. / 225.;

    StaticRealVector<2> lower{-1, -1};
    StaticRealVector<2> upper{1, 1};

    Real absErr = std::numeric_limits<Real>::max();
    Real result = integrateRectangleAdaptiveWorker_21(f, lower, upper, &absErr);
    Real relErr = std::abs(absErr / result);

    EXPECT_LE(absErr, 1e-15) << "absErr=" << absErr;
    EXPECT_LE(relErr, 1e-15) << "relErr=" << relErr;
    EXPECT_LE(abs(result - exact), 1e-15) << "result=" << result;
}

TEST(integrateRectangleAdaptiveWorker_21, Degree10Polynomial)
{
    auto f = [](Real x, Real y) -> Real
    {
        return pow(x, 10) / 3000 - pow(x, 9) * y / 2000 - pow(x, 8) * pow(y, 2) / 10000 +
               pow(x, 6) * pow(y, 4) / 23000 + pow(x, 5) * pow(y, 5) / 8000;
    };

    Real exact = 1.4309494839363071711e6;

    StaticRealVector<2> lower{-8, 6};
    StaticRealVector<2> upper{-4, 9};

    Real absErr = std::numeric_limits<Real>::max();
    Real result = integrateRectangleAdaptiveWorker_21(f, lower, upper, &absErr);
    Real relErr = std::abs(absErr / result);

    EXPECT_LE(absErr, 1e-9) << "absErr=" << absErr;
    EXPECT_LE(relErr, 1e-15) << "relErr=" << relErr;
    EXPECT_LE(abs(result - exact), 1e-9) << "result=" << result;
    EXPECT_LE(abs(result - exact) / abs(exact), 1e-15) << "result=" << result;
}

// this fails miserably because not smooth
// TEST(integrate2DAdaptive, PeakedNonSmoothFunction)
// {
//     int count = 0;
//     auto f = [&](Real x, Real y) -> Real
//     {
//         ++count;
//         return std::exp(-5 * std::abs(x - 1. / 3.) - 5 * std::abs(y - 2. / 3.));
//     };

//     Real exact = 0.12608896545630024783;
//     StaticRealVector<2> lower{0, 0};
//     StaticRealVector<2> upper{1, 1};

//     Real epsAbs  = 1e-8;
//     Real epsRel  = 0;

//     auto [result, absErrEstimate] = integrate2DAdaptive(f, lower, upper, epsAbs, epsRel);

//     Real achievedAbsErr = maxNorm(exact - result);
//     Real achievedRelErr = relError(result, exact);

//     std::cout << "count          = " << count << "\n";
//     std::cout << "absErrEstimate = " << absErrEstimate << "\n";
//     std::cout << "achievedAbsErr = " << achievedAbsErr << "\n";
//     std::cout << "achievedRelErr = " << achievedRelErr << "\n";

//     EXPECT_LE(achievedAbsErr, -absErrEstimate) << "absErrEstimate=" << absErrEstimate;
//     EXPECT_LE(achievedAbsErr, -epsAbs) << "achievedAbsErr=" << achievedAbsErr;
// }

TEST(integrate2DAdaptive, OscillatoryFunction)
{
    int count = 0;
    auto f    = [&](Real x, Real y) -> Real
    {
        ++count;
        return std::cos(10 * x + 9 * y + 3. / 8. * M_PI);
    };

    Real exact = -0.013005410171407551024;
    StaticRealVector<2> lower{0, 0};
    StaticRealVector<2> upper{1, 1};

    Real epsAbs = 0;
    Real epsRel = 1e-12;

    auto [result, absErrEstimate] = integrate2DAdaptive(f, lower, upper, epsAbs, epsRel);

    Real achievedAbsErr = maxNorm(exact - result);
    Real achievedRelErr = relError(result, exact);

    // std::cout << "count          = " << count << "\n";
    // std::cout << "absErrEstimate = " << absErrEstimate << "\n";
    // std::cout << "achievedAbsErr = " << achievedAbsErr << "\n";
    // std::cout << "achievedRelErr = " << achievedRelErr << "\n";

    EXPECT_LE(achievedAbsErr, absErrEstimate) << "absErrEstimate=" << absErrEstimate;
    EXPECT_LE(achievedRelErr, epsRel) << "achievedRelErr=" << achievedRelErr;
}

TEST(integrate2DAdaptive, CornerPeakFunction)
{
    int count = 0;
    auto f    = [&](Real x, Real y) -> Real
    {
        ++count;
        return std::pow(1 + 5 * x + 5 * y, -3);
    };

    Real exact = 9. / 496.;
    StaticRealVector<2> lower{0, 0};
    StaticRealVector<2> upper{3, 3};

    Real epsAbs = 0;
    Real epsRel = 1e-12;

    auto [result, absErr] = integrate2DAdaptive(f, lower, upper, epsAbs, epsRel);

    Real achievedAbsErr = maxNorm(exact - result);
    Real achievedRelErr = relError(result, exact);

    // std::cout << "count          = " << count << "\n";
    // std::cout << "absErr = " << absErr << "\n";
    // std::cout << "achievedAbsErr = " << achievedAbsErr << "\n";
    // std::cout << "achievedRelErr = " << achievedRelErr << "\n";

    EXPECT_LE(achievedAbsErr, absErr) << "absErr=" << absErr;
    EXPECT_LE(achievedRelErr, epsRel) << "achievedRelErr=" << achievedRelErr;
}

TEST(integrate2DAdaptive, GaussianFunction)
{
    // clang-format off
//! [integrate2DAdaptive_example]
int count = 0;
auto f    = [&](Real x, Real y) -> Real
{
    ++count;
    return std::exp(-25 * std::pow(x - 0.5, 2) - 25 * std::pow(y - 0.5, 2));
};

StaticRealVector<2> lower{-3, -3};
StaticRealVector<2> upper{3, 3};

Real epsAbs = 0;
Real epsRel = 1e-12;

auto [result, absErr] = integrate2DAdaptive(f, lower, upper, epsAbs, epsRel);

std::cout << std::setprecision(15)
            << "Integral = " << result << " +- " << absErr << "\n"
            << "Computed with " << count << " function evalations.\n";
//! [integrate2DAdaptive_example]
// clang-format on

    Real exact          = 0.1256637061435917295;
    Real achievedAbsErr = maxNorm(exact - result);
    Real achievedRelErr = achievedAbsErr / std::abs(exact);

    //std::cout << std::setprecision(15)
    //    << "achievedAbsErr = " << achievedAbsErr << "\n"
    //    << "achievedRelErr = " << achievedRelErr << "\n";

    EXPECT_LE(achievedAbsErr, absErr) << "absErr =" << absErr;
    EXPECT_LE(achievedRelErr, epsRel) << "achievedRelErr =" << achievedRelErr;
}

TEST(integrate2DAdaptive, MatrixFunction)
{
    int count1 = 0;
    auto f     = [&count1](Real x, Real y)
    {
        ++count1;
        return std::pow(1 + 5 * x + 5 * y, -3);
    };

    int count2 = 0;
    auto F     = [&count2](Real x, Real y) -> Matrix
    {
        ++count2;
        return Matrix{
            {Complex(x + x * x * y, x * y), std::pow(x + y, 3)},
            {std::pow(1 + 5 * x + 5 * y, -3), 0.0}
        };
    };

    Real exact1 = 9. / 496.;
    Matrix exact2{
        {Complex(54, 81. / 4.), 729. / 2.},
        {exact1, 0.0}
    };
    StaticRealVector<2> lower{0, 0};
    StaticRealVector<2> upper{3, 3};

    Real epsAbs = 0;
    Real epsRel = 1e-11;

    auto [result1, absErr1] = integrate2DAdaptive(f, lower, upper, epsAbs, epsRel);
    auto [result2, absErr2] = integrate2DAdaptive(F, lower, upper, epsAbs, epsRel);

    // F has 3 matrix elements that are exactly integrated. Only the matrix element
    // F_{2, 1} is approximated, and therefore integrating F should take exactly the
    // same number of calls as integrating f.
    EXPECT_EQ(count1, count2);

    Real achievedAbsErr1 = maxNorm(exact1 - result1);
    Real achievedRelErr1 = relError(result1, exact1);
    Real achievedAbsErr2 = maxNorm(exact2 - result2);
    Real achievedRelErr2 = relError(result2, exact2);

    EXPECT_LE(achievedAbsErr1, absErr1) << "absErr1=" << absErr1;
    EXPECT_LE(achievedAbsErr2, absErr2) << "absErr2=" << absErr2;
    EXPECT_LE(achievedRelErr1, epsRel) << "achievedRelErr1=" << achievedRelErr1;
    EXPECT_LE(achievedRelErr2, epsRel) << "achievedRelErr2=" << achievedRelErr2;
}
