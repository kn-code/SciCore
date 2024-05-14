//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include <gtest/gtest.h>

#include <SciCore/IDECheb.h>

using namespace SciCore;

struct JCModel
{
    Real gamma   = 1.0;
    Real Gamma   = 13.0;
    Real epsilon = 20.0;

    Matrix minusIL() const
    {
        Complex I(0, 1);
        return Matrix{
            {0, 0,            0,           0},
            {0, 0,            0,           0},
            {0, 0, -I * epsilon,           0},
            {0, 0,            0, I * epsilon}
        };
    }

    Matrix minusIK(Real t) const
    {
        Complex gammaK = std::sqrt(Complex(gamma * gamma - 4 * gamma * Gamma, 0));
        Complex u      = gamma * Gamma * std::exp(-1.5 * gamma * t) *
                    (std::cosh(gammaK * t / 2.0) + Complex(gamma) / gammaK * std::sinh(gammaK * t / 2.0));

        return Matrix{
            {0, u, 0, 0},
            {0, -u, 0, 0},
            {0, 0, -0.5 * std::exp(-Complex(gamma, epsilon) * t) * gamma * Gamma, 0},
            {0, 0, 0, -0.5 * std::exp(-Complex(gamma, -epsilon) * t) * gamma * Gamma}
        };
    }

    Complex compute_g(Real t)
    {
        const Complex gammaPrime = std::sqrt(Complex(gamma * gamma - 2.0 * gamma * Gamma, 0));
        return (std::exp(-gamma * t / 2.0) *
                (std::cosh(gammaPrime * t / 2.0) + gamma / gammaPrime * std::sinh(gammaPrime * t / 2.0)))
            .real();
    };

    Matrix propagator(Real t)
    {
        using std::exp;

        Complex g = compute_g(t);

        const Complex I(0, 1);
        return Matrix{
            {1.0, 1.0 - g * g,                       0.0,                       0.0},
            {0.0,       g * g,                       0.0,                       0.0},
            {0.0,         0.0, exp(-I * epsilon * t) * g,                       0.0},
            {0.0,         0.0,                       0.0, exp(+I * epsilon * t) * g}
        };
    };
};

TEST(convolution, RealFunctions)
{
    auto f = [](Real x) -> Real
    {
        return exp(-x);
    };

    auto g = [](Real t) -> Real
    {
        return sin(t);
    };

    // Exact result of convolution
    auto conv = [](Real x) -> Real
    {
        return (-cos(x) + exp(2 - x) * (cos(2) - sin(2)) + sin(x)) / 2.;
    };

    Real t0   = 2;
    Real tMax = 5;
    Cheb cf(f, 0, tMax - t0, 128);
    Cheb cg(g, t0, tMax, 128);

    Cheb computedConv = convolution<Real>(cf, cg);

    RealVector tTest = RealVector::LinSpaced(100, t0, tMax);
    for (Real t : tTest)
    {
        EXPECT_LT(maxNorm(computedConv(t) - conv(t)), 1e-15) << "t =" << t;
    }
}

TEST(convolution, MatrixAndVectorFunction)
{
    Real t0   = 2;
    Real tMax = 5;

    auto f = [](Real x) -> RealMatrix
    {
        return RealMatrix{
            {exp(-x), sin(x)},
            {      0, cos(x)}
        };
    };

    auto g = [](Real t) -> RealVector
    {
        return RealVector{
            {t, t * t}
        };
    };

    auto conv = [](Real t) -> RealVector
    {
        return RealVector{
            {-3 - exp(2 - t) + t + pow(t, 2) - 2 * cos(2 - t) + 4 * sin(2 - t),
             -2 * (-t + 2 * cos(2 - t) + sin(2 - t))}
        };
    };

    Cheb cf(f, 0, tMax - t0, 128);
    Cheb cg(g, t0, tMax, 128);

    Cheb computedConv = convolution<StaticRealVector<2>>(cf, cg);

    RealVector tTest = RealVector::LinSpaced(100, t0, tMax);
    for (Real t : tTest)
    {
        EXPECT_LT(maxNorm(computedConv(t) - conv(t)), 1e-13) << "t =" << t;
    }
}

TEST(solveIdeCheb, simpleRealHomogeneousIDE)
{
    std::vector<int> nMinChebs{0, 64, 128, 256};
    for (size_t i = 0; i < nMinChebs.size(); ++i)
    {
        int nMinCheb = nMinChebs[i];
        Real t0      = 3;
        Real tMax    = 10;

        auto k = [](Real) -> Real
        {
            return -1;
        };
        Cheb ck(k, 0, tMax - t0, 128);

        auto sol = [](Real t) -> Real
        {
            return -(exp(1.5 - t / 2.) *
                     (-3 * cos((sqrt(3) * (-3 + t)) / 2.) + sqrt(3) * sin((sqrt(3) * (-3 + t)) / 2.))) /
                   3.;
        };

        Real g    = -1.0;
        Real f0   = 1.0;
        Cheb cSol = solveIdeCheb(g, ck, f0, t0, tMax, nMinCheb);
        EXPECT_LE(nMinCheb, (int)cSol.coefficients().size());

        RealVector tTest = RealVector::LinSpaced(100, t0, tMax);
        for (Real t : tTest)
        {
            EXPECT_LT(maxNorm(sol(t) - cSol(t)), 1e-15) << "t =" << t;
        }
    }
}

TEST(solveIdeCheb, simpleRealInhomogeneousIDE)
{
    Real t0   = 3;
    Real tMax = 8;

    auto h = [](Real t) -> Real
    {
        return sin(t + 1);
    };
    Cheb ch(h, t0, tMax, 128);

    auto k = [](Real) -> Real
    {
        return -1;
    };
    Cheb ck(k, 0, tMax - t0, 128);

    auto sol = [](Real t) -> Real
    {
        return -(-2 * sqrt(3) * exp(1.5) * cos((sqrt(3) * (-3 + t)) / 2.) -
                 exp(1.5) * cos((8 + 3 * sqrt(3) - sqrt(3) * t) / 2.) +
                 exp(1.5) * cos((8 - 3 * sqrt(3) + sqrt(3) * t) / 2.) + 2 * exp(1.5) * sin((sqrt(3) * (-3 + t)) / 2.) -
                 2 * sqrt(3) * exp(t / 2.) * sin(1 + t) - 2 * exp(1.5) * sin((8 + 3 * sqrt(3) - sqrt(3) * t) / 2.) +
                 sqrt(3) * exp(1.5) * sin((8 + 3 * sqrt(3) - sqrt(3) * t) / 2.) +
                 2 * exp(1.5) * sin((8 - 3 * sqrt(3) + sqrt(3) * t) / 2.) +
                 sqrt(3) * exp(1.5) * sin((8 - 3 * sqrt(3) + sqrt(3) * t) / 2.)) /
               (2. * sqrt(3) * exp(t / 2.));
    };

    Real g       = -1.0;
    Real f0      = 1.0;
    int nMinCheb = 0;
    Cheb cSol    = solveIdeCheb(g, ck, ch, f0, t0, tMax, nMinCheb);

    RealVector tTest = RealVector::LinSpaced(100, t0, tMax);
    for (Real t : tTest)
    {
        EXPECT_LT(maxNorm(sol(t) - cSol(t)), 1e-15) << "t =" << t;
    }
}

TEST(solveIdeCheb, exactSolutionDissipativeJaynesCummings)
{
    Real t0   = 0;
    Real tMax = 10;

    JCModel jc;

    auto minusIK = [&](Real t) -> Matrix
    {
        return jc.minusIK(t);
    };

    Vector rho0{
        {0.1, 0.1, 0.1, 0.9}
    };
    Cheb cMinusIK(minusIK, 0, tMax - t0, 256);
    int nMinCheb = 0;
    Cheb cSol    = solveIdeCheb(jc.minusIL(), cMinusIK, rho0, t0, tMax, nMinCheb);

    RealVector tTest = RealVector::LinSpaced(100, t0, tMax);
    for (Real t : tTest)
    {
        EXPECT_LT(maxNorm(jc.propagator(t) * rho0 - cSol(t)), 1e-14)
            << "t = " << t << "\nexact = " << (jc.propagator(t) * rho0).transpose()
            << "\ncomputed = " << cSol(t).transpose();
    }
}

TEST(solveIdeCheb, withSections)
{
    Real t0   = 3;
    Real tMax = 10;

    std::vector<RealVector> sections{
        RealVector{{0, tMax - t0}},
        RealVector{{0, 0.13847, tMax - t0}},
        RealVector{{0, (tMax - t0) / 2, tMax - t0}},
        RealVector{{0, 0.25 * (tMax - t0), 0.5 * (tMax - t0), tMax - t0}},
        RealVector{{0, (tMax - t0) / 3, tMax - t0}},
        RealVector{{0, (tMax - t0) / 5, tMax - t0}},
        RealVector{{0, (tMax - t0) / 5, 0.3 * (tMax - t0), 0.7 * (tMax - t0), tMax - t0}},
    };

    auto k = [](Real) -> Real
    {
        return -1;
    };

    auto sol = [](Real t) -> Real
    {
        return -(exp(1.5 - t / 2.) * (-3 * cos((sqrt(3) * (-3 + t)) / 2.) + sqrt(3) * sin((sqrt(3) * (-3 + t)) / 2.))) /
               3.;
    };

    RealVector absoluteErrorGoals{
        {1e-4, 1e-6, 1e-10, 1e-10}
    };
    for (int i = 0; i < sections.size(); ++i)
    {
        RealVector section = sections[i];
        for (Real epsAbs : absoluteErrorGoals)
        {
            Real epsRel = 0;
            Real hMin   = 0;
            ChebAdaptive ck(k, section, epsAbs, epsRel, hMin);

            Real g            = -1.0;
            Real f0           = 1.0;
            int nMinCheb      = 128;
            ChebAdaptive cSol = solveIdeCheb(g, ck, f0, t0, tMax, epsAbs, epsRel, nMinCheb);
            EXPECT_EQ(cSol.lowerLimit(), t0);
            EXPECT_EQ(cSol.upperLimit(), tMax);

            RealVector tTest = RealVector::LinSpaced(200, t0, tMax);
            for (Real t : tTest)
            {
                EXPECT_LT(maxNorm(sol(t) - cSol(t)), epsAbs) << "i = " << i << ", t = " << t;
            }
        }
    }
}

TEST(solveIdeCheb, JaynesCummingsWithSections)
{
    Real t0   = 0;
    Real tMax = 10;

    std::vector<RealVector> sections{
        RealVector{{0, tMax - t0}},
        RealVector{{0, 0.13847, tMax - t0}},
        RealVector{{0, (tMax - t0) / 2, tMax - t0}},
        RealVector{{0, 0.25 * (tMax - t0), 0.5 * (tMax - t0), tMax - t0}},
        RealVector{{0, (tMax - t0) / 3, tMax - t0}},
        RealVector{{0, (tMax - t0) / 5, tMax - t0}},
        RealVector{{0, (tMax - t0) / 5, 0.3 * (tMax - t0), 0.7 * (tMax - t0), tMax - t0}},
    };

    for (int i = 0; i < sections.size(); ++i)
    {
        RealVector section = sections[i];

        Real epsAbs = 1e-14;
        Real epsRel = 0;
        Real hMin   = 0;

        JCModel jc;

        StaticMatrix<4, 4> minusIL = jc.minusIL();
        auto minusIK               = [&](Real t) -> StaticMatrix<4, 4>
        {
            return jc.minusIK(t);
        };

        Vector rho0{
            {0.1, 0.1, 0.1, 0.9}
        };

        ChebAdaptive cMinusIK(minusIK, section, epsAbs, epsRel, hMin);
        int nMinCheb      = 256;
        ChebAdaptive cSol = solveIdeCheb(minusIL, cMinusIK, rho0, t0, tMax, epsAbs, epsRel, nMinCheb);

        RealVector tTest = RealVector::LinSpaced(200, t0, tMax);
        for (Real t : tTest)
        {
            EXPECT_LT(maxNorm(jc.propagator(t) * rho0 - cSol(t)), epsAbs)
                << "section = " << section.transpose() << "\nt = " << t
                << "\nexact = " << (jc.propagator(t) * rho0).transpose() << "\ncomputed = " << cSol(t).transpose();
        }
    }
}

/*TEST(computePropagatorIde, JaynesCummings)
{
    Real t0   = 0;
    Real tMax = 10;

    JCModel jc;

    StaticMatrix<4,4> minusIL = jc.minusIL();
    auto minusIK = [&](Real t) -> StaticMatrix<4,4>
    {
        return jc.minusIK(t);
    };

    int nMinCheb = 0;
    Cheb cMinusIK(minusIK, 0, tMax - t0, 256);
    Cheb cPropagator = computePropagatorIde(minusIL, cMinusIK, t0, tMax, nMinCheb);

    RealVector tTest = RealVector::LinSpaced(200, t0, tMax);
    for (Real t : tTest)
    {
        EXPECT_LT(maxNorm(jc.propagator(t) - cPropagator(t)), 1e-14) << "t =" << t;
    }
}*/

TEST(computePropagatorIde, JaynesCummingsWithSections)
{
    Real t0   = 0;
    Real tMax = 10;

    std::vector<RealVector> sections{
        RealVector{{0, tMax - t0}},
        RealVector{{0, 0.13847, tMax - t0}},
        RealVector{{0, (tMax - t0) / 2, tMax - t0}},
        RealVector{{0, 0.25 * (tMax - t0), 0.5 * (tMax - t0), tMax - t0}},
        RealVector{{0, (tMax - t0) / 3, tMax - t0}},
        RealVector{{0, (tMax - t0) / 5, tMax - t0}},
        RealVector{{0, (tMax - t0) / 5, 0.3 * (tMax - t0), 0.7 * (tMax - t0), tMax - t0}}};

    for (int i = 0; i < sections.size(); ++i)
    {
        RealVector section = sections[i];

        Real epsAbs = 1e-14;
        Real epsRel = 0;
        Real hMin   = 0;

        JCModel jc;

        StaticMatrix<4, 4> minusIL = jc.minusIL();
        auto minusIK               = [&](Real t) -> StaticMatrix<4, 4>
        {
            return jc.minusIK(t);
        };

        ChebAdaptive cMinusIK(minusIK, section, epsAbs, epsRel, hMin);
        int nMinCheb             = 256;
        ChebAdaptive cPropagator = computePropagatorIde(minusIL, cMinusIK, t0, tMax, epsAbs, epsRel, nMinCheb);

        RealVector tTest = RealVector::LinSpaced(200, t0, tMax);
        for (Real t : tTest)
        {
            EXPECT_LT(maxNorm(jc.propagator(t) - cPropagator(t)), epsAbs) << "t = " << t;
        }
    }
}
