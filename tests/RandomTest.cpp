//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include <gtest/gtest.h>

#include <SciCore/Random.h>

using namespace SciCore;

TEST(Xoshiro256, Test)
{
// clang-format off
//! [random_example]
// Create seeded random number generator
uint64_t seed = 1234;
Xoshiro256 rng(seed);

// Uniform distribution between 0 and 100
std::uniform_real_distribution<> dis(0.0, 100.0);

// Create a single random real number
Real randomRealScalar = dis(rng);

// Create a RealVector with 100000 entries uniformly random in [0,100]
auto randomRealVector = randomVector<RealVector>(100'000, dis, rng);
//! [random_example]
// clang-format on

        EXPECT_LE(randomRealScalar, 99);
        EXPECT_GE(randomRealScalar, 0);
        EXPECT_GE(randomRealVector.sum() / randomRealVector.size(), 49.0);
        EXPECT_LE(randomRealVector.sum() / randomRealVector.size(), 51.0);
}
