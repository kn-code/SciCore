//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include "SciCore/Random.h"

#include <limits>

namespace SciCore
{

uint64_t splitMix64(uint64_t x) noexcept
{
    uint64_t z = (x += 0x9e3779b97f4a7c15);
    z          = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z          = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

Xoshiro256::Xoshiro256() noexcept
{
    std::random_device rd;
    s[0] = splitMix64(rd());
    s[1] = splitMix64(s[0]);
    s[2] = splitMix64(s[1]);
    s[3] = splitMix64(s[2]);
}

Xoshiro256::Xoshiro256(uint64_t seed) noexcept
{
    s[0] = splitMix64(seed);
    s[1] = splitMix64(s[0]);
    s[2] = splitMix64(s[1]);
    s[3] = splitMix64(s[2]);
}

uint64_t Xoshiro256::operator()() noexcept
{
    uint64_t result_starstar = Xoshiro256::rotl(s[1] * 5, 7) * 9;
    uint64_t t               = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;

    s[3] = Xoshiro256::rotl(s[3], 45);

    return result_starstar;
}

void Xoshiro256::jump() noexcept
{
    constexpr uint64_t JUMP[] = {0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c};

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    for (int i = 0; i < static_cast<int>(sizeof JUMP / sizeof *JUMP); ++i)
    {
        for (int b = 0; b < 64; b++)
        {
            if (JUMP[i] & uint64_t(1) << b)
            {
                s0 ^= s[0];
                s1 ^= s[1];
                s2 ^= s[2];
                s3 ^= s[3];
            }
            this->operator()();
        }
    }

    s[0] = s0;
    s[1] = s1;
    s[2] = s2;
    s[3] = s3;
}

void Xoshiro256::longJump() noexcept
{
    constexpr uint64_t LONG_JUMP[] = {0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635};

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    for (int i = 0; i < static_cast<int>(sizeof LONG_JUMP / sizeof *LONG_JUMP); i++)
    {
        for (int b = 0; b < 64; b++)
        {
            if (LONG_JUMP[i] & uint64_t(1) << b)
            {
                s0 ^= s[0];
                s1 ^= s[1];
                s2 ^= s[2];
                s3 ^= s[3];
            }
            this->operator()();
        }
    }

    s[0] = s0;
    s[1] = s1;
    s[2] = s2;
    s[3] = s3;
}

} // namespace SciCore
