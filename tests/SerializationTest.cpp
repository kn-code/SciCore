//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#include <fstream>

#include <gtest/gtest.h>

#include <SciCore/Random.h>
#include <SciCore/Serialization.h>

using namespace SciCore;

TEST(serialization, RealVector)
{
    uint64_t seed = 1234;
    Xoshiro256 rng(seed);
    std::uniform_real_distribution<> dis(-1000.0, 1000.0);

    std::string jsonFilename = "serialization_test_1.json";
    std::remove(jsonFilename.c_str());
    RealVector data = randomVector<RealVector>(1000, dis, rng);

    // Write
    {
        std::ofstream os(jsonFilename);
        cereal::JSONOutputArchive archive(os);
        archive(cereal::make_nvp("data", data));
    }

    // Read
    RealVector readData;
    {
        std::ifstream is(jsonFilename);
        cereal::JSONInputArchive archive(is);
        archive(cereal::make_nvp("data", readData));
    }

    EXPECT_EQ(data, readData);
}

TEST(serialization, Vector)
{
    uint64_t seed = 1234;
    Xoshiro256 rng(seed);
    std::uniform_real_distribution<> dis(-1000.0, 1000.0);

    std::string jsonFilename = "serialization_test_2.json";
    std::remove(jsonFilename.c_str());
    Vector data = randomVector<Vector>(1000, dis, rng);

    // Write
    {
        std::ofstream os(jsonFilename);
        cereal::JSONOutputArchive archive(os);
        archive(cereal::make_nvp("data", data));
    }

    // Read
    Vector readData;
    {
        std::ifstream is(jsonFilename);
        cereal::JSONInputArchive archive(is);
        archive(cereal::make_nvp("data", readData));
    }

    EXPECT_EQ(data, readData);
}

TEST(serialization, RowVector)
{
    uint64_t seed = 1234;
    Xoshiro256 rng(seed);
    std::uniform_real_distribution<> dis(-1000.0, 1000.0);

    std::string jsonFilename = "serialization_test_3.json";
    std::remove(jsonFilename.c_str());
    RowVector data = randomVector<RowVector>(1000, dis, rng);

    // Write
    {
        std::ofstream os(jsonFilename);
        cereal::JSONOutputArchive archive(os);
        archive(cereal::make_nvp("data", data));
    }

    // Read
    RowVector readData;
    {
        std::ifstream is(jsonFilename);
        cereal::JSONInputArchive archive(is);
        archive(cereal::make_nvp("data", readData));
    }

    EXPECT_EQ(data, readData);
}

TEST(serialization, RealMatrix)
{
    uint64_t seed = 1234;
    Xoshiro256 rng(seed);
    std::uniform_real_distribution<> dis(-1000.0, 1000.0);

    std::string jsonFilename = "serialization_test_4.json";
    std::remove(jsonFilename.c_str());
    RealMatrix data = randomMatrix<RealMatrix>(10, 15, dis, rng);

    // Write
    {
        std::ofstream os(jsonFilename);
        cereal::JSONOutputArchive archive(os);
        archive(cereal::make_nvp("data", data));
    }

    // Read
    RealMatrix readData;
    {
        std::ifstream is(jsonFilename);
        cereal::JSONInputArchive archive(is);
        archive(cereal::make_nvp("data", readData));

    }

    EXPECT_EQ(data, readData);
}

TEST(serialization, Matrix)
{
    uint64_t seed = 1234;
    Xoshiro256 rng(seed);
    std::uniform_real_distribution<> dis(-1000.0, 1000.0);

    std::string jsonFilename = "serialization_test_5.json";
    std::remove(jsonFilename.c_str());
    Matrix data = randomMatrix<Matrix>(3, 7, dis, rng);

    // Write
    {
        std::ofstream os(jsonFilename);
        cereal::JSONOutputArchive archive(os);
        archive(cereal::make_nvp("data", data));
    }

    // Read
    Matrix readData;
    {
        std::ifstream is(jsonFilename);
        cereal::JSONInputArchive archive(is);
        archive(cereal::make_nvp("data", readData));
    }

    EXPECT_EQ(data, readData);
}

TEST(serialization, EmptyObjects)
{
    std::string jsonFilename = "serialization_test_6.json";
    std::remove(jsonFilename.c_str());
    RealVector emptyRealVector;
    Vector emptyVector;
    RealRowVector emptyRealRowVector;
    RealMatrix emptyRealMatrix;
    Matrix emptyMatrix;

    // Write
    {
        std::ofstream os(jsonFilename);
        cereal::JSONOutputArchive archive(os);
        archive(cereal::make_nvp("emptyRealVector", emptyRealVector));
        archive(cereal::make_nvp("emptyVector", emptyVector));
        archive(cereal::make_nvp("emptyRealRowVector", emptyRealRowVector));
        archive(cereal::make_nvp("emptyRealMatrix", emptyRealMatrix));
        archive(cereal::make_nvp("emptyMatrix", emptyMatrix));
    }

    // Read
    RealVector readEmptyRealVector{
        {1, 2, 3}
    };
    Vector readEmptyVector{
        {1.0, 2.0, Complex(0, 3.0)}
    };
    RealRowVector readEmptyRealRowVector{
        {1, 2, 3}
    };
    RealMatrix readEmptyRealMatrix{
        {1, 2, 3}
    };
    Matrix readEmptyMatrix{
        {1.0, 2.0, Complex(0, 3.0)}
    };

    {
        std::ifstream is(jsonFilename);
        cereal::JSONInputArchive archive(is);
        archive(cereal::make_nvp("emptyRealVector", readEmptyRealVector));
        archive(cereal::make_nvp("emptyVector", readEmptyVector));
        archive(cereal::make_nvp("emptyRealRowVector", readEmptyRealRowVector));
        archive(cereal::make_nvp("emptyRealMatrix", readEmptyRealMatrix));
        archive(cereal::make_nvp("emptyMatrix", readEmptyMatrix));
    }

    EXPECT_EQ(emptyRealVector, readEmptyRealVector);
    EXPECT_EQ(emptyVector, readEmptyVector);
    EXPECT_EQ(emptyRealRowVector, readEmptyRealRowVector);
    EXPECT_EQ(emptyRealMatrix, readEmptyRealMatrix);
    EXPECT_EQ(emptyMatrix, readEmptyMatrix);
}
