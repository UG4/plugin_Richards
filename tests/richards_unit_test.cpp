
// SPDX-FileCopyrightText: 2025 Goethe Universitaet Frankfurt
// SPDX-License-Identifier: EUPL-1.2
// SPDX-FileContributor: Arne Naegel 
// SPDX-FileType: SOURCE

#include <gtest/gtest.h> // Google test framework.
// Build with the CMake configuration documented in README.md.

#include "../van_genuchten.h"

namespace {

struct RichardsFixtures : public ::testing::Test
{
    RichardsFixtures()
        : SiltLoamParams{0.396, 0.131, 0.423, 2.06, 4.96e-2},
          model(SiltLoamParams)
    {}

    ug::Richards::VanGenuchtenParameters SiltLoamParams;
    ug::Richards::VanGenuchtenModel model;
};

TEST(RichardsTests, DefaultMIsDerivedFromN)
{
    ug::Richards::VanGenuchtenParameters params;

    EXPECT_DOUBLE_EQ(params.n, 2.0);
    EXPECT_NEAR(params.m, 1.0 - (1.0 / params.n), 1e-12);
}

TEST_F(RichardsFixtures, SaturationForNegativeHeadIsFinite)
{
    EXPECT_GE(model.Saturation(-1.0), 0.0);
    EXPECT_LE(model.Saturation(-1.0), 1.0);
}

}  // namespace
