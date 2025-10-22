// SPDX-FileCopyrightText: 2025 Goethe Universitaet Frankfurt
// SPDX-License-Identifier: EUPL-1.2
// SPDX-FileContributor: Arne Naegel 
// SPDX-FileType: SOURCE

#include <boost/test/unit_test.hpp> 		// Boost test framework
//g++ unit_test.cpp -I/Users/anaegel/Software/third-party/nlohmann/json/include/ -I/Users/anaegel/Software/ug4-git/externals/BoostForUG4/ -std=c++11 -I/Users/anaegel/adolc_base/include -L /Users/anaegel/adolc_base/lib64/ -ladolc


#include "../van_genuchten.h"      			// Plugin dependencies

// Define common variables etc
struct RichardsFixtures
{
	RichardsFixtures()
	: SiltLoamParams ({ 0.396, 0.131, 0.423, 2.06, 4.96e-2}),
	  model(SiltLoamParams)
	{};

	~RichardsFixtures() {};

	ug::Richards::VanGenuchtenParameters SiltLoamParams;
	ug::Richards::VanGenuchtenModel model;
};


// Boost start here
BOOST_AUTO_TEST_SUITE(RichardsTests_Suite_No1)

    BOOST_AUTO_TEST_CASE(RichardsTests_simple_checks)
    {
        BOOST_WARN(false);
        BOOST_CHECK(false);
        BOOST_MESSAGE("If the next test fails, this test case will abort");
        BOOST_REQUIRE(false);
    }


    BOOST_FIXTURE_TEST_CASE(PluginTests_simple_fixture, RichardsFixtures)
    {
    	RichardsFixtures fixtures;
    	ug::Richards::VanGenuchtenModel &model =  fixtures.model;
        BOOST_WARN_EQUAL(1,2);

        BOOST_WARN_GT(1,2);

        BOOST_CHECK_LE(0.0, model.Saturation(-1.0));
    }

BOOST_AUTO_TEST_SUITE_END()
