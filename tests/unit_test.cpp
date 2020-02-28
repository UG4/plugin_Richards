/*
 * unit_test.cpp
 *
 *  Created on: 21.12.2019
 *      Author: anaegel
 */
//g++ unit_test.cpp -I/Users/anaegel/Software/third-party/nlohmann/json/include/ -I/Users/anaegel/Software/ug4-git/externals/BoostForUG4/ -std=c++11 -I/Users/anaegel/adolc_base/include -L /Users/anaegel/adolc_base/lib64/ -ladolc

#include <iostream>
#include "../van_genuchten.h"

using namespace ug::Richards;


int main(int argc, char **argv)
{
	/*VanGenuchtenParameters SiltLoamParams = { 0.396, 0.131, 0.423, 2.06, 4.96e-2};

	nlohmann::json j = SiltLoamParams;
	std::cout << j << std::endl;


	VanGenuchtenModel model(SiltLoamParams);

	std::cout << "Saturations:"<< std::endl;
	std::cout << model.Saturation(-1.0) << std::endl;
	std::cout << model.Saturation(-0.5) << std::endl;
	std::cout << model.Saturation(0.0) << std::endl;
	std::cout << model.Saturation(0.5) << std::endl;
	std::cout << model.Saturation(1.0)<< std::endl;
	std::cout << model.Saturation(2.0)<< std::endl;
	std::cout << model.Saturation(3.0)<< std::endl;
	std::cout << model.Saturation(5.0)<< std::endl;

	std::cout << "DSaturations:"<< std::endl;
	std::cout << model.dSaturation_dH(-1.0) << std::endl;
	std::cout << model.dSaturation_dH(-0.5) << std::endl;
	std::cout << model.dSaturation_dH(0.0) << std::endl;
	std::cout << model.dSaturation_dH(0.5) << std::endl;
	std::cout << model.dSaturation_dH(1.0) << std::endl;
	std::cout << model.dSaturation_dH(2.0) << std::endl;
	std::cout << model.dSaturation_dH(3.0) << std::endl;
	std::cout << model.dSaturation_dH(5.0) << std::endl;


	std::cout << "Conductivities:"<< std::endl;
	std::cout << model.Conductivity(-1.0) << std::endl;
	std::cout << model.Conductivity(-0.5) << std::endl;
	std::cout << model.Conductivity(0.0) << std::endl;
	std::cout << model.Conductivity(0.5) << std::endl;
	std::cout << model.Conductivity(1.0) << std::endl;
	std::cout << model.Conductivity(2.0) << std::endl;
	std::cout << model.Conductivity(3.0) << std::endl;
	std::cout << model.Conductivity(5.0) << std::endl;

*/


	return 0;
}

