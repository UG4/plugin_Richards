/*
 * unit_test.cpp
 *
 *  Created on: 21.12.2019
 *      Author: anaegel
 */
//g++ unit_test.cpp -I/Users/anaegel/Software/third-party/nlohmann/json/include/ -I/Users/anaegel/Software/ug4-git/externals/BoostForUG4/ -std=c++11 -I/Users/anaegel/adolc_base/include -L /Users/anaegel/adolc_base/lib64/ -ladolc

#include <iostream>
#include <sstream>

#include "van_genuchten.h"

using namespace ug::Richards;

double VanGenuchtenModel::one = 1.0;

namespace ug {
namespace Richards {

/// JSON serialize.
void to_json(JSONType& j, const VanGenuchtenParameters& p) {
	j = JSONType{{"thetaS", p.thetaS}, {"thetaR", p.thetaR}, {"alpha", p.alpha}, {"n", p.n}, {"Ksat", p.Ksat}};
}


/// JSON de-serialize.
void from_json(const JSONType& jobject, VanGenuchtenParameters& p) {

	jobject.at("alpha").get_to(p.alpha);
	jobject.at("n").get_to(p.n);

     // Optional parameters
	{
		p.Ksat = 1.0;
		auto it = jobject.find("Ksat");
		if (it != jobject.end()) {it->get_to(p.Ksat);}
	}

	{
		p.thetaS = 1.0;
		auto it = jobject.find("thetaS");
		if (it != jobject.end()) {it->get_to(p.thetaS);}
	}

	{
		p.thetaR = 0.0;
		auto it = jobject.find("thetaR");
		if (it != jobject.end()) {it->get_to(p.thetaR);}
	}
}


void CreateJSONMap(const JSONType &array, std::map<std::string, JSONType> &map)
{

	for (JSONType::const_iterator it = array.cbegin(); it!=array.cend(); ++it)
	{
		if (!it->is_array()) continue;

		JSONType jelem = *it;
		jelem.erase("uid");
		map.insert( std::pair<std::string, JSONType>(it->at("uid"), *it) );
	}

	for (std::map<std::string, JSONType>::iterator it=map.begin(); it!=map.end(); ++it)
	    std::cout << it->first << " => " << it->second << '\n';
}


std::string VanGenuchtenModel::config_string() const {
	nlohmann::json j2 = m_param;
	std::stringstream ss;
	ss << j2;
	return ss.str();
}

// Factory function
SmartPtr<VanGenuchtenModel> CreateVanGenuchtenModel(const char *jstring) {

	SmartPtr<VanGenuchtenModel> inst = SPNULL;
	try
	{
			nlohmann::json j = nlohmann::json::parse(jstring);
			VanGenuchtenParameters p = j.get<VanGenuchtenParameters>();
			inst = make_sp(new VanGenuchtenModel(p));
	}
	catch (...)
	{
			std::cout << "Construction failed!" << std::endl;
	}
	return inst;
};


/*SmartPtr<VanGenuchtenModel> VanGenuchtenModelFactory::create_default()
{
	// VanGenuchtenParameters SiltLoamParams = { 0.396, 0.131, 0.423, 2.06, 4.96e-2};
	VanGenuchtenParameters BeitNetofaClayParams = { 0.446, 0.0, 0.152, 1.17, 8.2e-4};
	SmartPtr<VanGenuchtenModel> inst = make_sp(new VanGenuchtenModel(BeitNetofaClayParams));
	return inst;
};
*/

SmartPtr<VanGenuchtenModel> VanGenuchtenModelFactory::create(const char *jstring)
{
	return CreateVanGenuchtenModel(jstring);
};

}
}
