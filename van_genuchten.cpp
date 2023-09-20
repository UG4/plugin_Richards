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


namespace ug {
namespace Richards {

//////////////////////////////////////////////
// VanGenuchtenParameters
//////////////////////////////////////////////

#ifdef UG_JSON
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

	{
		p.m =  1.0 - (1.0/p.n);
		auto it = jobject.find("m");
		if (it != jobject.end()) {it->get_to(p.m);}
	}
}


/// JSON serialize.
void to_json(JSONType& j, const HaverkampParameters& p) {
	j = JSONType{{"thetaS", p.thetaS}, {"thetaR", p.thetaR}, {"alpha", p.alpha}, {"n", p.n}, {"beta", p.beta}, {"m", p.m}, {"Ksat", p.Ksat}};
}


/// JSON de-serialize.
void from_json(const JSONType& jobject, HaverkampParameters& p) {

	jobject.at("alpha").get_to(p.alpha);
	jobject.at("n").get_to(p.n);

	jobject.at("beta").get_to(p.beta);
	jobject.at("m").get_to(p.m);

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

//////////////////////////////////////////////
// VanGenuchtenModel
//////////////////////////////////////////////




//! Create a map uid -> json
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






/*SmartPtr<VanGenuchtenModel> VanGenuchtenModelFactory::create_default()
{
	// VanGenuchtenParameters SiltLoamParams = { 0.396, 0.131, 0.423, 2.06, 4.96e-2};
	VanGenuchtenParameters BeitNetofaClayParams = { 0.446, 0.0, 0.152, 1.17, 8.2e-4};
	SmartPtr<VanGenuchtenModel> inst = make_sp(new VanGenuchtenModel(BeitNetofaClayParams));
	return inst;
};
*/

//! Factory functions.
SmartPtr<VanGenuchtenModel> CreateVanGenuchtenModel(const char *json)
{ return CreateModel<VanGenuchtenModel>(json); }


SmartPtr<VanGenuchtenModel> RichardsModelFactory::create_van_genuchten(const char *jstring)
{ return CreateVanGenuchtenModel(jstring); }

SmartPtr<HaverkampModel> RichardsModelFactory::create_haverkamp(const char *jstring)
{ return CreateModel<HaverkampModel>(jstring); }

#endif

}
}
