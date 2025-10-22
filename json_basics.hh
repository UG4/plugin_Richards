// JSON lib.

#pragma once

#ifdef UG_JSON
#include <nlohmann/json.hpp>
#endif

namespace ug {

// #define UG4_WITH_JSON
#ifdef UG_JSON
typedef nlohmann::json JSONType;
using JSONPointer = nlohmann::json::json_pointer;
#endif

#ifdef UG_JSON

#ifndef UG_HAS_JSON_BASICS
// These structs are defined in "bindings/json/json_basics.hh"
template <class T>
struct is_json_constructible
{
	const static bool value = std::is_default_constructible<T>::value;
};

/*template <class T>
inline constexpr bool is_json_constructible_v = is_json_constructible<T>::value;
*/
#endif


//! This constructs an object from JSON.
template <typename P>
SmartPtr<P> JSONSerializer(nlohmann::json j)
{
	UG_LOG("JSONSerializer will become deprecated. Directly use JSONBuilder instead!")
	UG_COND_THROW(! is_json_constructible<P>::value, "ERROR: Type is not constructible!")
	SmartPtr<P> data = make_sp(new P());
	j.get_to<P>(*data);
	return data;
};

//! This constructs an object from a string.
template <typename P>
SmartPtr<P> JSONSerializer(const char *jstring)
{
	nlohmann::json j = nlohmann::json::parse(jstring);
	return JSONSerializer<P>(j);

};

#endif //UG_JSON


} // namespace ug
