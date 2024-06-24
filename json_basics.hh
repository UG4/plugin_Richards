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


template <class T>
struct is_json_constructible
{
	const static bool value = std::is_default_constructible<T>::value;
};

template <class T>
inline constexpr bool is_json_constructible_v = is_json_constructible<T>::value;



//! This constructs an object from JSON.
template <typename P>
SmartPtr<P> JSONSerializer(nlohmann::json j)
{
	UG_COND_THROW(! is_json_constructible_v<P>, "ERROR: Type is not constructible!")
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

#endif


} // namespace ug
