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
}
