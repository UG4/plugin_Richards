################################################################################
# SAMPLE
#
# This script shows how easy it is to add a build script for your plugins.
# Simply specify the plugins name and a list of your source files.
################################################################################

cmake_minimum_required(VERSION 2.6)

# Name of your plugin and sources.
set(pluginName	Richards)

project(UG_TEST_${pluginName})

# Include the definitions and dependencies for ug-plugins.
# include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)


# Library.
set(SOURCES		van_genuchten.cpp
				richards_plugin.cpp)

# Dependencies: JSON, ADOL-C
include_directories(AFTER ${UG_ROOT_PATH}/externals/JSONForUG4/json-cxx/single_include)  # JSON

# include_directories(AFTER /Users/anaegel/adolc_base/include)
link_directories(${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})   # ADOL-C

# This plugin requires C++ 11.
set (CMAKE_CXX_STANDARD 11)

# Avoid linking error: undefined reference to boost::system::system_category()
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)

# Dependencies.
add_dependencies(${pluginName} adolc)
