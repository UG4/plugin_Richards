# SPDX-FileCopyrightText: 2025 Goethe Universitaet Frankfurt
# SPDX-License-Identifier: EUPL-1.2
# SPDX-FileContributor: Arne Naegel 
# SPDX-FileType: SOURCE

################################################################################
# SAMPLE
#
# This script shows how easy it is to add a build script for your plugins.
# Simply specify the plugins name and a list of your source files.
################################################################################

cmake_minimum_required(VERSION 2.8.12)

# Name of your plugin and sources.
set(pluginName	Richards)

project(UG_PLUGIN_${pluginName})

# Include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)


# Library.
set(SOURCES		van_genuchten.cpp
				richards_plugin.cpp)

# Unit tests.
set(SOURCES_TEST	van_genuchten.cpp
					tests/richards_unit_test.cpp)


# Dependencies: JSON, autodiff

if(NOT USE_JSON)
	MESSAGE(WARNING "This plugin requires JSON. Trying default location.")
	include_directories(AFTER ${UG_ROOT_PATH}/externals/JSONForUG4/json-cxx/single_include)  # JSON
else()
	MESSAGE(MESSAGE "Richards: JSON Found.")
endif ()


if(NOT USE_AUTODIFF)
	MESSAGE(WARNING "This plugin requires Autodiff. Trying default location.")
	include_directories(AFTER ${UG_ROOT_PATH}/externals/AutodiffForUG4/autodiff)  # JSON
else ()
	MESSAGE(MESSAGE "Richards: Autodiff found: ${autodiff_DIR}")
endif ()

# This plugin requires C++ 11.
set (CMAKE_CXX_STANDARD 17) # for autodiff
# add_compile_definitions(_HAS_AUTO_PTR_ETC 1)

# Avoid linking error: undefined reference to boost::system::system_category()
# add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)


if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} ug4 ConvectionDiffusion)
	# target_link_libraries (${pluginName} ug4 ConvectionDiffusion adolc)
endif(buildEmbeddedPlugins)

# Dependencies.
# add_dependencies(${pluginName} adolc)


# Test environment (OPTIONAL).
# add_executable(richards_test ${SOURCES_TEST})
# target_link_libraries (richards_test adolc ug4)