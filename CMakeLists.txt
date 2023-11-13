# Copyright (c) 2011-2020:  G-CSC, Goethe University Frankfurt
# Author: Arne Naegel
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

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


# include_directories(AFTER ${UG_ROOT_PATH}/externals/plugin_ADTools/include)
# include_directories(AFTER ${UG_ROOT_PATH}/plugins/autodiff)
# include_directories(AFTER ${UG_ROOT_PATH}/plugins/ADTools/autodiff)
# link_directories(${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})   # ADOL-C


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