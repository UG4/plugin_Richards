/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */


#if __cplusplus >= 201703L
#define _HAS_AUTO_PTR_ETC 0 // activate for boost, if c++17 and beyond
#endif

#include <string>

#include "bridge/util.h"

// replace this with util_domain_dependent.h or util_algebra_dependent.h
// to speed up compilation time
#include "bridge/util_domain_algebra_dependent.h"

#include "common/ug_config.h"
#include "common/error.h"


#include "van_genuchten.h"
#include "richards_linker.h"
#include "richards_elem_disc.h"

#ifdef UG_PARALLEL
#include "pcl/pcl_util.h"
#endif



using namespace std;
using namespace ug::bridge;

namespace ug{
namespace Richards{

/** 
 *  \defgroup sample_plugin Sample
 *  \ingroup plugins
 *  This is a small sample plugin.
 *  \{
 */

void RichardsSaysHello()
{
#ifdef UG_PARALLEL
	pcl::SynchronizeProcesses();
	cout << "Hello, I'm your plugin on proc " <<
				pcl::ProcRank() << "." << endl;
	pcl::SynchronizeProcesses();
#else
	UG_LOG("Hello, I'm your personal plugin in serial environment!\n");
#endif
}

void CrashFct(const string& reason)
{
	UG_THROW("I crash because: "<< reason);
}

void CrashFctFatal(const string& reason)
{
	UG_THROW("I crash fatal because: "<< reason);
}

void RichardsCrashes()
{
	try{
		CrashFct("Some funny reason");
	}
	catch(bad_alloc& err)
	{
		UG_LOG("Bad Alloc");
	}
	UG_CATCH_THROW("Something wrong in PluginCrashes");
}

void RichardsCrashesFatal()
{
	try{
		CrashFctFatal("Some fatal reason");
	}
	catch(bad_alloc& err)
	{
		UG_LOG("Bad Alloc");
	}
	UG_CATCH_THROW("Something wrong in PluginCrashesFatal");
}

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts
 * of the plugin. All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
//	useful defines
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

}

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
//	useful defines
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

	{
		string name = string("RichardsElemDisc").append(suffix);
		typedef RichardsElemDisc<TDomain> T;
		typedef typename RichardsElemDisc<TDomain>::base_type TBase;

		reg.add_class_<T,TBase>(name, grp)
					.template add_constructor<void (*)(char*, char*) >("")
					.add_method("set_storage_data", &T::set_storage_data)
					.add_method("get_storage_data", &T::get_storage_data)
					.add_method("set_flux_data", &T::set_flux_data)
					.add_method("get_flux_data", &T::get_flux_data)
					.add_method("set_conductivity_data", &T::set_conductivity_data)
					.add_method("get_conductivity_data", &T::get_conductivity_data)
					// . add_method("set_conductivity", &T::set_conductivity)
					.set_construct_as_smart_pointer(true);

				reg.add_class_to_group(name, "RichardsElemDisc", tag);

	}

	{
			string name = string("RichardsElemDiscFactory").append(suffix);
			typedef RichardsElemDiscFactory<TDomain> T;

			reg.add_class_<T>(name, grp)
						.template add_constructor<void (*)() >("")
						.add_method("create_richards", &T::create_richards)
						.add_method("set_saturation", &T::set_saturation)
						.add_method("set_conductivity", &T::set_conductivity)
						.add_method("set_abs_conductivity", &T::set_abs_conductivity)
						.add_method("set_function", &T::set_function)
						.add_method("set_gravity", &T::set_function)
						.set_construct_as_smart_pointer(true);

					reg.add_class_to_group(name, "RichardsElemDiscFactory", tag);

		}
	{
		string name = string("SelectVTKData").append(suffix);
		reg.add_function(name, &SelectVTKData<TDomain>, grp);
	}
}

/**
 * Function called for the registration of Dimension dependent parts
 * of the plugin. All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename T, int dim>
static void add_user_data(const char* idstring, Registry& reg, const string grp, const string suffix, const string tag)
{
	string name = string(idstring).append(suffix);
	//typedef RichardsSaturation<dim> T;
	typedef typename T::user_data_base_type TUserData;
	typedef typename T::base_type::TModel TModel;
	typedef IRichardsLinker<dim> TCap;

	reg.add_class_<T,TUserData,TCap>(name, grp)
		.template add_constructor<void (*)(const TModel&) >("")
		.add_method("set_capillary", &T::set_capillary)
		.set_construct_as_smart_pointer(true);

	reg.add_class_to_group(name, idstring, tag);
}

template <int dim>
static void Dimension(Registry& reg, string grp)
{
//	useful defines
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

	{
		// This is an abstract base class for pressure dependent quantities, e.g. saturations.
		string name = string("IRichardsLinker").append(suffix);
		typedef IRichardsLinker<dim> T;
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IRichardsLinker", tag);

	}
	// user data
	{
		// van Genuchten
		typedef RichardsSaturation<dim> S;
		add_user_data<S, dim>("RichardsSaturation", reg, grp, suffix, tag);

		typedef RichardsConductivity<dim> C;
		add_user_data<C, dim>("RichardsConductivity", reg, grp, suffix, tag);
	}

	{
		// Haverkamp
		typedef HaverkampSaturation<dim> S;
		add_user_data<S, dim>("HaverkampSaturation", reg, grp, suffix, tag);

		typedef HaverkampConductivity<dim> C;
		add_user_data<C, dim>("HaverkampConductivity", reg, grp, suffix, tag);
	}

	{
		// Exponential
		typedef ExponentialSaturation<dim> S;
		add_user_data<S, dim>("ExponentialSaturation", reg, grp, suffix, tag);

		typedef ExponentialConductivity<dim> C;
		add_user_data<C, dim>("ExponentialConductivity", reg, grp, suffix, tag);
	}




	{
		string name = string("UserDataFactory").append(suffix);
		typedef UserDataFactory<dim> T;
		reg.add_class_<T>(name, grp)
		.template add_constructor<void (*)() >("")
		// Exponential
		.add_method("create_saturation", static_cast<typename T::return_type (T::*)(const ExponentialModel &)>(&T::create_saturation) )
		.add_method("create_conductivity", static_cast<typename T::return_type (T::*)(const ExponentialModel &)>(&T::create_conductivity) )
		// van Genuchten-Mualem
		.add_method("create_saturation", static_cast<typename T::return_type (T::*)(const VanGenuchtenModel &)>(&T::create_saturation) )
		.add_method("create_conductivity", static_cast<typename T::return_type (T::*)(const VanGenuchtenModel &)>(&T::create_conductivity) )

		.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "UserDataFactory", tag);
	}


	{
		string name = string("OnSurfaceCondition").append(suffix);
		typedef OnSurfaceCondition<dim> T;
		typedef typename T::user_data_base_type TUserData;

		reg.add_class_<T,TUserData>(name, grp)
		   .template add_constructor<void (*)(SmartPtr<typename T::TVectorData>) >("")
		   .set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "OnSurfaceCondition", tag);

	}

}

/**
 * Function called for the registration of Algebra dependent parts
 * of the plugin. All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
//	useful defines
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

}

/**
 * Function called for the registration of Domain and Algebra independent parts
 * of the plugin. All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{
	reg.add_function("RichardsSaysHello", &RichardsSaysHello, grp)
		.add_function("RichardsCrashes", &RichardsCrashes, grp)
		.add_function("RichardsCrashesFatal", &RichardsCrashesFatal, grp)
	    .add_function("CreateVanGenuchtenModel", &CreateVanGenuchtenModel, grp);


	{
		 typedef ExponentialModel T;
		 typedef typename T::parameter_type P;
		 reg.add_class_<P>(std::string("ExponentialModelParameters"), grp)
			.add_constructor<void (*)() >("json-string containing the parameters", "", "", "");


		 string name = string("ExponentialModel");
		 reg.add_class_<T>(name, grp)
			.add_constructor<void (*)(const P&) >("json-string containing the parameters", "", "", "")
			.add_method("config_string", &T::config_string)
			.add_method("saturation", &T::Saturation)
			.add_method("saturation_deriv", &T::dSaturation_dH)
			.add_method("permeability", &T::Saturation)
			.add_method("permeability_deriv", &T::dSaturation_dH)
			.add_method("conductivity", &T::Conductivity)
		   	.add_method("conductivity_deriv", &T::dConductivity_dH);
	}

	{
	   	 string name = string("VanGenuchtenModel");
	   	 typedef VanGenuchtenModel T;
	   	 reg.add_class_<T>(name, grp)
			.add_constructor<void (*)(const char*) >("json-string containing the parameters", "", "", "")
			.add_method("config_string", &T::config_string)
		    .add_method("saturation", &T::Saturation)
			.add_method("saturation_deriv", &T::dSaturation_dH)
			.add_method("conductivity", &T::Conductivity)
	   	 	.add_method("conductivity_deriv", &T::dConductivity_dH);
	}

	{
	   	 string name = string("HaverkampModel");
	   	 typedef HaverkampModel T;
	   	 reg.add_class_<T>(name, grp)
			.add_constructor<void (*)(const char*) >("json-string containing the parameters", "", "", "")
			.add_method("config_string", &T::config_string)
		    .add_method("saturation", &T::Saturation)
			.add_method("saturation_deriv", &T::dSaturation_dH)
			.add_method("conductivity", &T::Conductivity)
	   	 	.add_method("conductivity_deriv", &T::dConductivity_dH);
	}

	{
		 string name = string("RichardsModelFactory");
		 typedef RichardsModelFactory T;
		 reg.add_class_<T>(name, grp)
			.template add_constructor<void (*)() >("")
			.add_method("create_exponential", &T::create_van_genuchten)
			.add_method("create_van_genuchten", &T::create_van_genuchten)
			.add_method("create_haverkamp", &T::create_haverkamp)
			.set_construct_as_smart_pointer(true);
	}



}

}; // end Functionality

// end group sample_plugin
/// \}

}// end of namespace Richards


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_Richards(Registry* reg, string grp)
{
	grp.append("/Richards");
	typedef Richards::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(*reg,grp);
		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

extern "C" UG_API void
FinalizeUGPlugin_Richards()
{
}

}//	end of namespace ug
