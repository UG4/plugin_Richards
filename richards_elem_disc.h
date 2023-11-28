#/*
 * richards_equation.h
 *
 *  Created on: 28.01.2020
 *      Author: anaegel
 */


#ifndef RICHARDSEQPLUGIN_RICHARDS_EQUATION_H_
#define RICHARDSEQPLUGIN_RICHARDS_EQUATION_H_

// UG4 lib
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#include "lib_disc/spatial_disc/user_data/linker/scale_add_linker.h"
#include "lib_disc/io/vtkoutput.h"

// Plugins
#include "../ConvectionDiffusion/convection_diffusion_base.h"
#include "../ConvectionDiffusion/fv1/convection_diffusion_fv1.h"

namespace ug{
namespace Richards {


template <typename TDomain>
class RichardsElemDisc : public ug::ConvectionDiffusionPlugin::ConvectionDiffusionFV1<TDomain>
{
public:
	static const int dim = TDomain::dim;
	typedef ug::ConvectionDiffusionPlugin::ConvectionDiffusionFV1<TDomain> base_type;

	typedef CplUserData<number, dim> TNumberData;
	typedef CplUserData<MathVector<dim>, dim> TVectorData;

	typedef SmartPtr<CplUserData<number, dim> > NumberExport;
	typedef SmartPtr<CplUserData<MathVector<dim>, dim> > VectorExport;

	/// CTOR
	RichardsElemDisc (const char *functions, const char* subsets)
	: base_type(functions, subsets), m_spTotalConductivity(SPNULL)
	{
		// m_exDarcyVel = make_sp(new DataExport<MathVector<dim>, dim>(functions));
		// m_exStorage = make_sp(new DataExport<number, dim>(functions));
	}


	///	the export of storage terms
	SmartPtr<TNumberData> get_storage_data() {return m_spStorage;}
	void set_storage_data(SmartPtr<TNumberData> storage) {
		m_spStorage = storage;
		this->set_mass(m_spStorage);
		this->set_mass_scale(0.0);
	}


	///	returns the export of the darcy velocity
	SmartPtr<TVectorData>  get_flux_data() {return m_spFlux;}
	void set_flux_data(SmartPtr<TVectorData> myFlux)
	{
		m_spFlux = myFlux;
		this->set_flux(myFlux);
		this->set_diffusion(0.0);
	}

	/// conductivity
	SmartPtr<TNumberData> get_conductivity_data() {return m_spTotalConductivity;}
	void set_conductivity_data(SmartPtr<TNumberData> myCond)
	{
		m_spTotalConductivity = myCond;
	}


protected:
	///	Export for the Darcy velocity
	SmartPtr<TVectorData> m_spFlux;

	///	Export for the brine mass fraction
	SmartPtr<TNumberData> m_spStorage;

	/// Export for absolute conductivity
	SmartPtr<TNumberData> m_spTotalConductivity;

};

template <typename TDomain>
void SelectVTKData(SmartPtr<VTKOutput<TDomain::dim> > spVtk, SmartPtr<RichardsElemDisc<TDomain> > spElemDisc, const char *myTag)
{
	typedef UserData<number, TDomain::dim> TUserDataNumber;
	typedef UserData<MathVector<TDomain::dim>, TDomain::dim> TUserDataVector;

	spVtk->select_element(spElemDisc->get_flux_data().template cast_static<TUserDataVector> (), std::string("flux").append(myTag).c_str());
	spVtk->select_element(spElemDisc->get_storage_data().template cast_static<TUserDataNumber> (), std::string("saturation").append(myTag).c_str());
};



///! Factory for creating elem discs.
template <typename TDomain>
class RichardsElemDiscFactory
{
public:
	static const int dim = TDomain::dim;
	typedef IElemDisc<TDomain>  TElemDisc;
	typedef ug::ConvectionDiffusionPlugin::ConvectionDiffusionFV1<TDomain> TConvDiff;
	typedef RichardsElemDisc<TDomain> TRichards;

	typedef ConstUserVector<dim> TConstUserVector;
	typedef ConstUserNumber<dim> TConstUserNumber;
	typedef ScaleAddLinker<MathVector<dim>, dim, number> TScaleAddLinkerVector;
	typedef ScaleAddLinker<number, dim, number> TScaleAddLinkerNumber;

	typedef CplUserData<MathVector<dim>, dim> TVectorData;
	typedef CplUserData<number, dim> TNumberData;

	RichardsElemDiscFactory() :
		m_spSaturation(SPNULL), m_spRelConductivity(SPNULL), m_spAbsConductivity(SPNULL),
		m_functions("h"), m_dGravity(-1.0)
	{}

	//! Create classic Richards equation.
	/*!
	 *  Element discretization for:
	 *	\frac{\partial S(-h)}{\partial t}  + \nabla \cdot [- C \nabla (h+z)] = 0
	 *
	 * Input parameters: saturation S, conductivity C = Cabs*Crel, gravity $\nabla z$
	 */
	SmartPtr<TRichards> create_richards(const char* subsets)
	{
		//! Create elem disc.
		SmartPtr<TRichards> base = make_sp<TRichards> (new TRichards(m_functions.c_str(), subsets));

		// Saturation and conductivity are functions of capillary head.
		SmartPtr<TScaleAddLinkerNumber>  myCapHead = make_sp<TScaleAddLinkerNumber> (new TScaleAddLinkerNumber());
		myCapHead->add(-1.0, base->value());

		m_spSaturation->set_capillary(myCapHead);
		m_spRelConductivity->set_capillary(myCapHead);

		// Gravity (= \nabla z).
		SmartPtr<TConstUserVector> myGravity = make_sp<TConstUserVector> (new TConstUserVector(0.0));
		myGravity->set_entry(dim-1, -m_dGravity);

		// Conductivity C=Cabs*Crel.
		SmartPtr<TNumberData> spTotalConductivity;
		SmartPtr<TNumberData> spRelConductivity = m_spRelConductivity.template cast_dynamic<TNumberData> ();
		if (m_spAbsConductivity.valid())
		{
			auto spLinker = make_sp<TScaleAddLinkerNumber> (new TScaleAddLinkerNumber());
			spLinker->add(m_spAbsConductivity, spRelConductivity);
			spTotalConductivity = spLinker;

		} else
		{
			spTotalConductivity = m_spRelConductivity.template cast_dynamic<typename TRichards::TNumberData> ();
		}
		base->set_conductivity_data(spTotalConductivity);


		// Flux C*\grad (h + z).
		SmartPtr<TScaleAddLinkerVector> myFlux = make_sp<TScaleAddLinkerVector> (new TScaleAddLinkerVector());
		myFlux->add(spTotalConductivity, base->gradient());
		myFlux->add(spTotalConductivity, myGravity);

		SmartPtr<TScaleAddLinkerVector> myFlux2 = make_sp<TScaleAddLinkerVector> (new TScaleAddLinkerVector());
		myFlux2->add(-1.0, myFlux);

		base->set_flux_data(myFlux2);
		base->set_storage_data(m_spSaturation.template cast_dynamic<typename TRichards::TNumberData> ());
		UG_ASSERT(m_spSaturation.template cast_dynamic<typename TRichards::TNumberData> ().valid(),
				 "Error: Failed to interprete saturation as number data!")


		return base;
	}

	// Classic Richards.
	typedef IRichardsLinker<dim> TSaturation;
	typedef IRichardsLinker<dim> TConductivity;
	//typedef typename RichardsConductivity<dim>::richards_base_type TConductivity;

	void set_saturation(SmartPtr<TSaturation> data) {m_spSaturation = data;}
	void set_conductivity(SmartPtr<TConductivity> data) {m_spRelConductivity = data;}

	// Extensions.
	void set_abs_conductivity(SmartPtr<TNumberData> data) {m_spAbsConductivity = data;}

	void set_function(const char * functions) { m_functions = functions; }
	void set_gravity(number g) { m_dGravity = g; }

	const char* get_disc_type() {return "fv1";}

protected:
	SmartPtr<TSaturation> m_spSaturation;
	SmartPtr<TConductivity> m_spRelConductivity;

	SmartPtr<TNumberData> m_spAbsConductivity;

	std::string m_functions;
	number m_dGravity;
};

}
}
#endif /* RICHARDSEQPLUGIN_RICHARDS_EQUATION_H_ */
