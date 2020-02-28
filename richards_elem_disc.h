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


	RichardsElemDisc (const char *functions, const char* subsets)
	: base_type(functions, subsets)
	{
		// m_exDarcyVel = make_sp(new DataExport<MathVector<dim>, dim>(functions));
		// m_exStorage = make_sp(new DataExport<number, dim>(functions));
	}


	///	returns the export of storage terms
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


protected:
	///	Export for the Darcy velocity
	SmartPtr<TVectorData> m_spFlux;

	///	Export for the brine mass fraction
	SmartPtr<TNumberData> m_spStorage;

};

template <typename TDomain>
void SelectVTKData(SmartPtr<VTKOutput<TDomain::dim> > spVtk, SmartPtr<RichardsElemDisc<TDomain> > spElemDisc, const char *myTag)
{
	typedef UserData<number, TDomain::dim> TUserDataNumber;
	typedef UserData<MathVector<TDomain::dim>, TDomain::dim> TUserDataVector;

	spVtk->select_element(spElemDisc->get_flux_data().template cast_static<TUserDataVector> (), std::string("flux").append(myTag).c_str());
	spVtk->select_element(spElemDisc->get_storage_data().template cast_static<TUserDataNumber> (), std::string("saturation").append(myTag).c_str());
};

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

	RichardsElemDiscFactory() {}

	//! Create classic Richards equation
	SmartPtr<TRichards> create_richards(const char *functions, const char* subsets)
	{
		SmartPtr<TRichards> base = make_sp<TRichards> (new TRichards(functions, subsets));

		// Saturation and conductivity are functions of capillary head.
		SmartPtr<TScaleAddLinkerNumber>  myCapHead = make_sp<TScaleAddLinkerNumber> (new TScaleAddLinkerNumber());
		myCapHead->add(-1.0, base->value());

		m_spSaturation->set_capillary(myCapHead);
		m_spConductivity->set_capillary(myCapHead);

		// Gravity.
		SmartPtr<TConstUserVector> myGravity = make_sp<TConstUserVector> (new TConstUserVector(0.0));
		myGravity->set_entry(dim-1, 1.0);

		// Flux C*\grad (h + z)
		SmartPtr<TScaleAddLinkerVector> myFlux = make_sp<TScaleAddLinkerVector> (new TScaleAddLinkerVector());;
		myFlux->add(m_spConductivity, base->gradient());
		myFlux->add(m_spConductivity, myGravity);

		SmartPtr<TScaleAddLinkerVector> myFlux2 = make_sp<TScaleAddLinkerVector> (new TScaleAddLinkerVector());
		myFlux2->add(-1.0, myFlux);

		base->set_flux_data(myFlux2);
		base->set_storage_data(m_spSaturation);

		return base;
	}

	void set_saturation(SmartPtr<RichardsSaturation<dim> > data) {m_spSaturation = data;}
	void set_conductivity(SmartPtr<RichardsConductivity<dim> > data) {m_spConductivity = data;}

	const char* get_disc_type() {return "fv1";}

protected:
	SmartPtr<RichardsSaturation<dim> > m_spSaturation;
	SmartPtr<RichardsConductivity<dim> > m_spConductivity;
};

}
}
#endif /* RICHARDSEQPLUGIN_RICHARDS_EQUATION_H_ */
