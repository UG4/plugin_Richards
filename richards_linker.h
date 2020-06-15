/*
 * Copyright (c) 2020-2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel
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

#ifndef __H__UG__RICHARDS_PLUGIN__LINKER_H__
#define __H__UG__RICHARDS_PLUGIN__LINKER_H__

// UG4
#include "lib_disc/spatial_disc/user_data/linker/linker.h"

// Plugin
#include "van_genuchten.h"

namespace ug{
namespace Richards{

/*
template <typename T>
struct crtp
{
    T& underlying() { return static_cast<T&>(*this); }
    T const& underlying() const { return static_cast<T const&>(*this); }
};

class LucasBase : public crtp<LucasBase>
{


	double evalOutside(double x)
	{
		 return underlying().evalBase()
	};


//	double evalBase(double x)
//	{ return x*x}
}

class Lucas1 : public LucasBase{
	double evalBase(double x)
	{ return x*x}
}
*/
////////////////////////////////////////////////////////////////////////////////
// Richards linker
////////////////////////////////////////////////////////////////////////////////

//! Prototype for a vanGenuchten Linker
template <int dim, class TFunctor>
class RichardsLinker
	: public StdDataLinker< RichardsLinker<dim, TFunctor>, number, dim>
{
	///	Base class type
		typedef StdDataLinker< RichardsLinker<dim, TFunctor>, number, dim> base_type;
		typedef number data_type;
		typedef typename TFunctor::TModel TModel;

		enum inputs { _H_=0, _SIZE_};


	public:
		typedef typename base_type::base_type user_data_base_type;

		RichardsLinker(const TModel &model) :
			m_spCapillary(NULL), m_spDCapillary(NULL), m_model(model)
		{
			//	this linker needs exactly one input
			this->set_num_input(_SIZE_);
			UG_LOG("RichardsLinker::RichardsLinker" << std::endl);
		}


		inline void evaluate (data_type& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si) const
		{
		// 	UG_LOG("RichardsLinker::evaluate1: " << std::endl);
			number cap;
			double dummy;
			(*m_spCapillary)(cap, globIP, time, si);
			TFunctor::get_func_values(m_model, &cap, &value, &dummy, 1);
		}

		template <int refDim>
		inline void evaluate(data_type vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     GridObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     LocalVector* u,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
		// 	UG_LOG("RichardsLinker::evaluate2: " << std::endl);
			number vCapVal[nip];
			(*m_spCapillary)(&vCapVal[0], vGlobIP, time, si,
								elem, vCornerCoords, vLocIP, nip, u, vJT);

			number vdSdH[nip]; // dummies
			TFunctor::get_func_values(m_model,vCapVal, vValue, vdSdH, nip);

		}

		template <int refDim>
		void eval_and_deriv(data_type vValue[],
		                    const MathVector<dim> vGlobIP[],
		                    number time, int si,
		                    GridObject* elem,
		                    const MathVector<dim> vCornerCoords[],
		                    const MathVector<refDim> vLocIP[],
		                    const size_t nip,
		                    LocalVector* u,
		                    bool bDeriv,
		                    int s,
		                    std::vector<std::vector<number > > vvvDeriv[],
		                    const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			// UG_LOG("RichardsLinker::eval_and_deriv: " << std::endl);

			// Checks.
			UG_ASSERT(s >=0, "Huhh: Requires non-negative s");
			UG_ASSERT(static_cast<size_t>(s) < m_spCapillary->num_series(), "Huhh: Requires data m_spCapillary!");
			UG_ASSERT(static_cast<size_t>(s)  < m_spDCapillary->num_series(), "Huhh: Requires data m_spCapillary!");
			/*UG_ASSERT(nip == m_spDCapillary->num_ip(s), "Huhh: Requires data m_spCapillary:"
								<< nip << "!=" << m_spCapillary->num_ip(s));

			UG_ASSERT(nip == m_spCapillary->num_ip(s), "Huhh: Requires data m_spCapillary:"
					<< nip << "!=" << m_spCapillary->num_ip(s));*/

			//	Get the data from ip series.
			const number* vH = m_spCapillary->values(s);
			number vdSdH[nip];
			TFunctor::get_func_values(m_model, vH, vValue, vdSdH, nip);


			// Compute the derivatives at all ips.

			// Check, if something to do.
			if(!bDeriv || this->zero_derivative()) return;

			// Clear all derivative values.
			this->set_zero(vvvDeriv, nip);

			//	Derivatives w.r.t height
			if( m_spDCapillary.valid() && !m_spDCapillary->zero_derivative())
			{



			for(size_t ip = 0; ip < nip; ++ip)
				for(size_t fct = 0; fct < m_spDCapillary->num_fct(); ++fct)
				{
				//	get derivative of  w.r.t. to all functions
					const number* vDHeight = m_spDCapillary->deriv(s, ip, fct);

				//	get common fct id for this function
					const size_t commonFct = this->input_common_fct(_H_, fct);

				//	loop all shapes and set the derivative

					if (this->num_sh(commonFct) == nip)
					{
						// FV1 mass lumping
						vvvDeriv[ip][commonFct][ip] += vdSdH[ip] * vDHeight[ip];
					}
					else
					{
						for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
						{
							UG_ASSERT(commonFct < vvvDeriv[ip].size(), commonFct<<", "<<vvvDeriv[ip].size());
							vvvDeriv[ip][commonFct][sh] += vdSdH[ip] * vDHeight[sh];
						}
					}
				}

			}

		}

	public:
	///	set height import
		void set_capillary(SmartPtr<CplUserData<number, dim> > data)
		{
			m_spCapillary = data;  // for evaluation
			m_spDCapillary = data.template cast_dynamic<DependentUserData<number, dim> >(); // for derivatives
			base_type::set_input(_H_, data, data);
		}

		TModel& model() { return m_model; }
		const TModel& model() const { return m_model; }

	protected:
		///	import for height
		SmartPtr<CplUserData<number, dim> > m_spCapillary;
		SmartPtr<DependentUserData<number, dim> > m_spDCapillary;

		TModel m_model;
};


//! Returns saturations of VanGenuchtenModel
struct SaturationAdapter
{
	typedef VanGenuchtenModel TModel;
	static void get_func_values(const TModel &model, const double *h, double *s, double *dsdh, size_t n)
	{ model.get_saturations(h, s, dsdh, n); }

};

template <int dim>
class RichardsSaturation : public RichardsLinker<dim,SaturationAdapter>
{
public:
	typedef RichardsLinker<dim,SaturationAdapter> base_type;
	RichardsSaturation(const typename SaturationAdapter::TModel &model) : RichardsLinker<dim,SaturationAdapter> (model)
	{}
};


//! Returns conductivities of VanGenuchtenModel
struct ConductivityAdapter
{
	typedef VanGenuchtenModel TModel;
	static void get_func_values(const TModel &model, const double *h, double *k, double *dkdh, size_t n)
	{ model.get_conductivities(h, k, dkdh, n);}
};


template <int dim>
class RichardsConductivity : public RichardsLinker<dim,ConductivityAdapter>
{
public:
	typedef RichardsLinker<dim,ConductivityAdapter> base_type;
	RichardsConductivity(const typename ConductivityAdapter::TModel &model) : RichardsLinker<dim,ConductivityAdapter> (model)
	{}
};






//! Returns saturations of Gardner
struct GardnerSaturationAdapter
{
	typedef GardnerModel TModel;
	static void get_func_values(const TModel &model, const double *h, double *s, double *dsdh, size_t n)
	{ model.get_saturations(h, s, dsdh, n); }
};

template <int dim>
class GardnerSaturation : public RichardsLinker<dim,GardnerSaturationAdapter>
{
public:
	typedef RichardsLinker<dim,GardnerSaturationAdapter> base_type;
	GardnerSaturation(const typename GardnerSaturationAdapter::TModel &model) : RichardsLinker<dim,SaturationAdapter> (model)
	{}
};


//! Returns conductivities of Gardner
struct GardnerConductivityAdapter
{
	typedef GardnerModel TModel;
	static void get_func_values(const TModel &model, const double *h, double *k, double *dkdh, size_t n)
	{ model.get_conductivities(h, k, dkdh, n);}
};


template <int dim>
class GardnerConductivity : public RichardsLinker<dim,GardnerConductivityAdapter>
{
public:
	typedef RichardsLinker<dim,GardnerConductivityAdapter> base_type;
	GardnerConductivity(const typename GardnerConductivityAdapter::TModel &model) : RichardsLinker<dim,GardnerConductivityAdapter> (model)
	{}
};


template <int dim>
class OnSurfaceCondition : public UserData<number, dim, bool>
{

public:
	typedef number TData;
	typedef bool TRet;

	typedef UserData<number, dim, bool> user_data_base_type;
	typedef MathVector<dim> math_vector_type;
	typedef UserData<MathVector<dim>, dim> TVectorData;

	OnSurfaceCondition(SmartPtr<TVectorData> spFlux): m_spFlux(spFlux) {};
	virtual ~OnSurfaceCondition(){}


protected:
	SmartPtr<TVectorData> m_spFlux;

public:
	///	returns if provided data is continuous over geometric object boundaries
			bool continuous() const {return true;}

		///	returns if grid function is needed for evaluation
			 bool requires_grid_fct() const {return false;}

		public:
		///	returns value for a global position
			virtual TRet operator() (TData& value,
									 const MathVector<dim>& globIP,
									 number time, int si) const
			{
				math_vector_type fluxVec;
				m_spFlux->operator()(fluxVec, globIP, time, si);
			//	return m_spFlux->operator(value);
				value = 47.11;
				return (fluxVec[dim]>0) ? true : false;
			}

		///	returns values for global positions
			virtual void operator()(TData vValue[],
									const MathVector<dim> vGlobIP[],
									number time, int si, const size_t nip) const
			{}

			virtual void operator()(TData vValue[],
				                        const MathVector<dim> vGlobIP[],
				                        number time, int si,
				                        GridObject* elem,
				                        const MathVector<dim> vCornerCoords[],
				                        const MathVector<1> vLocIP[],
				                        const size_t nip,
				                        LocalVector* u,
				                        const MathMatrix<1, dim>* vJT = NULL) const
			{}

				virtual void operator()(TData vValue[],
				                        const MathVector<dim> vGlobIP[],
				                        number time, int si,
				                        GridObject* elem,
				                        const MathVector<dim> vCornerCoords[],
				                        const MathVector<2> vLocIP[],
				                        const size_t nip,
				                        LocalVector* u,
				                        const MathMatrix<2, dim>* vJT = NULL) const
				                        {}

				virtual void operator()(TData vValue[],
				                        const MathVector<dim> vGlobIP[],
				                        number time, int si,
				                        GridObject* elem,
				                        const MathVector<dim> vCornerCoords[],
				                        const MathVector<3> vLocIP[],
				                        const size_t nip,
				                        LocalVector* u,
				                        const MathMatrix<3, dim>* vJT = NULL) const
				{}

};
// TODO: Should be replaced!
//! Returns relative permeability of VanGenuchtenModel
/* struct RelativePermeabilityAdapter
{
	typedef VanGenuchtenModel TModel;
	static void get_func_values(const TModel &model, const double *h, double *k, double *dkdh, size_t n)
	{ model.get_relative_permeabilities(h, k, dkdh, n); }
};


template <int dim>
class RichardsRelativePermeability : public RichardsLinker<dim,RelativePermeabilityAdapter>
{
public:
	typedef RichardsLinker<dim,RelativePermeabilityAdapter> base_type;
	RichardsRelativePermeability(const typename RelativePermeabilityAdapter::TModel &model) : RichardsLinker<dim, RelativePermeabilityAdapter> (model)
	{}
};
*/


} // namespace Richards
} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DARCY_VELOCITY_LINKER__ */
