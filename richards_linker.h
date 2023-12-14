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
#include "lib_disc/spatial_disc/user_data/const_user_data.h"

// Plugin
#include "van_genuchten.h"

#include "json_basics.hh"

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


//! Model traits define a contract mapping TModel to saturations/conductivities
template <typename TModel, int dim>
struct model_traits
{
	typedef void saturation_type;
	typedef void conductivity_type;
};

////////////////////////////////////////////////////////////////////////////////
// Richards linker
////////////////////////////////////////////////////////////////////////////////

//! This is a 'dummy' base class. It indicates a pressure dependent linker.
template <int dim>
struct IRichardsLinker
{
	IRichardsLinker() : m_spCapillary(NULL), m_spDCapillary(NULL) {};
	virtual ~IRichardsLinker(){};

	typedef CplUserData<number, dim> TCplUserData;
	SmartPtr<TCplUserData> m_spCapillary;
	SmartPtr<DependentUserData<number, dim> > m_spDCapillary;

	virtual void set_capillary(SmartPtr<CplUserData<number, dim> > data) = 0;
};



//! Prototype of a linker. Returns values depending on the Functor class.
template <int dim, class TFunctor>
class RichardsLinker
	: public StdDataLinker< RichardsLinker<dim, TFunctor>, number, dim>,
	  public IRichardsLinker<dim>
{
public:
	///	Base class type
		typedef StdDataLinker< RichardsLinker<dim, TFunctor>, number, dim> base_type;
		typedef number data_type;
		typedef typename TFunctor::TModel TModel;
		typedef IRichardsLinker<dim> richards_base_type;

		enum inputs { _H_=0, _SIZE_};


	public:
		typedef typename base_type::base_type user_data_base_type;

		RichardsLinker(const TModel &model) :
			richards_base_type(), m_model(model)
		{
			// UG_LOG("RichardsLinker::RichardsLinker" << std::endl);
			//this linker needs exactly one input
			this->set_num_input(_SIZE_);
		}


		inline void evaluate (data_type& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si) const
		{
		// 	UG_LOG("RichardsLinker::evaluate1: " << std::endl);
			number cap;
			double dummy;
			(*richards_base_type::m_spCapillary)(cap, globIP, time, si);
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
			(*richards_base_type::m_spCapillary)(&vCapVal[0], vGlobIP, time, si,
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
			size_t sid=this->series_id(_H_,s);

			// Checks.
			UG_ASSERT(s >=0, "Huhh: Requires non-negative s");
			UG_ASSERT(static_cast<size_t>(sid) < richards_base_type::m_spCapillary->num_series(), "Huhh: Requires data m_spCapillary!");
			UG_ASSERT(static_cast<size_t>(sid)  < richards_base_type::m_spDCapillary->num_series(), "Huhh: Requires data m_spCapillary!");
			UG_ASSERT(nip == this->m_spCapillary->num_ip(sid),
					"Huhh: Requires data m_spCapillary:" << nip << "!=" << this->m_spCapillary->num_ip(sid));
			UG_ASSERT(nip == this->m_spDCapillary->num_ip(sid),
					"Huhh: Requires data m_spCapillary:" << nip << "!=" << this->m_spDCapillary->num_ip(sid));



			//	Get the data from ip series.
			const number* vH = richards_base_type::m_spCapillary->values(sid);
			number vdSdH[nip];
			TFunctor::get_func_values(m_model, vH, vValue, vdSdH, nip);


			// Compute the derivatives at all ips.

			// Check, if something to do.
			if((!bDeriv) || this->zero_derivative()) return;

			// Clear all derivative values.
			this->set_zero(vvvDeriv, nip);

			//	Derivatives w.r.t height
			if( richards_base_type::m_spDCapillary.valid() && !richards_base_type::m_spDCapillary->zero_derivative())
			{

			for(size_t ip = 0; ip < nip; ++ip)
				for(size_t fct = 0; fct < richards_base_type::m_spDCapillary->num_fct(); ++fct)
				{
				//	get derivative of  w.r.t. to all functions
					const number* vDHeight = richards_base_type::m_spDCapillary->deriv(sid, ip, fct);

				//	get common fct id for this function
					const size_t commonFct = this->input_common_fct(_H_, fct);

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
		//!	Implements IRichardsLinker interface.
		void set_capillary(SmartPtr<CplUserData<number, dim> > data)
		{
			richards_base_type::m_spCapillary = data;  // for evaluation
			richards_base_type::m_spDCapillary = data.template cast_dynamic<DependentUserData<number, dim> >(); // for derivatives
			base_type::set_input(_H_, data, data);
		}


	public:
		TModel& model() { return m_model; }
		const TModel& model() const { return m_model; }

	protected:
		TModel m_model; //! The underlying model.
};


/*! Maps M::get_saturations -> get_func_values */
template <typename M>
struct SaturationAdapter
{
	typedef M TModel;
	static void get_func_values(const TModel &model, const double *h, double *s, double *dsdh, size_t n)
	{ model.get_saturations(h, s, dsdh, n); }
};

//! Returns conductivities (note: corresponds to relative permeability, iff Ksat=1.0).
/*! Maps M::get_conductivities -> get_func_values */
template <typename M>
struct ConductivityAdapter
{
	typedef M TModel;
	static void get_func_values(const M& model, const double *h, double *k, double *dkdh, size_t n)
	{ model.get_conductivities(h, k, dkdh, n);}
};


/**********************************************
 * Simple exponential model.
 **********************************************/

template <int dim>
struct ExponentialSaturation : public RichardsLinker<dim, SaturationAdapter<ExponentialModel> >
{
	typedef RichardsLinker<dim, SaturationAdapter<ExponentialModel> > base_type;
	ExponentialSaturation(const ExponentialModel &m) : base_type (m) {}
};


template <int dim>
struct ExponentialConductivity : public RichardsLinker<dim, ConductivityAdapter<ExponentialModel> >
{
	typedef RichardsLinker<dim,ConductivityAdapter<ExponentialModel> > base_type;
	ExponentialConductivity(const ExponentialModel &m) : base_type (m) {}
};



template<int dim>
struct model_traits<ExponentialModel, dim>
{
	typedef ExponentialSaturation<dim> saturation_type;
	typedef ExponentialConductivity<dim> conductivity_type;
};

/**********************************************
 * Van Genuchten - Mualem
 **********************************************/

//! Shortcuts for van Genuchten
typedef SaturationAdapter<VanGenuchtenModel> vanGenuchtenSaturationAdapter;
typedef ConductivityAdapter<VanGenuchtenModel> vanGenuchtenConductivityAdapter;

//! van Genuchten classes. (ideally, those could be type-def'ed as well..)
template <int dim>
struct RichardsSaturation : public RichardsLinker<dim,vanGenuchtenSaturationAdapter>
{
	typedef RichardsLinker<dim,vanGenuchtenSaturationAdapter> base_type;
	RichardsSaturation(const typename base_type::TModel &model) : base_type (model) {}
};

template <int dim>
struct RichardsConductivity : public RichardsLinker<dim, vanGenuchtenConductivityAdapter>
{
	typedef RichardsLinker<dim,vanGenuchtenConductivityAdapter > base_type;
	RichardsConductivity(const typename base_type::TModel &model) : base_type(model) {}
};


// TODO: Can we define the  classes as follows?
template <int dim>
using VanGenuchtenSaturation = RichardsLinker<dim,vanGenuchtenSaturationAdapter>;

template <int dim>
using VanGenuchtenConductivity = RichardsLinker<dim,vanGenuchtenConductivityAdapter>;


template<int dim>
struct model_traits<VanGenuchtenModel, dim>
{
	typedef VanGenuchtenSaturation<dim> saturation_type;
	typedef VanGenuchtenConductivity<dim> conductivity_type;
};



/**********************************************
 * Extended van-Genuchten-Mualem
 **********************************************/

typedef SaturationAdapter<ExtendedVanGenuchtenModel> IppischSaturationAdapter;
typedef ConductivityAdapter<ExtendedVanGenuchtenModel> IppischConductivityAdapter;


template <int dim>
using ExtendedVanGenuchtenSaturation = RichardsLinker<dim,IppischSaturationAdapter>;

template <int dim>
using ExtendedVanGenuchtenConductivity  = RichardsLinker<dim,IppischConductivityAdapter>;


template<int dim>
struct model_traits<ExtendedVanGenuchtenModel, dim>
{
	typedef ExtendedVanGenuchtenSaturation<dim> saturation_type;
	typedef ExtendedVanGenuchtenConductivity<dim> conductivity_type;
};


/**********************************************
 * Gardner
 **********************************************/
typedef SaturationAdapter<GardnerModel> GardnerSaturationAdapter;
typedef ConductivityAdapter<GardnerModel> GardnerConductivityAdapter;

template <int dim>
class GardnerSaturation : public RichardsLinker<dim,GardnerSaturationAdapter>
{
public:
	typedef RichardsLinker<dim,GardnerSaturationAdapter> base_type;
	GardnerSaturation(const typename base_type::TModel &model) : base_type(model) {}
};

template <int dim>
class GardnerConductivity : public RichardsLinker<dim,GardnerConductivityAdapter>
{
public:
	typedef RichardsLinker<dim,GardnerConductivityAdapter> base_type;
	GardnerConductivity(const typename base_type::TModel &model) : base_type (model) {}
};


template<int dim>
struct model_traits<GardnerModel, dim>
{
	typedef GardnerSaturation<dim> saturation_type;
	typedef GardnerConductivity<dim> conductivity_type;
};

/**********************************************
 * Haverkamp
 **********************************************/
typedef SaturationAdapter<HaverkampModel> HaverkampSaturationAdapter;
typedef ConductivityAdapter<HaverkampModel> HaverkampConductivityAdapter;

template <int dim>
class HaverkampSaturation : public RichardsLinker<dim,HaverkampSaturationAdapter>
{
public:
	typedef RichardsLinker<dim,HaverkampSaturationAdapter> base_type;
	HaverkampSaturation(const typename base_type::TModel &model) : base_type(model) {}
};

template <int dim>
class HaverkampConductivity : public RichardsLinker<dim,HaverkampConductivityAdapter>
{
public:
	typedef RichardsLinker<dim,HaverkampConductivityAdapter> base_type;
	HaverkampConductivity(const typename base_type::TModel &model) : base_type (model) {}
};


template<int dim>
struct model_traits<HaverkampModel, dim>
{
	typedef HaverkampSaturation<dim> saturation_type;
	typedef HaverkampConductivity<dim> conductivity_type;
};


/**********************************************
 * Haverkamp
 **********************************************/


//! Factory class. This constructs appropriate "UserData" from suitable models.

#define RUDF_SMARTPTR_TYPE(TModel,func,dim)  SmartPtr<typename model_traits<TModel,dim>::func ## _type>


#define RUDF_DEFINE_SATURATION(TModel, TReturn, dim) TReturn\
	create_saturation(const TModel &m)\
	{ auto sat = make_sp(new typename model_traits<TModel,dim>::saturation_type(m)); sat->set_capillary(m_capillary); return sat;}

#define RUDF_DEFINE_CONDUCTIVITY(TModel, TReturn, dim) TReturn\
	create_conductivity(const TModel &m)\
	{ auto sat = make_sp(new typename model_traits<TModel,dim>::conductivity_type(m)); sat->set_capillary(m_capillary); return sat;}

#define RUDF_DEFINE_CREATE(TModel, func, dim) SmartPtr<typename model_traits<TModel,dim>::func ## _type>\
	create_ ## func (const TModel &m)\
	{ auto sat = make_sp(new typename model_traits<TModel,dim>::func ## _type(m)); sat->set_capillary(m_capillary); return sat;}


/*
#define RUDF_DEFINE_CONDUCTIVITY(TModel,dim) SmartPtr<typename model_traits<TModel,dim>::conductivity_type>\
	create_conductivity(const TModel &m)\
	{ auto sat = make_sp(new typename model_traits<TModel,dim>::conductivity_type(m)); sat->set_capillary(m_capillary); return sat;} */

template <int dim>
class UserDataFactory {


public:
	typedef CplUserData<number, dim> TUserDataNumber;
	typedef SmartPtr<TUserDataNumber> return_type;

	typedef SmartPtr<CplUserData<number, dim> > input_type;

	UserDataFactory() : m_capillary(SPNULL) {}
	UserDataFactory(input_type capillary) : m_capillary(capillary) {}

	//! Exponential models.
	/* return_type create_saturation(const ExponentialModel &m)
		{ auto sat = make_sp(new typename model_traits<ExponentialModel, dim>::saturation_type(m)); sat->set_capillary(m_capillary); return sat; }

		return_type create_conductivity(const ExponentialModel &m)
		{ auto cond = make_sp(new typename model_traits<ExponentialModel, dim>::conductivity_type(m)); cond->set_capillary(m_capillary); return cond;}

	 */

	RUDF_DEFINE_CREATE(ExponentialModel,saturation,dim);
	RUDF_DEFINE_CREATE(ExponentialModel,conductivity,dim);


	//! Haverkamp models.
	/*return_type create_saturation(const HaverkampModel &m)
	{ auto sat = make_sp(new HaverkampSaturation<dim>(m)); sat->set_capillary(m_capillary); return sat; }

	return_type create_conductivity(const HaverkampModel &m)
	{ auto cond = make_sp(new HaverkampConductivity<dim>(m)); cond->set_capillary(m_capillary); return cond;}*/

	RUDF_DEFINE_CREATE(HaverkampModel,saturation,dim);
	RUDF_DEFINE_CREATE(HaverkampModel,conductivity,dim);

	// van Genuchten models.
	/*return_type create_saturation(const VanGenuchtenModel &m)
	{ auto sat = make_sp(new VanGenuchtenSaturation<dim>(m)); sat->set_capillary(m_capillary); return sat; }

	return_type create_conductivity(const VanGenuchtenModel &m)
	{ auto cond = make_sp(new VanGenuchtenConductivity<dim>(m)); cond->set_capillary(m_capillary); return cond;}*/

	RUDF_DEFINE_CREATE(VanGenuchtenModel,saturation,dim);
	RUDF_DEFINE_CREATE(VanGenuchtenModel,conductivity,dim);

	// Extended van Genuchten models.
	/*return_type create_saturation(const ExtendedVanGenuchtenModel &m)
	{ auto sat = make_sp(new ExtendedVanGenuchtenSaturation<dim>(m)); sat->set_capillary(m_capillary); return sat; }

	return_type create_conductivity(const ExtendedVanGenuchtenModel &m)
	{ auto cond = make_sp(new ExtendedVanGenuchtenConductivity<dim>(m)); cond->set_capillary(m_capillary); return cond;}
*/

	RUDF_DEFINE_CREATE(ExtendedVanGenuchtenModel,saturation,dim);
	RUDF_DEFINE_CREATE(ExtendedVanGenuchtenModel,conductivity,dim);

protected:
	input_type m_capillary;
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


/*


template<TUserData>
SmartPtr<TUserData> CreateUserData(nlohmann::json j, nlohmann::json jparams)
{
	SmartPtr<TUserDataNumber> inst = SPNULL;
	const int dim = TUserData::dim;

	auto jtype = j.at("type");
	auto jvalue = j.at("value");


	if (jtype.is_string())
	{
		std::string mystring;
		jtype.get_to(mystring);

		// JSON: "value" : {"$ref"  : "#/bar"}
		// LUA value = { "$ref" = "#/bar" }

		if (std::string("RichardsSaturation").compare(mystring))
		{
			typedef RichardsSaturation<dim> TRichardsSaturation;
			typedef typename TRichardsSaturation::parameter_type TParameter;
			TParameter p;

			make_sp(new TRichardsSaturation(p));
		}


	} else if (jtype.is_number())
	{
		double myval;
		jtype.get_to(myval);
		inst = make_sp(new ConstUserNumber<dim>(myval));
	}

	return inst;

};

template <int dim>
SmartPtr<CplUserData<number, dim> > CreateUserDataNumber(const char *jstring) {

	typedef CplUserData<number, dim> TUserDataNumber;

	SmartPtr<TUserDataNumber> inst = SPNULL;
	try
	{
			nlohmann::json j = nlohmann::json::parse(jstring);
			VanGenuchtenParameters p = j.get<VanGenuchtenParameters>();
			inst = make_sp(new VanGenuchtenModel(p));
	}
	catch (...)
	{
			std::cout << "Construction failed!" << std::endl;
	}
	return inst;
};
*/




} // namespace Richards
} // end namespace ug

#ifdef UG_JSON
namespace nlohmann {

   // Construct RichardsSaturation from JSON
   template <int dim>
   struct adl_serializer<ug::Richards::RichardsSaturation<dim> >
   {
	   typedef typename ug::Richards::RichardsSaturation<dim> TRichardsLinker;
	   typedef typename ug::Richards::RichardsSaturation<dim>::model_type TModel;

       static TRichardsLinker from_json(const json& j)
       { return j.get<TModel>(); } // initialize from model }

       static void to_json(json& j, TRichardsLinker s)
       { j = s.model(); }
   };

   // Construct RichardsConductivity from JSON
   template <int dim>
     struct adl_serializer<ug::Richards::RichardsConductivity<dim> >
     {
  	   typedef typename ug::Richards::RichardsConductivity<dim> TRichardsLinker;
  	   typedef typename ug::Richards::RichardsConductivity<dim>::model_type TModel;

         static TRichardsLinker from_json(const json& j)
         { return j.get<TModel>(); } // initialize from model }

         static void to_json(json& j, TRichardsLinker s)
         { j = s.model(); }
     };
};
#endif


#endif /* __H__UG__RICHARDS_PLUGIN__LINKER_H__ */
