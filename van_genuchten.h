/*
 * richards.h
 *
 *  Created on: 21.12.2019
 *      Author: anaegel
 */


#ifndef RICHARDS__VAN_GENUCHTEN_H_
#define RICHARDS__VAN_GENUCHTEN_H_

// Standard lib.
#include <cmath>

// Use ADOL-C for automatic differentiation.
// #include <adolc/adtl.h>

// Use autodiff for automatic differentiation.
#include <autodiff/forward/dual.hpp>



// UG4 lib.
#include "common/util/smart_pointer.h"
#include "common/assert.h"
#include "registry/class.h"

// My libs.
#include <nlohmann/json.hpp>
#include "json_basics.hh"



namespace ug{
namespace Richards{

typedef autodiff::dual dual;

// const double one = 1.0;



// Parameters for a van Genuchten model.
struct VanGenuchtenParameters
{
	double alpha;
	double n;
	double m; 		// default: 1.0 - (1.0/n);}

	double thetaS;  // default: 1.0
	double thetaR;  // default: 0.0

	double Ksat;    // saturated conductivity K_s (optional)
};


// Parameters for Haverkamp model.
struct HaverkampParameters
{
	double alpha;
	double n;

	double beta;
	double m;

	double thetaS;  // default: 1.0
	double thetaR;  // default: 0.0

	double Ksat;    // saturated conductivity K_s (optional)
};



#ifdef UG_JSON
/// JSON serialize.
void to_json(JSONType& j, const VanGenuchtenParameters& p);

/// JSON de-serialize.
 void from_json(const JSONType& j, VanGenuchtenParameters& p);

 /// JSON serialize.
 void to_json(JSONType& j, const HaverkampParameters& p);

 /// JSON de-serialize.
  void from_json(const JSONType& j, HaverkampParameters& p);
#endif



#ifdef UG_JSON
 /// Some JSON
 void CreateJSONMap(const JSONType &array, std::map<std::string, JSONType> &map);
#endif

/*
class BrooksCorey
{
public:
	BrooksCorrey()
	: m_pd(1.0), m_lambda(1.0)
	{}

	double capillary_pressure(double Sw_eff)
	{
		return pow(m_pd*Sw_eff, 1.0/m_lambda);
	}



protected:
	double m_pd; 			// entry pressure
	double m_lambda;		// ~pore size distribution


};
*/


inline void UGCheckValues(const dual &val)
{
/*	if (!std::isfinite(val.getValue()) || std::isnan(val.getValue()))
	{
		std::cerr << "WARNING: Invalid value:" << val.getValue() << std::endl;
		UG_ASSERT(std::isfinite(val.getValue()), "Function is not bounded");
	}

	if (!std::isfinite(*(val.getADValue())) || std::isnan(*val.getADValue()) )
	{
		std::cerr << "WARNING: Invalid value:" << val.getADValue() << std::endl;
		UG_ASSERT(std::isfinite(*(val.getADValue())), "Derivative is not bounded");
	}
	*/
}



//! This is the interface for a Richards-type model. All derived classes use CRTP for evaluation.
/*! Two non-dimensional quantities are returned
 *  - Saturation $ 0 \le S(\psi) \le 1$
 *  - Relative permeability 0 \le k_r(\psi) \le 1
 *
 *  For convenience, we also return
 *  - Conductitivity K = K_s *k_r(\psi)
 */
template <typename TDerived>
class IRichardsModel
{

protected:
	// CRTP functions.
	const TDerived* me() const
	{ return static_cast<const TDerived*>(this); }

	TDerived* me()
	{ return static_cast<TDerived*>(this); }

	template <typename TFunc>
	void get_value_and_deriv(TFunc F, double H, double &f, double &df) const
	{
		dual h= H;
		f = F(h).val;
		df = derivative(F, wrt(h), at(h));
	};

public:
	//! Calls Saturation_
	void get_saturation(double H, double &S, double &dSdH) const
	{
		auto sat = [this](dual h){ return me()->Saturation_(h);};
		get_value_and_deriv(sat, H, S, dSdH);
	}

	void get_saturations(const double *H, double *S, double *dSdH, size_t n) const
	{
		for (size_t i=0; i<n; ++i)
		{ this->get_saturation(H[i], S[i], dSdH[i]); }
	}

	//! Calls RelativePermeability_
	void get_relative_permeability(double H, double &Kr, double &dKrdH) const
	{
		auto perm = [this](dual h){ return me()->RelativePermeability_(h);};
		get_value_and_deriv(perm, H, Kr, dKrdH);
	}

	void get_relative_permeabilities(const double *H, double *K, double *dKdH, size_t n) const
	{
		for (size_t i=0; i<n; ++i)
		{ this->get_relative_permeability(H[i], K[i], dKdH[i]); }
	}

	//! Calls Conductivity_
	void get_conductivity(double H, double &K, double &dKdH) const
	{
		auto cond = [this](dual h){ return me()->Conductivity_(h);};
		get_value_and_deriv(cond, H, K, dKdH);

	}

	void get_conductivities(const double *H, double *K, double *dKdH, size_t n) const
	{
		for (size_t i=0; i<n; ++i)
		{ this->get_conductivity(H[i], K[i], dKdH[i]); }
	}


public: // The following functions can be used from LUA

	// Saturation
	double Saturation(double H) const
	{
		auto func = [this](dual h){ return me()->Saturation_(h);};
		dual h=H;
		return func(h).val;
	}
	double dSaturation_dH(double H) const
	{
		auto func = [this](dual h){ return me()->Saturation_(h);};
		dual h=H;
		return derivative(func, wrt(h), at(h));
	}

	/// Conductivity K = Ks*kr
	double Conductivity(double H)
	{
		auto func = [this](dual h){ return me()->Conductivity_(h);};
		dual h=H;
		return func(h).val;
	}

	double dConductivity_dH(double H)
	{
		dual h=H;
		auto func = [this](dual h){ return me()->Conductivity_(h);};
		return derivative(func, wrt(h), at(h));
	}

};



struct BrooksCorreyFunctions
{
	static dual ComputeEffectiveSaturation(dual pc, double lambda, double pb)
	{
		dual val = pow(pb/pc, lambda);
		UGCheckValues(val);
		return val;
	};
};

struct VanGenuchtenFunctions
{
	/// Effective (reduced) saturation $0 \le \hat S \le 1$
	/** \hat S:= (1/1+(alpha * psi )^n)^m,  */
	static dual ComputeEffectiveSaturation(dual psi_, double alpha, double n, double m)
	{
			UGCheckValues(psi_);
			dual apn = alpha*psi_;
			UGCheckValues(apn);

			apn = pow(apn, n);
			UGCheckValues(apn);

			dual val = pow(1.0/(1.0+apn), m);
			UGCheckValues(val);
			return val;
	}

	//! Two argument version.
	static dual ComputeEffectiveSaturation(dual psi_, double alpha, double n)
	{ return ComputeEffectiveSaturation(psi_, alpha, n, 1.0-(1.0/n)); }
};

//! Base class for a parameterized model. Provides serialization.
template <typename TParameter>
struct IParameterizedModel
{
	IParameterizedModel() {}

	IParameterizedModel(const TParameter &p) : m_param(p) {}

	void set_parameters(const TParameter &p)
	{ m_param = p; }

	const TParameter& get_parameters() const
	{ return m_param; }

	std::string config_string() const
	{
		nlohmann::json jaux = m_param;
		return jaux.dump();
	}

	TParameter m_param;
};



/// Implements a van Genuchten-Mualem model.
class VanGenuchtenModel
: public IRichardsModel<VanGenuchtenModel>,
  public IParameterizedModel<VanGenuchtenParameters>
{
public:
	typedef IRichardsModel<VanGenuchtenModel> base_type;
	typedef IParameterizedModel<VanGenuchtenParameters> parameterized_model_type;
	typedef VanGenuchtenParameters parameter_type;

	friend base_type;
	friend parameterized_model_type;

#ifdef UG_JSON
	// TODO: Move to Factory
	VanGenuchtenModel(const char* json)
	{
		nlohmann::json j = nlohmann::json::parse(json);
		VanGenuchtenParameters p = j.get<VanGenuchtenParameters>();
		set_parameters(p);
	}
#endif

	VanGenuchtenModel(const parameter_type &p) :  parameterized_model_type(p)
	{}


protected:

	// p_c = p_a - p_w         : Kapillardruck
	// psi = p_c / (rho_w * g) : Kapillardruckhoehe

	/// Effective saturation: $0 \le \hat S \le 1$.
	inline dual EffSaturation_(dual psi_) const
	{
		const VanGenuchtenParameters &p = m_param;
		return VanGenuchtenFunctions::ComputeEffectiveSaturation(psi_, p.alpha, p.n, p.m);
	}

	/// Rescaled Saturation: $$ S:=\theta_r+ (\theta_s - \theta_r) * \hat S$$.
	inline dual Saturation_(dual psi_) const
	{
		const VanGenuchtenParameters &p = m_param;
		UGCheckValues(psi_);

		if (psi_ <= 0) return p.thetaS;
		else
		{
			dual Seff = EffSaturation_(psi_);
			UGCheckValues(Seff);
			return p.thetaR + (p.thetaS-p.thetaR)*Seff;
		}
	}

	/// Relative permeability:  $0 \le k_r \le 1$
	inline dual RelativePermeability_(dual psi_) const
	{
		const VanGenuchtenParameters &p = m_param;
		UGCheckValues(psi_);

		if (psi_<= 0) return 1.0;

		dual Seff = EffSaturation_(psi_); 	// Seff -> 1 is dangerous!
		if (fabs(Seff.val-1.0)<1e-15) return 1.0;
		UGCheckValues(Seff);

		double m = p.m;
		dual brackS = 1.0-pow(Seff,1.0/m);  // => brackS -> 0.0
		UGCheckValues(brackS);

		dual auxS = pow(brackS,m); 			// => this derivative may explode!
		UGCheckValues(auxS);

		return sqrt(Seff)*(1.0-auxS)*(1.0-auxS); // check!
	}

	// Conductivity C:=K_{sat}*k_r
	inline dual Conductivity_(dual psi_) const
	{
		const VanGenuchtenParameters &p = m_param;
		return p.Ksat*RelativePermeability_(psi_);
	}

};



/// Implements a Haverkamp model.
class HaverkampModel  : public IRichardsModel<HaverkampModel>, public IParameterizedModel<HaverkampParameters>
{
public:
	typedef IRichardsModel<HaverkampModel> base_type;
	typedef IParameterizedModel<HaverkampParameters> parameterized_model_type;
	typedef HaverkampParameters parameter_type;

	friend base_type;
	friend parameterized_model_type;

	HaverkampModel(const parameter_type &p) : parameterized_model_type(p)
	{}

#ifdef UG_JSON
	// TODO: Move to Factory
	HaverkampModel(const char* json)
	{
		nlohmann::json j = nlohmann::json::parse(json);
		parameter_type p = j.get<parameter_type>();
		set_parameters(p);
	}
#endif



protected:


	// p_c = p_a - p_w         : Kapillardruck
	// psi = p_c / (rho_w * g) : Kapillardruckhoehe


	/// Rescaled Saturation $$ S:=\theta_r+ (\theta_s - \theta_r) * \hat S $$
	dual Saturation_(dual psi_) const
	{
		const parameter_type &p = get_parameters();
		if (psi_ <= 0.0) return p.thetaS;

		dual Seff = VanGenuchtenFunctions::ComputeEffectiveSaturation(psi_, p.alpha, p.n, 1.0);
		return p.thetaR + (p.thetaS-p.thetaR)*Seff;

	}

	/// Relative permeability  $0 \le k_r \le 1$
	dual RelativePermeability_(dual psi_) const
	{
		const parameter_type &p = get_parameters();
		if (psi_<= 0.0) return 1.0;

		return VanGenuchtenFunctions::ComputeEffectiveSaturation(psi_, p.beta, p.m, 1.0);
	}

	// Conductivity K:=K_{sat}*k_r
	dual Conductivity_(dual psi_) const
	{
		const parameter_type &p = get_parameters();
		return p.Ksat*RelativePermeability_(psi_);
	}

};


/// Implements a van Genuchten model.
class GardnerModel :
		public IParameterizedModel<HaverkampParameters>,
		public IRichardsModel<GardnerModel>
{
public:
	typedef IRichardsModel<GardnerModel> base_type;
	typedef IParameterizedModel<HaverkampParameters> parameterized_model_type;
	typedef VanGenuchtenParameters parameter_type;

	GardnerModel(const VanGenuchtenParameters &p) : m_param(p) {}

	// Saturation
	void get_saturation(double H, double &S, double &dSdH) const
	{
		auto sat = [this](dual h){ return me()->Saturation_(h, m_param);};
		base_type::get_value_and_deriv(sat, H, S, dSdH);

		/*dual head_ = H;
		auto mySat = [this](dual h){ return Saturation_(h, m_param);};

		dual sat = mySat(head_);
		S = sat.val;
		dSdH = derivative(mySat, wrt(head_), at(head_));*/

		//UGCheckValues();
	}

	/*
	void get_saturations(const double *H, double *S, double *dSdH, size_t n) const
	{
		for (size_t i=0; i<n; ++i)
		{ get_saturation(H[i], S[i], dSdH[i]); }
	}
*/


	/// Conductivity K
	void get_conductivity(double H, double &K, double &dKdH) const
	{
		dual h0 = H;
		auto cond = [this](dual h){ return Conductivity_(h, m_param);};
		base_type::get_value_and_deriv(cond, H, K, dKdH);
	}

/*	void get_conductivities(const double *H, double *K, double *dKdH, size_t n) const
	{
		for (size_t i=0; i<n; ++i)
		{ get_conductivity(H[i], K[i], dKdH[i]); }
	}*/


	std::string config_string() const;

protected:
	static double one;
	VanGenuchtenParameters m_param;


	// p_c = p_a - p_w         : Kapillardruck
	// hpc = p_c / (rho_w * g) : Kapillardruckhoehe

	/// Normalisierter Wassergehalt (reduzierte Sättigung):
	/** 0 \le (1/1+(alpha * psi )^n)^m \le 1*/
	static dual EffSaturation_(dual psi_, const VanGenuchtenParameters &p)
	{ return exp(-p.alpha*psi_); }

	/// Rescaled Saturation
	static dual Saturation_(dual psi_, const VanGenuchtenParameters &p)
	{
		if (psi_ <= 0) return p.thetaS;
		else
		{
			dual Seff = EffSaturation_(psi_, p);
			return p.thetaR + (p.thetaS-p.thetaR)*Seff;
		}
	}

	// Conductivity K_s*k_r
	static dual Conductivity_(dual psi_ , const VanGenuchtenParameters &p)
	{ return p.Ksat* EffSaturation_(psi_, p); }

};




#ifdef UG_JSON


//! Generic model creation.
template <typename TModel, typename TParameter = typename TModel::parameter_type>
SmartPtr<TModel> CreateModel(const char *jstring) {

	SmartPtr<TModel> inst = SPNULL;
	try
	{
			nlohmann::json j = nlohmann::json::parse(jstring);
			TParameter p = j.get<TParameter>();
			inst = make_sp(new TModel(p));
	}
	catch (...)
	{
			std::cout << "Construction failed!" << std::endl;
	}
	return inst;
};



//! Factory function.
SmartPtr<VanGenuchtenModel> CreateVanGenuchtenModel(const char *json);


//! Factory class.
struct RichardsModelFactory {
	SmartPtr<VanGenuchtenModel> create_van_genuchten(const char *json);
	SmartPtr<HaverkampModel> create_haverkamp(const char *json);
};

#endif




}
}



#ifdef UG_JSON
// Finally, we need to define how to construct from JSON.
namespace nlohmann {

  template <>
  struct adl_serializer<ug::Richards::VanGenuchtenModel>
  {
       static ug::Richards::VanGenuchtenModel from_json(const json& j)
       { return j.get<typename ug::Richards::VanGenuchtenParameters>(); }; // initialize from model }

       static void to_json(json& j, typename ug::Richards::VanGenuchtenModel m)
       { j = m.get_parameters();}
   };
};
#endif

#endif /* RICHARDS__VAN_GENUCHTEN_H_ */
