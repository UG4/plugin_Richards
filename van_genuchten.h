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

inline void UGCheckValues(const dual &number)
{
	if (!std::isfinite(number.val) || std::isnan(number.val))
	{
		std::cerr << "WARNING: Invalid value:" << number.val << std::endl;
		UG_ASSERT(std::isfinite(number.val), "Function is not bounded");
	}

	if (!std::isfinite(number.grad) || std::isnan(number.grad))
	{
		std::cerr << "WARNING: Invalid value:" << number.grad << std::endl;
		UG_ASSERT(std::isfinite(number.grad), "Derivative is not bounded");
	}

}




// Parameters for a van Genuchten model.
struct VanGenuchtenParameters
{
	double alpha=1.0;

	double n=2.0;
	double m=0.5; 		// default: 1.0 - (1.0/n);}

	double thetaR=0.0;  // default: 0.0
	double thetaS=1.0;  // default: 1.0

	double Ksat=1.0;    // saturated conductivity K_s (optional)

};

#ifdef UG_JSON
	NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(VanGenuchtenParameters, alpha, n, m, thetaS, thetaR, Ksat);
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
		/*f = F(h).val; df = derivative(F, wrt(h), at(h));*/
		auto [F0, Df0] = derivatives(F, wrt(h), at(h));
		f = F0; df = Df0;
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
		dual h=H;
		auto func = [this](dual h){ return me()->Saturation_(h);};
		return func(h).val;
	}
	double dSaturation_dH(double H) const
	{
		dual h=H;
		auto func = [this](dual h){ return me()->Saturation_(h);};
		return derivative(func, wrt(h), at(h));
	}

	/// Conductivity K = Ks*kr
	double Conductivity(double H)
	{
		dual h=H;
		auto func = [this](dual h){ return me()->Conductivity_(h);};
		return func(h).val;
	}

	double dConductivity_dH(double H)
	{
		dual h=H;
		auto func = [this](dual h){ return me()->Conductivity_(h);};
		return derivative(func, wrt(h), at(h));
	}

};



struct BrooksCoreyFunctions
{
	static dual ComputeEffectiveSaturation(dual pc, double lambda, double pb)
	{
		dual val = pow(pb/pc, lambda);
		UGCheckValues(val);
		return val;
	};
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


/*********************************************
 * Exponential model.
 *********************************************/

struct ExponentialModelParameters
{
	double pentry = 1.0;
	double alpha = 1.0;
	double beta = 1.0;

	double thetaR = 0.0;  // default: 0.0
	double thetaS = 1.0;  // default: 1.0

	double Ksat = 1.0;    // saturated conductivity K_s (optional)
};

#ifdef UG_JSON
	NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(ExponentialModelParameters, pentry, alpha, beta, thetaR, thetaS, Ksat)
#endif


/// Implement a simple exponential model.
class ExponentialModel :
	public IRichardsModel<ExponentialModel>,
	public IParameterizedModel<ExponentialModelParameters>
{
public:
	typedef IRichardsModel<ExponentialModel> base_type;
	typedef ExponentialModelParameters parameter_type;
	typedef IParameterizedModel<ExponentialModelParameters> parameterized_model_type;

	friend base_type;
	friend parameterized_model_type;

	ExponentialModel(const parameter_type &p) : parameterized_model_type(p) {}

protected:

	// p_c = p_a - p_w         : Kapillardruck >=0
	// psi = p_c / (rho_w * g) : Kapillardruckhoehe >=0

	/// Effective saturation: $0 \le \hat S \le 1$.
	inline dual EffSaturation_(dual pc) const
	{
		const parameter_type &p = this->get_parameters();
		return exp(-p.alpha * (pc/p.pentry));
	}

	/// Rescaled Saturation: $$ S:=\theta_r+ (\theta_s - \theta_r) * \hat S$$.
	inline dual Saturation_(dual pc) const
	{
		const parameter_type &p = this->get_parameters();
		UGCheckValues(pc);

		if (pc.val <= 0) return p.thetaS; // Medium is saturated
		const dual Seff = EffSaturation_(pc); UGCheckValues(Seff);
		//UG_LOG("seff:" << (pc.val/p.pentry) << " -> " << Seff.val << std::endl);
		return p.thetaR + (p.thetaS-p.thetaR)*Seff;

	}

	/// Relative permeability:  $0 \le k_r \le 1$
	inline dual RelativePermeability_(dual pc) const
	{
		const parameter_type &p = this->get_parameters();
		UGCheckValues(pc);

		if (pc.val<= 0) return 1.0;
		const dual kval = exp(-p.beta * (pc/p.pentry));

		//UG_LOG("k:" << (pc.val/p.pentry) << " -> " << kval.val << std::endl);
		return kval;

	}

	// Conductivity C:=K_{sat}*k_r
	inline dual Conductivity_(dual pc) const
	{
		return get_parameters().Ksat*RelativePermeability_(pc);
	}

};

/*********************************************
 * Van Genuchten-Mualem
 *********************************************/

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



/// Implements a van Genuchten-Mualem model.
class VanGenuchtenModel
: public IRichardsModel<VanGenuchtenModel>, public IParameterizedModel<VanGenuchtenParameters>
{
public:
	typedef IRichardsModel<VanGenuchtenModel> base_type;
	typedef IParameterizedModel<VanGenuchtenParameters> parameterized_model_type;
	typedef VanGenuchtenParameters parameter_type;

	friend base_type;
	friend parameterized_model_type;

	VanGenuchtenModel(const parameter_type &p) :  parameterized_model_type(p) {}


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


/*********************************************
 * Haverkamp model.
 *********************************************/



// Parameters for Haverkamp model.
struct HaverkampParameters
{
	double alpha;
	double n;

	double beta;
	double m;

	double thetaR=0.0;  // default: 0.0
	double thetaS=1.0;  // default: 1.0

	double Ksat=1.0;    // saturated conductivity K_s (optional)
};



#ifdef UG_JSON

 /// JSON serialize.
 void to_json(JSONType& j, const HaverkampParameters& p);

 /// JSON de-serialize.
  void from_json(const JSONType& j, HaverkampParameters& p);
#endif




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

	/// Saturation
	void get_saturation(double H, double &S, double &dSdH) const
	{
		auto sat = [this](dual h){ return me()->Saturation_(h, m_param);};
		base_type::get_value_and_deriv(sat, H, S, dSdH);
	}

	/// Conductivity K
	void get_conductivity(double H, double &K, double &dKdH) const
	{
		auto cond = [this](dual h){ return Conductivity_(h, m_param);};
		base_type::get_value_and_deriv(cond, H, K, dKdH);
	}

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
			UG_THROW ("Construction failed: " << jstring);
	}
	return inst;
};



//! Factory function.
SmartPtr<VanGenuchtenModel> CreateVanGenuchtenModel(const char *json);


//! Factory class.
struct RichardsModelFactory {
	SmartPtr<ExponentialModel> create_exponential(const char *json);
	SmartPtr<VanGenuchtenModel> create_van_genuchten(const char *json);
	SmartPtr<HaverkampModel> create_haverkamp(const char *json);
};

#endif




}


#ifdef UG_JSON
	template <> struct is_json_constructible<typename Richards::ExponentialModelParameters> { const static bool value = true;};
#endif
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
