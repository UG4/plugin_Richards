/*
 * richards.h
 *
 *  Created on: 21.12.2019
 *      Author: anaegel
 */



#ifndef SAMPLE_RICHARDS_RADU_RICHARDS_H_
#define SAMPLE_RICHARDS_RADU_RICHARDS_H_

// Standard lib.
#include <cmath>

// ADOL-C lib for automatic differentiation.
#include <adolc/adtl.h>


// JSON lib.
#include <nlohmann/json.hpp>

// UG4 lib.
#include "common/util/smart_pointer.h"
#include "common/assert.h"



namespace ug{
namespace Richards{

typedef adtl::adouble adouble;
typedef nlohmann::json JSONType;


struct VanGenuchtenParameters
{
	double thetaS;
	double thetaR;
	double alpha;
	double n;
	double Ksat; // optional

	double m() const {return 1.0 - (1.0/n);}
};

/// JSON serialize.
void to_json(JSONType& j, const VanGenuchtenParameters& p);

/// JSON de-serialize.
 void from_json(const JSONType& j, VanGenuchtenParameters& p);

 /// Some JSON
 void CreateJSONMap(const JSONType &array, std::map<std::string, JSONType> &map);


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

inline void UGCheckValues(const adouble &val)
{
	if (!isnormal(val.getValue()))
		UG_ASSERT(isnormal(val.getValue()), "Value is not bounded");

	if (!isfinite(*(val.getADValue())) )
		UG_ASSERT(isfinite(*(val.getADValue())), "Derivative is not bounded");

	UG_ASSERT(!isnan(val.getValue()), "Value is not bounded");
	UG_ASSERT(!isnan(*(val.getADValue())), "Derivative is not bounded");
}

/// Implements a van Genuchten model.
class VanGenuchtenModel
{
public:
	VanGenuchtenModel(const VanGenuchtenParameters &p) : m_param(p)
	{}

	// Saturation
	double Saturation(double H) const
	{
		adouble head_(H);
		return Saturation_(head_, this->m_param).getValue();
	}
	double dSaturation_dH(double H) const
	{
		adouble head_(H, &one);
		return *Saturation_(head_, this->m_param).getADValue();
	}

	void get_saturation(double H, double &S, double &dSdH) const
	{
		adouble head_(H, &one);
		adouble sat = Saturation_(head_, this->m_param);
		S = sat.getValue();
		dSdH = *(sat.getADValue());
		// std::cerr <<  dSdH << ";";
		// dSdH = -0.017;

	}

	void get_saturations(const double *H, double *S, double *dSdH, size_t n) const
	{
		//std::cerr <<  "dSdH=";
		for (size_t i=0; i<n; ++i)
		{ get_saturation(H[i], S[i], dSdH[i]); }
		// std::cerr <<  std::endl;
	}



	void get_relative_permeability(double H, double &Kr, double &dKrdH) const
	{
			adouble head_(H, &one);
			adouble kr = RelativePermeability_(head_, this->m_param);
			Kr = kr.getValue();
			dKrdH = *(kr.getADValue());
	}

	void get_relative_permeabilities(const double *H, double *K, double *dKdH, size_t n) const
	{
		for (size_t i=0; i<n; ++i)
		{ get_relative_permeability(H[i], K[i], dKdH[i]); }
	}


	/// Conductivity K
		double Conductivity(double H)
		{
			adouble head_(H);
			return Conductivity_(head_, this->m_param).getValue();
		}

		double dConductivity_dH(double H)
		{
			adouble head_(H, &one);
			return *Conductivity_(head_, this->m_param).getADValue();
		}

		void get_conductivity(double H, double &K, double &dKdH) const
		{
			adouble head_(H, &one);
			adouble k = Conductivity_(head_, this->m_param);
			K = k.getValue();
			dKdH = *(k.getADValue());
		}

		void get_conductivities(const double *H, double *K, double *dKdH, size_t n) const
		{
			for (size_t i=0; i<n; ++i)
			{ get_conductivity(H[i], K[i], dKdH[i]); }
		}


	std::string config_string() const;

protected:
	static double one;
	VanGenuchtenParameters m_param;


	// p_c = p_a - p_w         : Kapillardruck
	// hpc = p_c / (rho_w * g) : Kapillardruckhoehe

	/// Normalisierter Wassergehalt (reduzierte Sättigung):
	/** 0 \le (1/1+(alpha * psi )^n)^m \le 1*/
	static adouble EffSaturation_(adouble psi_, const VanGenuchtenParameters &p)
	{
		UGCheckValues(psi_);

		adouble apn = p.alpha*psi_;
		apn = pow(apn, p.n);
		adouble val = pow(1.0/(1.0+apn), p.m());

		UGCheckValues(val);

		return val;
	}

	/// Rescaled Saturation
	static adouble Saturation_(adouble psi_, const VanGenuchtenParameters &p)
	{
		if (psi_ <= 0) return p.thetaS;
		else
		{
			adouble Seff = EffSaturation_(psi_, p);
			return p.thetaR + (p.thetaS-p.thetaR)*Seff;
		}
	}

	/// Relative permeability k_r
	static adouble RelativePermeability_(adouble psi_ , const VanGenuchtenParameters &p)
	{
		if (psi_<= 0) return 1.0;

		adouble Seff = EffSaturation_(psi_, p);
		double m = p.m();
		adouble brackS = 1.0-pow(Seff,1.0/m);
		adouble auxS = pow(brackS,m);
		return sqrt(Seff)*(1.0-auxS)*(1.0-auxS); // check!
	}

	// Conductivity K_s*k_r
	static adouble Conductivity_(adouble psi_ , const VanGenuchtenParameters &p)
	{
		return p.Ksat*RelativePermeability_(psi_,p); // check!
	}



};


/// Implements a van Genuchten model.
class GardnerModel
{
public:
	GardnerModel(const VanGenuchtenParameters &p) : m_param(p)
	{}

	// Saturation
	void get_saturation(double H, double &S, double &dSdH) const
	{
		adouble head_(H, &one);
		adouble sat = Saturation_(head_, this->m_param);
		S = sat.getValue();
		dSdH = *(sat.getADValue());

		UGCheckValues(sat);
	}

	void get_saturations(const double *H, double *S, double *dSdH, size_t n) const
	{
		for (size_t i=0; i<n; ++i)
		{ get_saturation(H[i], S[i], dSdH[i]); }
	}



	/// Conductivity K
	void get_conductivity(double H, double &K, double &dKdH) const
	{
		adouble head_(H, &one);
		adouble k = Conductivity_(head_, this->m_param);
		K = k.getValue();
		dKdH = *(k.getADValue());

		UGCheckValues(k);
	}

	void get_conductivities(const double *H, double *K, double *dKdH, size_t n) const
	{
		for (size_t i=0; i<n; ++i)
		{ get_conductivity(H[i], K[i], dKdH[i]); }
	}


	std::string config_string() const;

protected:
	static double one;
	VanGenuchtenParameters m_param;


	// p_c = p_a - p_w         : Kapillardruck
	// hpc = p_c / (rho_w * g) : Kapillardruckhoehe

	/// Normalisierter Wassergehalt (reduzierte Sättigung):
	/** 0 \le (1/1+(alpha * psi )^n)^m \le 1*/
	static adouble EffSaturation_(adouble psi_, const VanGenuchtenParameters &p)
	{ return exp(-p.alpha*psi_); }

	/// Rescaled Saturation
	static adouble Saturation_(adouble psi_, const VanGenuchtenParameters &p)
	{
		if (psi_ <= 0) return p.thetaS;
		else
		{
			adouble Seff = EffSaturation_(psi_, p);
			return p.thetaR + (p.thetaS-p.thetaR)*Seff;
		}
	}

	// Conductivity K_s*k_r
	static adouble Conductivity_(adouble psi_ , const VanGenuchtenParameters &p)
	{ return p.Ksat* EffSaturation_(psi_, p); }

};



//! Factory function.
SmartPtr<VanGenuchtenModel> CreateVanGenuchtenModel(const char *json);

//! Factory method.
struct VanGenuchtenModelFactory {
	//SmartPtr<VanGenuchtenModel> create_default();
	SmartPtr<VanGenuchtenModel> create(const char *);
};





}
}

#endif /* SAMPLE_RICHARDS_RADU_RICHARDS_H_ */
