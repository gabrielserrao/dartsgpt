//--------------------------------------------------------------------------
#ifndef OPENDARTS_FLASH_EOS_AQ_ZIABAKHSH_H
#define OPENDARTS_FLASH_EOS_AQ_ZIABAKHSH_H
//--------------------------------------------------------------------------

#include "dartsflash/eos/aq/aq.hpp"

namespace ziabakhsh {
	extern std::vector<double> Psw;
	extern std::unordered_map<std::string, std::vector<double>> labda, ksi;
	extern std::unordered_map<std::string, double> eta, tau, beta, Gamma;
	extern double R;
}

class Ziabakhsh2012 : public AQBase
{
private:
	double m_c, m_ac;
	double V_H2O{ 18.1 }, Mw{ 18.0152 };
	double K0_H2O, lnKw;
	double rho0_H2O, drho0_H2OdT, d2rho0_H2OdT2, drho0_H2OdP, f0_H2O, df0_H2OdT, d2f0_H2OdT2, df0_H2OdP;
	std::vector<double> lnk_H, labda, ksi, dmcdxj, dmacdxj;

public:
	Ziabakhsh2012(CompData& comp_data);

	AQBase* getCopy() override { return new Ziabakhsh2012(*this); }
	
	virtual void init_PT(double p_, double T_, AQEoS::CompType component) override;
	virtual void solve_PT(std::vector<double>& x_, bool second_order, AQEoS::CompType comp_type) override;

	virtual double lnphii(int i) override;
	virtual double dlnphii_dP(int i) override;
	virtual double dlnphii_dT(int i) override;
	virtual double d2lnphii_dT2(int i) override;
	virtual double dlnphii_dxj(int i, int j) override;
	virtual std::vector<double> d2lnphii_dTdxj(int i) override;

	virtual double G_pure(double X, double T_, bool pt=true) override;

	// virtual int derivatives_test(double p_, double T_, std::vector<double>& n_, double tol, bool verbose=false) override;
};

//--------------------------------------------------------------------------
#endif // OPENDARTS_FLASH_EOS_AQ_ZIABAKHSH_H
//--------------------------------------------------------------------------
