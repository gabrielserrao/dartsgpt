//--------------------------------------------------------------------------
#ifndef OPENDARTS_FLASH_EOS_HELMHOLTZ_MIX_H
#define OPENDARTS_FLASH_EOS_HELMHOLTZ_MIX_H
//--------------------------------------------------------------------------

#include "dartsflash/global/global.hpp"
#include "dartsflash/global/components.hpp"

class Mix
{
protected:
	int nc;

	double omegaA, omegaB;
	std::vector<double> Pc, Tc, kappa, kij;

	double N;
	std::vector<double> a_c, a_i, b_i, a_ij, b_ij;

    double B_, D_;
    std::vector<double> B_i, D_i, B_ij, D_ij, D_iT;
    double D_T, D_TT;

public:
	Mix(CompData& data, double omegaa, double omegab, std::vector<double>& kappa_);
    virtual ~Mix() = default;

    Mix* getCopy() { return new Mix(*this); }

    void zeroth_order(std::vector<double>::iterator n_it);
    void first_order(std::vector<double>::iterator n_it);
    void second_order(double T, std::vector<double>::iterator n_it);

    double B(std::vector<double>::iterator n_it);
    double D(std::vector<double>::iterator n_it);

    double B() { return B_; }
    double D() { return D_; }
    double Bi(int i) { return B_i[i]; }
    double Di(int i) { return D_i[i]; }
    double Bij(int i, int j) { return B_ij[i*nc + j]; }
    double Dij(int i, int j) { return D_ij[i*nc + j]; }
    double DiT(int i) { return D_iT[i]; }
    double DT() { return D_T; }
    double DTT() { return D_TT; }

    double ai(double T, int i);
    double aij(int i, int j);
	std::vector<double> ai(double T);
    std::vector<double> aij(double T);

    double ac(int i);
    std::vector<double> ac();
    double alpha(double T, int i);
    double dalpha_dT(double T, int i);
    double d2alpha_dT2(double T, int i);
    
    double bi(int i);
    double bij(int i, int j);
    std::vector<double> bi();
    std::vector<double> bij();

    double dai_dT(double T, int i);
    double d2ai_dT2(double T, int i);
    double daij_dT(double T, int i, int j);
    double d2aij_dT2(double T, int i, int j);
};

//--------------------------------------------------------------------------
#endif // OPENDARTS_FLASH_EOS_HELMHOLTZ_MIX_H
//--------------------------------------------------------------------------
