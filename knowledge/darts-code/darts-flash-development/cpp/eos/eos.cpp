#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cassert>

#include "dartsflash/global/global.hpp"
#include "dartsflash/global/units.hpp"
#include "dartsflash/eos/eos.hpp"
#include <Eigen/Dense>

EoS::EoS(CompData& comp_data) 
{
    this->compdata = comp_data;
    this->nc = comp_data.nc;
    this->ni = comp_data.ni;
    this->ns = nc + ni;
    this->n.resize(ns);
	
    this->dlnphidP.resize(ns);
    this->dlnphidT.resize(ns);
    this->d2lnphidT2.resize(ns);
    this->dlnphidn.resize(ns*ns);
    this->d2lnphidTdn.resize(ns*ns);
    this->gpure.resize(ns);

    this->units = comp_data.units;
}

void EoS::set_eos_range(int i, const std::vector<double>& range)
{
    // Define range of applicability for specific EoS
    this->eos_range[i] = range;
    return;
}
bool EoS::eos_in_range(std::vector<double>::iterator n_it)
{
    // Check if EoS is applicable for specific state according to specified ranges
    for (auto& it: this->eos_range) {
        this->N = std::accumulate(n_it, n_it + this->ns, 0.);
        double xi = *(n_it + it.first) / this->N;
		if (xi < it.second[0] || xi > it.second[1])
        {
            return false;
        }
	}
    return true;
}

void EoS::solve_PT(double p_, double T_, std::vector<double>& n_, int start_idx, bool second_order)
{
    // Calculate composition-independent and composition-dependent EoS-parameters
    this->init_PT(p_, T_);
    this->solve_PT(n_.begin() + start_idx, second_order);
    return;
}
void EoS::solve_VT(double V_, double T_, std::vector<double>& n_, int start_idx, bool second_order)
{
    // Calculate composition-independent and composition-dependent EoS-parameters
    this->init_VT(V_, T_);
    this->solve_VT(n_.begin() + start_idx, second_order);
    return;
}

EoS::RootFlag EoS::is_root_type(std::vector<double>::iterator n_it, bool pt)
{
    if (pt)
    {
        this->solve_PT(n_it, false); 
    }
    else
    {
        this->solve_VT(n_it, false); 
    }
    return this->is_root_type(); 
}
EoS::RootFlag EoS::is_root_type(double X, double T_, std::vector<double>& n_, int start_idx, bool pt) 
{ 
	if (pt)
	{
		this->solve_PT(X, T_, n_, start_idx, false);
    }
    else
    {
        this->solve_VT(X, T_, n_, start_idx, false);
    }
    return this->is_root_type(); 
}

std::vector<double> EoS::lnphi()
{
    // Return vector of lnphi for each component
    std::vector<double> ln_phi(ns);
    for (int i = 0; i < ns; i++)
    {
        ln_phi[i] = this->lnphii(i);
    }
    return ln_phi;
}
std::vector<double> EoS::lnphi(double p_, double T_, std::vector<double>& n_)
{
    this->solve_PT(p_, T_, n_, 0, false);
    return this->lnphi();
}
std::vector<double> EoS::dlnphi_dn()
{
    for (int i = 0; i < ns; i++)
    {
        for (int j = 0; j < ns; j++)
        {
            dlnphidn[i * ns + j] = this->dlnphii_dnj(i, j);
        }
    }
    return dlnphidn;
}
std::vector<double> EoS::dlnphi_dP()
{
    for (int i = 0; i < ns; i++)
    {
        dlnphidP[i] = this->dlnphii_dP(i);
    }
    return dlnphidP;
}
std::vector<double> EoS::dlnphi_dT()
{
    for (int i = 0; i < ns; i++)
    {
        dlnphidT[i] = this->dlnphii_dT(i);
    }
    return dlnphidT;
}
std::vector<double> EoS::d2lnphi_dT2()
{
    for (int i = 0; i < ns; i++)
    {
        d2lnphidT2[i] = this->d2lnphii_dT2(i);
    }
    return d2lnphidT2;
}
std::vector<double> EoS::d2lnphi_dTdn()
{
    for (int i = 0; i < ns; i++)
    {
        for (int j = 0; j < ns; j++)
        {
            d2lnphidTdn[i*ns + j] = this->d2lnphii_dTdnj(i, j);
        }
    }
    return d2lnphidTdn;
}

std::vector<double> EoS::fugacity(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Calculate component fugacities in mixture
    if (this->eos_in_range(n_.begin() + start_idx))
    {
	    // Calculate fugacity coefficients and fugacity
        std::vector<double> fi(ns);
        double nT_inv = 1./std::accumulate(n_.begin() + start_idx, n_.begin() + start_idx + ns, 0.);

        if (pt)
        {
            this->solve_PT(X, T_, n_, start_idx, false);
        }
        else
        {
            this->solve_VT(X, T_, n_, start_idx, false);
        }
        
        for (int i = 0; i < nc; i++)
        {
            fi[i] = std::exp(this->lnphii(i)) * n_[start_idx + i] * nT_inv * this->p;
        }
        return fi;
    }
    else
    {
        return std::vector<double>(ns, NAN);
    }
}

Eigen::MatrixXd EoS::calc_hessian(std::vector<double>::iterator n_it)
{
    // Calculate dlnphii/dnj
    this->solve_PT(n_it);
    this->dlnphi_dn();

    // Construct Hessian matrix of Gibbs energy surface
    Eigen::MatrixXd H_ = Eigen::MatrixXd(nc, nc);
    for (int j = 0; j < nc; j++)
    {
        for (int i = j; i < nc; i++)
        {
            H_(i, j) = dlnphidn[i*nc + j];  // PHI contribution
            H_(j, i) = dlnphidn[i*nc + j];  // PHI contribution
        }
        H_(j, j) += 1. / *(n_it + j);  // U contribution
    }

    return H_;
}
bool EoS::is_convex(std::vector<double>::iterator n_it)
{
    // Check if GE/TPD surface at composition n is convex
    Eigen::MatrixXd H_ = this->calc_hessian(n_it);

    // Perform Cholesky decomposition; only possible with positive definite matrix
    // Positive definite Hessian matrix corresponds to local minimum
    Eigen::LLT<Eigen::MatrixXd> lltOfA(H_);
    if(lltOfA.info() == Eigen::NumericalIssue)
    {
        return false;
    }
    else
    {
        return true;
    }
}
bool EoS::is_convex(double p_, double T_, std::vector<double>& n_, int start_idx)
{
    // Check if GE/TPD surface at P, T and composition n is convex
    this->init_PT(p_, T_);
    return this->is_convex(n_.begin() + start_idx);
}
double EoS::calc_condition_number(double p_, double T_, std::vector<double>& n_, int start_idx)
{
    if (this->eos_in_range(n_.begin() + start_idx))
    {
        // Calculate condition number of Hessian to evaluate curvature of GE/TPD surface at P, T and composition n
        this->init_PT(p_, T_);
        Eigen::MatrixXd H_ = this->calc_hessian(n_.begin() + start_idx);

        // Get the condition number for stability if Newton was used
        Eigen::VectorXd eigen = H_.eigenvalues().real();
        double min_eigen = *std::min_element(eigen.begin(), eigen.end());
        double max_eigen = *std::max_element(eigen.begin(), eigen.end());
        return (min_eigen > 0.) ? max_eigen/min_eigen : NAN;
    }
    else
    {
        return NAN;
    }
}

double EoS::cpi(double T_, int i)
{
    // cpi/R: Ideal gas heat capacity at constant pressure of component i
    return this->compdata.cpi[i][0]
         + this->compdata.cpi[i][1] * T_
         + this->compdata.cpi[i][2] * std::pow(T_, 2)
         + this->compdata.cpi[i][3] * std::pow(T_, 3);
}
double EoS::Cpi(double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Cpi/R: Ideal gas heat capacity at constant pressure
    (void) pt;
    double Cpi_ = 0.;
    for (int i = 0; i < ns; i++)
    {
        Cpi_ += n_[start_idx+i] * this->cpi(T_, i);
    }
    return Cpi_;
}
double EoS::cvi(double T_, int i)
{
    // cvi/R: Ideal gas heat capacity at constant volume of component i
    return this->cpi(T_, i) - 1.;
}
double EoS::Cvi(double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Cvi/R: Ideal gas heat capacity at constant volume
    (void) pt;
    double Cvi_ = 0.;
    for (int i = 0; i < ns; i++)
    {
        Cvi_ += n_[start_idx+i] * this->cvi(T_, i);
    }
    return Cvi_;
}
double EoS::hi(double T_, int i)
{
    // hi/R: Ideal gas enthalpy of component i
    return this->compdata.cpi[i][0] * (T_-this->compdata.T_0)
            + 1. / 2 * this->compdata.cpi[i][1] * (std::pow(T_, 2)-std::pow(this->compdata.T_0, 2))
            + 1. / 3 * this->compdata.cpi[i][2] * (std::pow(T_, 3)-std::pow(this->compdata.T_0, 3))
            + 1. / 4 * this->compdata.cpi[i][3] * (std::pow(T_, 4)-std::pow(this->compdata.T_0, 4));  // hi/R
}
double EoS::Hi(double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Hi/R: Ideal gas enthalpy
    (void) pt;
    double Hi_ = 0.;
    for (int i = 0; i < ns; i++)
    {
        Hi_ += n_[start_idx+i] * this->hi(T_, i);
    }
    return Hi_;
}
double EoS::dHi_dT(double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // 1/R dHi/dT = Cpi/R
    return this->Cpi(T_, n_, start_idx, pt);
}
std::vector<double> EoS::dHi_dni(double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Derivatives of ideal gas enthalpy with respect to composition
    (void) pt;
    std::vector<double> dHidn(ns);
    for (int i = 0; i < ns; i++)
    {
        dHidn[i] = (n_[start_idx+i]/this->N > 0.) ? this->hi(T_, i) : 0.;
    }
    return dHidn;  // 1/R dHi/dn
}
std::vector<double> EoS::d2Hi_dTdni(double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Derivatives of ideal gas enthalpy with respect to composition
    (void) pt;
    std::vector<double> d2HidTdn(ns);
    for (int i = 0; i < ns; i++)
    {
        d2HidTdn[i] = (n_[start_idx+i]/this->N > 0.) ? this->cpi(T_, i) : 0.;
    }
    return d2HidTdn;  // 1/R dHi/dn
}
double EoS::si(double X, double T_, int i, bool pt)
{
    if (pt)
    {   // si(P,T)/R: ideal gas entropy at constant pressure of component i
        return this->compdata.cpi[i][0] * std::log(T_ / this->compdata.T_0)
             + this->compdata.cpi[i][1] * (T_ - this->compdata.T_0)
             + 1. / 2 * this->compdata.cpi[i][2] * (std::pow(T_, 2)-std::pow(this->compdata.T_0, 2))
             + 1. / 3 * this->compdata.cpi[i][3] * (std::pow(T_, 3)-std::pow(this->compdata.T_0, 3))
             - std::log(X / this->compdata.P_0);
    }
    else
    {   // si(V,T)/R: ideal gas entropy at constant volume of component i
        return this->compdata.cpi[i][0] * std::log(T_ / this->compdata.T_0)
             + this->compdata.cpi[i][1] * (T_ - this->compdata.T_0)
             + 1. / 2 * this->compdata.cpi[i][2] * (std::pow(T_, 2)-std::pow(this->compdata.T_0, 2))
             + 1. / 3 * this->compdata.cpi[i][3] * (std::pow(T_, 3)-std::pow(this->compdata.T_0, 3))
             - std::log(T_ / this->compdata.T_0) - std::log(X/N / this->compdata.V_0);
    }
}
double EoS::dsi_dT(double X, double T_, int i, bool pt)
{
    (void) X;
    if (pt)
    {   // 1/R dsi(P,T)/dT = cpi/T
        return this->cpi(T_, i) / T_;
    }
    else
    {   // 1/R dsi(V,T)/dT = cvi/T
        return this->cvi(T_, i) / T_;
    }
}
double EoS::Si(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Ideal gas entropy
    double Si_ = 0.;
    double nT_inv = 1./std::accumulate(n_.begin() + start_idx, n_.begin() + start_idx + ns, 0.);
    for (int i = 0; i < ns; i++)
    {
        double xi = n_[start_idx+i] * nT_inv;
        if (xi > 0.)
        {
            Si_ += n_[start_idx+i] * (this->si(X, T_, i, pt) - std::log(xi));
        }
    }
    return Si_;  // Si/R
}
double EoS::dSi_dT(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Ideal gas entropy
    double Si_ = 0.;
    double nT_inv = 1./std::accumulate(n_.begin() + start_idx, n_.begin() + start_idx + ns, 0.);
    for (int i = 0; i < ns; i++)
    {
        double xi = n_[start_idx+i] * nT_inv;
        if (xi > 0.)
        {
            Si_ += n_[start_idx+i] * this->dsi_dT(X, T_, i, pt);
        }
    }
    return Si_;  // Si/R
}
std::vector<double> EoS::dSi_dni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Calculate partial molar ideal gas entropy
    std::vector<double> dSidn(ns, 0.);

    double nT_inv = 1./std::accumulate(n_.begin() + start_idx, n_.begin() + start_idx + ns, 0.);
    for (int i = 0; i < ns; i++)
    {
        // Ideal contribution to Gibbs free energy
        double xi = n_[start_idx+i] * nT_inv;

        if (xi > 0.)
        {
            dSidn[i] += this->si(X, T_, i, pt) - std::log(xi);
            for (int j = 0; j < ns; j++)
            {
                if (i == j)
                {
                    dSidn[i] -= n_[start_idx + i] / xi * (nT_inv - n_[start_idx + j] * std::pow(nT_inv, 2));
                }
                else
                {
                    dSidn[i] -= n_[start_idx + i] / xi * -n_[start_idx + j] * std::pow(nT_inv, 2);
                }
            }
        }
    }
    return dSidn;  // 1/R dSi/dn
}
std::vector<double> EoS::d2Si_dTdni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Calculate partial molar ideal gas entropy
    std::vector<double> d2SidTdn(ns);

    for (int i = 0; i < ns; i++)
    {
        // Ideal contribution to Gibbs free energy
        double xi = n_[start_idx+i] / this->N;
        d2SidTdn[i] = (xi > 0) ? this->dsi_dT(X, T_, i, pt) : 0.;
    }
    return d2SidTdn;  // 1/R d2Si/dTdni
}
double EoS::Gi(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Ideal gas Gibbs free energy
    return this->Hi(T_, n_, start_idx, pt) - T_ * this->Si(X, T_, n_, start_idx, pt);  // Gi/R
}
double EoS::dGi_dT(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    return this->dHi_dT(T_, n_, start_idx, pt) - this->Si(X, T_, n_, start_idx, pt) - T_ * this->dSi_dT(X, T_, n_, start_idx, pt);  // 1/R dGi/dT
}
std::vector<double> EoS::dGi_dni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Derivative of ideal gas Gibbs free energy with respect to composition
    std::vector<double> dHidn = this->dHi_dni(T_, n_, start_idx, pt);
    std::vector<double> dSidn = this->dSi_dni(X, T_, n_, start_idx, pt);
    std::vector<double> dGidn(ns, 0.);

    for (int i = 0; i < ns; i++)
    {
        if (n_[start_idx + i]/this->N > 0.)
        {
            dGidn[i] = dHidn[i] - T_ * dSidn[i];
        }
    }
    return dGidn;  // 1/R dGi/dni
}
std::vector<double> EoS::d2Gi_dTdni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Calculate partial molar ideal gas entropy
    std::vector<double> d2GidTdn(ns);

    for (int i = 0; i < ns; i++)
    {
        // Ideal contribution to Gibbs free energy
        double xi = n_[start_idx+i] / this->N;
        d2GidTdn[i] = (xi > 0) ? this->cpi(T_, i) - this->si(X, T_, i, pt) - T_ * this->dsi_dT(X, T_, i, pt) : 0.;
    }
    return d2GidTdn;  // 1/R d2Gi/dTdni
}
/*
double EoS::Ai(double X, double T_, std::vector<double>& n_, bool pt)
{
    // Ideal gas Helmholtz free energy
    return this->Gi(X, T_, n_, pt) + PV;
}
double EoS::Ui(double X, double T_, std::vector<double>& n_, bool pt)
{

}
*/

double EoS::Sr(double X, double T_, std::vector<double>& n_, int start_idx, bool pt) 
{
    // Residual entropy
    if (pt)
    {
        // Sr(PT)/R = Hr(PT)/RT - Gr(PT)/RT
        double Gr = this->Gr(X, T_, n_, start_idx, pt);  // Gr/RT
        double Hr = this->Hr(X, T_, n_, start_idx, pt);  // Hr/RT
        return Hr - Gr;  // Sr/R = Hr/RT - Gr/RT
    }
    else
    {
        // Sr(VT)/R = Ur(VT)/RT - Ar(VT)/RT
        // double Ar = this->Ar(X, T_, n_, start_idx, pt);  // Ar/RT
        // double Ur = this->Ur(X, T_, n_, start_idx, pt);  // Ur/RT
        // return Ur - Ar;  // Sr/R = Ur/RT - Ar/RT
        return NAN;
    }
}
double EoS::Gr(double X, double T_, std::vector<double>& n_, int start_idx, bool pt) 
{
    // Residual Gibbs free energy
    if (pt)
    {
        // Gr(PT)/RT = Σi ni lnphii
        if (this->eos_in_range(n_.begin() + start_idx))
        {
            this->solve_PT(X, T_, n_, start_idx, false);

            double Gr = 0.;
            bool nans = false;
            for (int i = 0; i < ns; i++)
            {
                if (!std::isnan(this->lnphii(i)))
                {
                    Gr += n_[start_idx + i] * this->lnphii(i);  // partial molar Gibbs energy
                }
                else
                {
                    nans = true;
                }
            }
            // If all NANs, G will be equal to zero, so return NAN; else return G
            return (nans && Gr == 0.) ? NAN : Gr;  // Gr/RT
        }
        else
        {
            // Outside of EoS range
            return NAN;
        }
    }
    else
    {
        return NAN;
    }
}
double EoS::Hr(double X, double T_, std::vector<double>& n_, int start_idx, bool pt) 
{
    // Residual enthalpy
    if (pt)
    {
        this->solve_PT(X, T_, n_, start_idx, true);
        dlnphidT = this->dlnphi_dT();

        double Hr = 0.;
        bool nans = false;
        for (int i = 0; i < ns; i++)
        {
            if (!std::isnan(dlnphidT[i]))
            {
                Hr -= n_[start_idx + i] * T * dlnphidT[i];  // partial molar enthalpy Hri/RT
            }
            else
            {
                nans = true;
            }
        }
        // If all NANs, H will be equal to zero, so return NAN; else return H
        return (nans && Hr == 0.) ? NAN : Hr;  // Hr/RT
    }
    else
    {
        return NAN;
    }
}
/*
double EoS::Ar(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Calculate residual Helmholtz free energy of mixture at P, T, x
}
double EoS::Ur(double p_, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Calculate residual internal energy of mixture at P, T, x
}
*/

std::vector<double> EoS::dSr_dni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Derivatives of residual entropy with respect to composition
    std::vector<double> dSrdni(ns);
    
    if (pt)
    {
        // Sr(PT)/R = Hr(PT)/RT - Gr(PT)/RT
        std::vector<double> dGrdni = this->dGr_dni(X, T_, n_, start_idx, pt);  // d(Gr/RT)/dni
        std::vector<double> dHrdni = this->dHr_dni(X, T_, n_, start_idx, pt);  // d(Hr/RT)/dni

        for (int i = 0; i < ns; i++)
        {
            dSrdni[i] = dHrdni[i] - dGrdni[i];
        }
        return dSrdni;  // d(Sr/R)/dni = d(Hr/RT)/dni - d(Gr/RT)/dni
    }
    else
    {
        // Sr(VT)/R = Ur(VT)/RT - Ar(VT)/RT
        // std::vector<double> dArdni = this->dAr_dni(X, T_, n_, start_idx, pt);  // d(Ar/RT)/dni
        // std::vector<double> dUrdni = this->dUr_dni(X, T_, n_, start_idx, pt);  // d(Ur/RT)/dni

        // for (int i = 0; i < nc; i++)
        // {
        //     dSrdni[i] = dUrdni[i] - dArdni[i];
        // }
        return dSrdni;  // d(Sr/R)/dni = d(Ur/RT)/dni - d(Ar/RT)/dni
    }
}
std::vector<double> EoS::dGr_dni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Derivatives of residual Gibbs free energy with respect to composition
    if (pt)
    {
        // Gr(PT)/RT = Σi ni lnphii
        std::vector<double> dGrdni(ns, 0.);
        if (this->eos_in_range(n_.begin() + start_idx))
        {
            this->solve_PT(X, T_, n_, start_idx, true);
            std::vector<double> ln_phi = this->lnphi();
            this->dlnphidn = this->dlnphi_dn();

            for (int i = 0; i < ns; i++)
            {
                dGrdni[i] += ln_phi[i];  // partial molar Gibbs energy
                for (int j = 0; j < ns; j++)
                {
                    dGrdni[j] += n_[start_idx + i] * this->dlnphidn[i * ns + j];
                }
            }
            return dGrdni;  // d(Gr/RT)/dni
        }
        else
        {
            // Outside of EoS range
            return std::vector<double>(ns, NAN);
        }
    }
    else
    {
        return std::vector<double>(ns, NAN);
    }
}
std::vector<double> EoS::dHr_dni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Calculate partial molar enthalpy at P,T,n
    std::vector<double> dHrdn(ns, 0.);

    if (pt)
    {
        this->solve_PT(X, T_, n_, start_idx, true);
    }
    else
    {
        this->solve_VT(X, T_, n_, start_idx, true);
    }

    dlnphidT = this->dlnphi_dT();
    d2lnphidTdn = this->d2lnphi_dTdn();
    for (int j = 0; j < ns; j++)
    {
        dHrdn[j] -= T * dlnphidT[j];
        for (int i = 0; i < ns; i++)
        {
            dHrdn[i] -= n_[start_idx + j] * T * d2lnphidTdn[j*ns + i];
        }
    }
    return dHrdn;
}
/*
double EoS::dAr_dni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Ideal + residual Helmholtz free energy
    double Ai = this->ideal.A(T_, n_.begin() + start_idx);
    double Ar = this->Ar(X, T_, n_, start_idx, pt);
    return Ai + Ar;
}
double EoS::dUr_dni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Ideal + residual internal energy
    double Ui = this->ideal.U(T_, n_.begin() + start_idx);
    double Ur = this->Ur(X, T_, n_, start_idx, pt);
    return Ui + Ur;
}
*/

double EoS::S(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Ideal + residual entropy
    if (pt)
    {
        this->solve_PT(X, T_, n_, start_idx, true);
    }
    else
    {
        this->solve_VT(X, T_, n_, start_idx, true);
    }
    return this->Si(X, T_, n_, start_idx, pt) + this->Sr(X, T_, n_, start_idx, pt);  // S/R
}
double EoS::G(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Ideal + residual Gibbs free energy
    if (pt)
    {
        this->solve_PT(X, T_, n_, start_idx, true);
    }
    else
    {
        this->solve_VT(X, T_, n_, start_idx, true);
    }
    return this->Gi(X, T_, n_, start_idx, pt) + this->Gr(X, T_, n_, start_idx, pt) * this->T;  // G/R = Gi/R + Gr/RT * T
}
double EoS::H(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Ideal + residual enthalpy
    if (pt)
    {
        this->solve_PT(X, T_, n_, start_idx, true);
    }
    else
    {
        this->solve_VT(X, T_, n_, start_idx, true);
    }
    return this->Hi(T_, n_, start_idx, pt) + this->Hr(X, T_, n_, start_idx, pt) * T_;  // H/R = Hi/R + Hr/RT * T
}
/*
double EoS::A(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Ideal + residual Helmholtz free energy
    double Ai = this->ideal.A(T_, n_.begin() + start_idx);
    double Ar = this->Ar(X, T_, n_, start_idx, pt);
    return Ai + Ar;
}
double EoS::U(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Ideal + residual internal energy
    double Ui = this->ideal.U(T_, n_.begin() + start_idx);
    double Ur = this->Ur(X, T_, n_, start_idx, pt);
    return Ui + Ur;
}
*/

std::vector<double> EoS::dS_dni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Derivatives of total entropy with respect to composition
    std::vector<double> dSidn = this->dSi_dni(X, T_, n_, start_idx, pt);
    std::vector<double> dSrdn = this->dSr_dni(X, T_, n_, start_idx, pt);
    std::vector<double> dSdn(ns);

    for (int i = 0; i < ns; i++)
    {
        dSdn[i] = dSidn[i] + dSrdn[i];
    }
    return dSdn;  // 1/R dS/dni
}
std::vector<double> EoS::dG_dni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Derivatives of total Gibbs free energy with respect to composition
    std::vector<double> dGidn = this->dGi_dni(X, T_, n_, start_idx, pt);
    std::vector<double> dGrdn = this->dGr_dni(X, T_, n_, start_idx, pt);
    std::vector<double> dGdn(ns);

    for (int i = 0; i < ns; i++)
    {
        dGdn[i] = dGidn[i] + dGrdn[i] * T_;
    }
    return dGdn;  // 1/R dG/dni
}
std::vector<double> EoS::dH_dni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Derivatives of total enthalpy with respect to composition
    std::vector<double> dHdn = this->dHr_dni(X, T_, n_, start_idx, pt);

    for (int i = 0; i < ns; i++)
    {
        dHdn[i] *= T_;

        // Ideal contribution to enthalpy
        dHdn[i] += this->hi(T_, i);
    }
    return dHdn;  // 1/R dH/dni
}
/*
double EoS::dA_dni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Ideal + residual Helmholtz free energy
    double Ai = this->ideal.A(T_, n_.begin() + start_idx);
    double Ar = this->Ar(X, T_, n_, start_idx, pt);
    return Ai + Ar;
}
double EoS::dU_dni(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Ideal + residual internal energy
    double Ui = this->ideal.U(T_, n_.begin() + start_idx);
    double Ur = this->Ur(X, T_, n_, start_idx, pt);
    return Ui + Ur;
}
*/

double EoS::Cpr(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Calculate residual heat capacity at constant pressure: Cp/R = 1/R (dHr/dT)_P
    if (pt)
    {
        this->solve_PT(X, T_, n_, start_idx, true);
    }
    else
    {
        this->solve_VT(X, T_, n_, start_idx, true);
    }
    this->dlnphidT = this->dlnphi_dT();
    this->d2lnphidT2 = this->d2lnphi_dT2();

    // Calculate derivative of residual enthalpy Hr w.r.t T
    double cpr = 0.;
    for (int i = 0; i < this->ns; i++)
    {
        cpr -= n_[start_idx + i] * T_ * (2. * this->dlnphidT[i] + T_ * this->d2lnphidT2[i]);
    }
    return cpr;  // Cpr/R
}
double EoS::Cp(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Calculate heat capacity at constant pressure: Cp = (dH/dT)_P
    // Calculate derivative of enthalpy Hr w.r.t T and add ideal part of heat capacity
    double cpr = this->Cpr(X, T_, n_, start_idx, pt);
    double cpi = 0.;
    for (int i = 0; i < this->ns; i++)
    {
        cpi += n_[start_idx + i] * this->cpi(T_, i);
    }
    return cpi + cpr;  // Cp/R
}
double EoS::Cv(double X, double T_, std::vector<double>& n_, int start_idx, bool pt)
{
    // Total heat capacity at constant volume Cv/R
    double cvr = this->Cvr(X, T_, n_, start_idx, pt);
    double cvi = 0.;
    for (int i = 0; i < ns; i++)
    {
        cvi += n_[start_idx + i] * this->cvi(T, i);
    }
    return cvi + cvr;
}

std::vector<double> EoS::G_pure(double X, double T_, bool pt)
{
    // Calculate pure component Gibbs energy at P, T
    if (pt && (X != p || T_ != T))
    {
        this->init_PT(X, T_);
        gpure = std::vector<double>(nc, NAN);

        // Loop over components
        for (int i = 0; i < nc; i++)
        {
            std::vector<double> n_(nc, 0.);
            n_[i] = 1.;

            if (this->eos_in_range(n_.begin()))
            {
                this->solve_PT(n_.begin(), false);

                if (this->select_root(n_.begin()))
                {
                    gpure[i] = this->lnphii(i);
                }
                else
                {
                    gpure[i] = NAN;
                }
            }
        }
    }
    return gpure;
}

std::vector<double> EoS::mix_min(std::vector<double>& x_, std::vector<double>& gpure_, double& gmix)
{
    // Find local minimum of Gibbs energy of mixing
    // Last component is dependent variable s.t. Σi ni - 1 = 0
    std::vector<double> xx = x_;

    double obj = 0.;
    Eigen::VectorXd g(ns-1), X(ns-1);
    Eigen::MatrixXd H_(ns-1, ns-1);

    // Initialize solution vector
    for (int i = 0; i < ns-1; i++)
    {
        X(i) = xx[i];
    }

    // Solve fugacity coefficient and derivatives
    this->solve_PT(xx.begin(), false);
    double sumX = std::accumulate(X.data(), X.data() + X.size(), 0.);
    double ln_1x = std::log(1.-sumX);
    obj = (1.-sumX) * (this->lnphii(ns-1) + ln_1x - gpure_[ns-1]);
    for (int i = 0; i < ns-1; i++)
    {
        obj += X[i] * (this->lnphii(i) + std::log(X[i]) - gpure_[i]);
    }
    double obj_init = obj;

    // Start Newton loop
    int iter = 1;
    while (iter < 10)
    {
        // Entries for xj
        this->solve_PT(xx.begin(), true);
        for (int j = 0; j < ns-1; j++)
        {
            // Gradient vector for xj entries dG/dxj
            g(j) = this->lnphii(j) + std::log(X[j]) - gpure_[j] - this->lnphii(ns-1) - ln_1x + gpure_[ns-1];

            for (int k = j; k < ns-1; k++)
            {
                // Hessian for d2G/dxjdxk
                double Hjk = this->dlnphii_dnj(j, k) - this->dlnphii_dnj(ns-1, k);
                H_(j, k) = Hjk;
                H_(k, j) = Hjk;
            }
            H_(j, j) += 1./xx[j] + 1./(1.-sumX);
        }

        if (g.squaredNorm() < 1e-3)
        {
            gmix = obj;
            return xx;
        }

        // Solve Newton step
        Eigen::LLT<Eigen::MatrixXd> lltOfH(H_);
        Eigen::VectorXd dX;
        if (lltOfH.info() == Eigen::NumericalIssue )
        {
            // print("LLT of H not possible", 1);
            gmix = (obj < obj_init) ? obj : obj_init;
            return (obj < obj_init) ? xx : x_;
        }
        else
        {
            dX = lltOfH.solve(g);
        }

        bool reduced = false;
        double lamb = 0.5;
        double obj_old = obj;
        Eigen::VectorXd X_old = X;
        int iter2 = 1;
        while (true)
        {
            X = X_old - lamb*dX;
            xx = std::vector<double>(X.data(), X.data() + X.size());
            sumX = std::accumulate(X.data(), X.data() + X.size(), 0.);
            ln_1x = std::log(1.-sumX);
            xx.push_back(1.-sumX);

            if (!std::count_if(xx.begin(), xx.end(), [](double xi) { return (xi < 0.); }))
            {
                // Calculate objective function
                this->solve_PT(xx.begin(), false);
                obj = (1.-sumX) * (this->lnphii(ns-1) + ln_1x - gpure_[ns-1]);
                for (int i = 0; i < ns-1; i++)
                {
                    obj += X[i] * (this->lnphii(i) + std::log(X[i]) - gpure_[i]);
                }

                if (obj < obj_old)
                {
                    reduced = true;
                }
            }

            if (iter2 == 10)
            {
                return (obj < obj_init) ? xx : x_;
            }
            else if (!reduced)
            {
                lamb *= 0.5;
                iter2++;
            }
            else
            {
                break;
            }
        }

        iter++;
    }
    gmix = obj;
    return xx;
}

double EoS::dxi_dnj(std::vector<double>& n_, int i, int j)
{
    // Derivative of mole fractions with respect to composition
    // xi = ni/nT
    double nT_inv = 1./std::accumulate(n_.begin(), n_.end(), 0.);
    if (i == j)
    {
        // dxi/dni = 1/nT - ni/nT^2
        return nT_inv - n_[i] * std::pow(nT_inv, 2);
    }
    else
    {   
        // dxi/dnj = -ni/nT^2
        return -n_[i] * std::pow(nT_inv, 2);
    }
}
double EoS::dxj_to_dnk(std::vector<double>& dlnphiidxj, std::vector<double>::iterator n_it, int k) {
    // Translate from dlnphii/dxj to dlnphii/dnk
	// dlnphii/dnk = 1/V * [dlnphii/dxk - sum_j xj dlnphii/dxj]
    double nT_inv = 1./std::accumulate(n_it, n_it + this->ns, 0.);
	double dlnphiidnk = dlnphiidxj[k];
	for (int j = 0; j < ns; j++)
	{
        double nj = *(n_it + j);
		dlnphiidnk -= nj * nT_inv * dlnphiidxj[j];
	}
	dlnphiidnk *= nT_inv;
    return dlnphiidnk;
}
std::vector<double> EoS::dxj_to_dnk(std::vector<double>& dlnphiidxj, std::vector<double>::iterator n_it)
{
    // Translate from dlnphii/dxj to dlnphii/dnk
	// dlnphii/dnj = 1/V * [dlnphii/dxj - Σk xk dlnphii/dxk]
    double nT_inv = 1./std::accumulate(n_it, n_it + this->ns, 0.);
    std::vector<double> dlnphiidnk(ns*ns);

    for (int i = 0; i < ns; i++)
    {
        double sum_xj = 0.;
        for (int j = 0; j < ns; j++)
        {
            double nj = *(n_it + j);
            sum_xj += nj * nT_inv * dlnphiidxj[i*ns + j];
        }
        for (int k = 0; k < ns; k++)
        {
            dlnphiidnk[i*ns + k] = nT_inv * (dlnphiidxj[i*ns + k] - sum_xj);
        }
    }
    return dlnphiidnk;
}

int EoS::dlnphi_test(double p_, double T_, std::vector<double>& n_, double tol, bool verbose) 
{
    // Analytical derivatives w.r.t. composition, pressure and temperature
	int error_output = 0;

    double p0 = p_;
    double T0 = T_;
    std::vector<double> n0 = n_;

	this->solve_PT(p0, T0, n0, 0, true);
	std::vector<double> lnphi0 = this->lnphi();
    std::vector<double> dlnphidn0 = this->dlnphi_dn();
	std::vector<double> dlnphidP0 = this->dlnphi_dP();
	std::vector<double> dlnphidT0 = this->dlnphi_dT();
	std::vector<double> d2lnphidT20 = this->d2lnphi_dT2();
    std::vector<double> d2lnphidTdn0 = this->d2lnphi_dTdn();

	// Calculate numerical derivatives w.r.t. composition
	double dn = 1e-4;
	std::vector<double> lnphi1, dlnphi1;
    std::vector<double> nn = n0;
    for (int j = 0; j < ns; j++)
	{
        // Numerical derivative of lnphi w.r.t. nj
        double dnj = n0[j] * dn;
	    nn[j] += dnj;
		this->solve_PT(p0, T0, nn, 0, false);
	    lnphi1 = this->lnphi();
		nn[j] -= dnj;
        for (int i = 0; i < ns; i++)
	    {
            double dlnphidn_an = dlnphidn0[i*ns+j];
            double dlnphi = (lnphi1[i] - lnphi0[i]);
	        double dlnphidn_num = dlnphi/dnj;
			double d = std::log(std::fabs(dlnphidn_num + 1e-15)) - std::log(std::fabs(dlnphidn_an + 1e-15));
            if (verbose || (std::fabs(d) > tol && (std::fabs(dlnphidn_an) > 1e-14 && std::fabs(dlnphidn_num) > 1e-14)))
	        {
	            print("comp", std::vector<int>{i, j});
        	    print("dlnphi/dn", std::vector<double>{dlnphidn_an, dlnphidn_num, d});
                error_output++;
	        }
        }
	}

	// Calculate numerical derivatives w.r.t. pressure
	double dp = 1e-6;
	this->solve_PT(p0 + dp, T0, n0, 0, false);
	lnphi1 = this->lnphi();

	// Compare analytical and numerical
	for (int i = 0; i < ns; i++)
	{
        double dlnphidP_num = (lnphi1[i] - lnphi0[i])/dp;
		// Use logarithmic scale to compare
		double d = std::log(std::fabs(dlnphidP_num + 1e-15)) - std::log(std::fabs(dlnphidP0[i] + 1e-15));
	    if (verbose || (std::fabs(d) > tol && (std::fabs(dlnphidP0[i]) > 1e-8 && std::fabs(dlnphidP_num) > 1e-8)))
        {
        	print("comp", i);
            print("dlnphi/dP", std::vector<double>{dlnphidP0[i], dlnphidP_num, d});
    	    error_output++;
	    }
    }

	// Calculate numerical derivatives w.r.t. temperature
	double dT = 1e-6;
	this->solve_PT(p0, T0 + dT, n0, 0, false);
	lnphi1 = this->lnphi();
	dlnphi1 = this->dlnphi_dT();

	// Compare analytical and numerical
	for (int i = 0; i < ns; i++)
	{
		double dlnphidT_num = (lnphi1[i] - lnphi0[i])/dT;
		double d = std::log(std::fabs(dlnphidT_num + 1e-15)) - std::log(std::fabs(dlnphidT0[i] + 1e-15));
		if (verbose || (std::fabs(d) > tol && (std::fabs(dlnphidT0[i]) > 1e-8 && std::fabs(dlnphidT_num) > 1e-8)))
		{
			print("comp", i);
			print("dlnphi/dT", std::vector<double>{dlnphidT0[i], dlnphidT_num, d});
			error_output++;
		}

        double dT2_an = d2lnphidT20[i];
		double d2lnphidT2_num = (dlnphi1[i] - dlnphidT0[i])/dT;
		d = std::log(std::fabs(d2lnphidT2_num + 1e-15)) - std::log(std::fabs(dT2_an + 1e-15));
		if (verbose || (std::fabs(d) > tol && (std::fabs(dT2_an) > 1e-8 && std::fabs(d2lnphidT2_num) > 1e-8)))
		{
			print("comp", i);
			print("d2lnphi/dT2", std::vector<double>{dT2_an, d2lnphidT2_num, d});
			error_output++;
		}

        // Second derivative of d2lnphi/dTdnj w.r.t. temperature and composition
        this->solve_PT(p0, T0 + dT, n0, 0, true);
        std::vector<double> dlnphidn1 = this->dlnphi_dn();
        for (int j = 0; j < ns; j++)
        {
            double dTdni_an = d2lnphidTdn0[i*ns+j];
            double dTdni = (dlnphidn1[i*ns+j] - dlnphidn0[i*ns+j]);
            double dTdni_num = dTdni / dT;
            d = std::log(std::fabs(dTdni_an + 1e-15)) - std::log(std::fabs(dTdni_num + 1e-15));
            if (verbose || (std::fabs(d) > tol && (std::fabs(dT2_an) > 1e-14 && std::fabs(dTdni) > 1e-14)))
            {
                print("i, j", {i, j});
                print("d2lnphi/dTdnj != d2lnphi/dTdnj", std::vector<double>{dTdni_an, dTdni_num, d});
                error_output++;
            }
        }
	}

	// Calculate lnphi and derivatives for -n
	std::vector<double> n_n = n0;
	std::transform(n0.begin(), n0.end(), n_n.begin(), [](double element) { return element *= -1; });
	this->solve_PT(p0, T0, n_n, 0, true);
	std::vector<double> lnphi_n = this->lnphi();
	std::vector<double> dlnphidP_n = this->dlnphi_dP();
	std::vector<double> dlnphidT_n = this->dlnphi_dT();
	std::vector<double> d2lnphidT2_n = this->d2lnphi_dT2();
	std::vector<double> dlnphidn_n = this->dlnphi_dn();
    std::vector<double> d2lnphidTdn_n = this->d2lnphi_dTdn();

	// Compare positive and negative
	for (int i = 0; i < ns; i++)
	{
		double d = (lnphi0[i] - lnphi_n[i])/lnphi0[i];
		if (verbose || std::fabs(d) > tol || (std::isnan(d) && lnphi0[i] > 0.))
		{
			print("comp", i);
			print("lnphi negative", std::vector<double>{lnphi0[i], lnphi_n[i], d});
			error_output++;
		}
		d = (dlnphidP0[i] - dlnphidP_n[i])/dlnphidP0[i];
		if (verbose || std::fabs(d) > tol || (std::isnan(d) && std::fabs(dlnphidP0[i]) > 0.))
		{
			print("comp", i);
			print("dlnphi/dP negative", std::vector<double>{dlnphidP0[i], dlnphidP_n[i], d});
			error_output++;
		}
		d = (dlnphidT0[i] - dlnphidT_n[i])/dlnphidT0[i];
		if (verbose || std::fabs(d) > tol || (std::isnan(d) && std::fabs(dlnphidT0[i]) > 0.))
		{
			print("comp", i);
			print("dlnphi/dT negative", std::vector<double>{dlnphidT0[i], dlnphidT_n[i], d});
			error_output++;
		}
		d = (d2lnphidT20[i] - d2lnphidT2_n[i])/d2lnphidT20[i];
		if (verbose || std::fabs(d) > tol || (std::isnan(d) && std::fabs(d2lnphidT20[i]) > 0.))
		{
			print("comp", i);
			print("d2lnphi/dT2 negative", std::vector<double>{d2lnphidT20[i], d2lnphidT2_n[i], d});
			error_output++;
		}
		for (int j = 0; j < ns; j++)
		{
			d = (dlnphidn0[i*ns + j] + dlnphidn_n[i*ns + j])/dlnphidn0[i*ns + j];
			if (verbose || std::fabs(d) > tol || (std::isnan(d) && std::fabs(dlnphidn0[i*ns+j]) > 0.))
			{
				print("comp", std::vector<int>{i, j});
				print("dlnphi/dn negative", std::vector<double>{dlnphidn0[i*ns + j], dlnphidn_n[i*ns + j], d});
				error_output++;
			}
		}
        for (int j = 0; j < ns; j++)
		{
			d = dlnphidn0[i*ns + j] + dlnphidn_n[i*ns + j];
			if (verbose || std::fabs(d) > tol || std::isnan(d))
			{
				print("comp", std::vector<int>{i, j});
				print("d2lnphi/dTdn negative", std::vector<double>{d2lnphidTdn0[i*ns + j], d2lnphidTdn_n[i*ns + j], d});
				error_output++;
			}
		}
	}

    return error_output;
}

int EoS::properties_test(double p_, double T_, std::vector<double>& n_, double tol, bool verbose)
{
    // Analytical derivatives w.r.t. composition, pressure and temperature
	int error_output = 0;

    double p0 = p_;
    double T0 = T_;
    std::vector<double> n0 = n_;

    // Thermal properties tests
    bool pt = true;
    double d, dT = T0 * 1e-6;
    double Hi0, Hi1, Hr0, Hr1, H0, H1, Gi0, Gi1, Gr0, Gr1, G0, G1, Si0, Si1, Sr0, Sr1, S0, S1;
    std::vector<double> dHi0, dHr0, dH0, dGi0, dGr0, dG0, dSi0, dSr0, dS0;

    // Ideal gas heat capacity at constant pressure Cp
    Hi0 = this->Hi(T0, n0, 0, pt);
    double CPi = this->Cpi(T0, n0, 0, pt);
    Hi1 = this->Hi(T0+dT, n0, 0, pt);
    double CPi_num = (Hi1-Hi0)/dT;
    d = std::log(std::fabs(CPi + 1e-15)) - std::log(std::fabs(CPi_num + 1e-15));
    if (verbose || std::fabs(d) > tol) { print("Cpi", std::vector<double>{CPi, CPi_num, d}); error_output++; }

    // Residual heat capacity at constant pressure Cpr
    Hr0 = this->Hr(p0, T0, n0, 0, pt) * T0;
    double CPr = this->Cpr(p0, T0, n0, 0, pt);
    Hr1 = this->Hr(p0, T0+dT, n0, 0, pt) * (T0+dT);
    double CPr_num = (Hr1-Hr0)/dT;
    d = std::log(std::fabs(CPr + 1e-15)) - std::log(std::fabs(CPr_num + 1e-15));
    if (verbose || std::fabs(d) > tol) { print("Cpr", std::vector<double>{CPr, CPr_num, d}); error_output++; }

    // Heat capacity at constant pressure Cp
    H0 = this->H(p0, T0, n0, 0, pt);
    double CP = this->Cp(p0, T0, n0, 0, pt);
    H1 = this->H(p0, T0+dT, n0, 0, pt);
    double CP_num = (H1-H0)/dT;
    d = std::log(std::fabs(CP + 1e-15)) - std::log(std::fabs(CP_num + 1e-15));
    if (verbose || std::fabs(d) > tol) { print("Cp", std::vector<double>{CP, CP_num, d}); error_output++; }

    // Derivative of (residual) Gibbs energy, enthalpy and entropy w.r.t. pressure, temperature and composition
    this->solve_PT(p0, T0, n0, 0, true);
    Gi0 = this->Gi(p0, T0, n0, 0, pt);
    Gr0 = this->Gr(p0, T0, n0, 0, pt);
    G0 = this->G(p0, T0, n0, 0, pt);
    Hi0 = this->Hi(T0, n0, 0, pt);
    Hr0 = this->Hr(p0, T0, n0, 0, pt);
    H0 = this->H(p0, T0, n0, 0, pt);
    Si0 = this->Si(p0, T0, n0, 0, pt);
    Sr0 = this->Sr(p0, T0, n0, 0, pt);
    S0 = this->S(p0, T0, n0, 0, pt);

    // double dGidP = this->dGi_dP(p0, T0, n0, 0, pt);
    // double dGrdP = this->dGr_dP(p0, T0, n0, 0, pt);
    // double dGdP = this->dG_dP(p0, T0, n0, 0, pt);
    // double dHidP = this->dHi_dP(T0, n0, 0, pt);
    // double dHrdP = this->dHr_dP(p0, T0, n0, 0, pt);
    // double dHdP = this->dH_dP(p0, T0, n0, 0, pt);
    // double dSidP = this->dSi_dP(p0, T0, n0, 0, pt);
    // double dSrdP = this->dSr_dP(p0, T0, n0, 0, pt);
    // double dSdP = this->dS_dP(p0, T0, n0, 0, pt);

    double dGidT = this->dGi_dT(p0, T0, n0, 0, pt);
    // double dGrdT = this->dGr_dT(p0, T0, n0, 0, pt);
    // double dGdT = this->dG_dT(p0, T0, n0, 0, pt);
    double dHidT = this->dHi_dT(T0, n0, 0, pt);
    // double dHrdT = this->dHr_dT(p0, T0, n0, 0, pt);
    // double dHdT = this->dH_dT(p0, T0, n0, 0, pt);
    double dSidT = this->dSi_dT(p0, T0, n0, 0, pt);
    // double dSrdT = this->dSr_dT(p0, T0, n0, 0, pt);
    // double dSdT = this->dS_dT(p0, T0, n0, 0, pt);

    std::vector<double> dGidn = this->dGi_dni(p0, T0, n0, 0, pt);
    std::vector<double> dGrdn = this->dGr_dni(p0, T0, n0, 0, pt);
    std::vector<double> dGdn = this->dG_dni(p0, T0, n0, 0, pt);
    std::vector<double> dHidn = this->dHi_dni(T0, n0, 0, pt);
    std::vector<double> dHrdn = this->dHr_dni(p0, T0, n0, 0, pt);
    std::vector<double> dHdn = this->dH_dni(p0, T0, n0, 0, pt);
    std::vector<double> dSidn = this->dSi_dni(p0, T0, n0, 0, pt);
    std::vector<double> dSrdn = this->dSr_dni(p0, T0, n0, 0, pt);
    std::vector<double> dSdn = this->dS_dni(p0, T0, n0, 0, pt);

    std::vector<double> d2GidTdn = this->d2Gi_dTdni(p0, T0, n0, 0, pt);
    // std::vector<double> d2GrdTdn = this->d2Gr_dTdni(p0, T0, n0, 0, pt);
    // std::vector<double> d2GdTdn = this->d2G_dTdni(p0, T0, n0, 0, pt);
    std::vector<double> d2HidTdn = this->d2Hi_dTdni(T0, n0, 0, pt);
    // std::vector<double> d2HrdTdn = this->d2Hr_dTdni(p0, T0, n0, 0, pt);
    // std::vector<double> d2HdTdn = this->d2H_dTdni(p0, T0, n0, 0, pt);
    std::vector<double> d2SidTdn = this->d2Si_dTdni(p0, T0, n0, 0, pt);
    // std::vector<double> d2SrdTdn = this->d2Sr_dTdni(p0, T0, n0, 0, pt);
    // std::vector<double> d2SdTdn = this->d2S_dTdni(p0, T0, n0, 0, pt);

    // Temperature derivatives
    this->solve_PT(p0, T0 + dT, n0, 0, true);

    // Derivative of ideal Gibbs free energy Gi w.r.t. temperature
    Gi1 = this->Gi(p0, T0 + dT, n0, 0, pt);
    double dGi_num = (Gi1-Gi0)/dT;
    d = std::log(std::fabs(dGidT + 1e-15)) - std::log(std::fabs(dGi_num + 1e-15));
    if (verbose || std::fabs(d) > tol) { print("dGi/dT", std::vector<double>{dGidT, dGi_num, d}); error_output++; }

    // // Derivative of residual Gibbs free energy Gr w.r.t. temperature
    // Gr1 = this->Gr(p0, T0 + dT, n0, 0, pt) * (T0+dT);
    // double dGr_num = (Gr1-Gr0)/dT;
    // d = std::log(std::fabs(dGrdT * T0 + 1e-15)) - std::log(std::fabs(dGr_num + 1e-15));
    // if (verbose || std::fabs(d) > tol) { print("dGr/dT", std::vector<double>{dGrdT * T0, dGr_num, d}); error_output++; }

    // // Derivative of Gibbs free energy G w.r.t. temperature
    // G1 = this->G(p0, T0 + dT, n0, 0, pt);
    // double dG_num = (G1-G0)/dT;
    // d = std::log(std::fabs(dGdT + 1e-15)) - std::log(std::fabs(dG_num + 1e-15));
    // if (verbose || std::fabs(d) > tol) { print("dG/dT", std::vector<double>{dGdT, dG_num, d}); error_output++; }
            
    // Derivative of ideal enthalpy Hi w.r.t. temperature
    Hi1 = this->Hi(T0 + dT, n0, 0, pt);
    double dHi_num = (Hi1-Hi0)/dT;
    d = std::log(std::fabs(dHidT + 1e-15)) - std::log(std::fabs(dHi_num + 1e-15));
    if (verbose || std::fabs(d) > tol) { print("dHi/dT", std::vector<double>{dHidT, dHi_num, d}); error_output++; }

    // // Derivative of residual enthalpy Hr w.r.t. temperature
    // Hr1 = this->Hr(p0, T0+dT, n0, 0, pt) * (T0+dT);
    // double dHr_num = (Hr1-Hr0)/dT;
    // d = std::log(std::fabs(dHrdT * T0 + 1e-15)) - std::log(std::fabs(dHr_num + 1e-15));
    // if (verbose || std::fabs(d) > tol) { print("dHr/dT", std::vector<double>{dHrdT * T0, dHr_num, d}); error_output++; }

    // // Derivative of enthalpy Hr w.r.t. temperature
    // H1 = this->H(p0, T0+dT, n0, 0, pt);
    // double dH_num = (H1-H0)/dT;
    // d = std::log(std::fabs(dHdT + 1e-15)) - std::log(std::fabs(dH_num + 1e-15));
    // if (verbose || std::fabs(d) > tol) { print("dH/dT", std::vector<double>{dHdT, dH_num, d}); error_output++; }
            
    // Derivative of residual entropy Sr w.r.t. temperature
    Si1 = this->Si(p0, T0+dT, n0, 0, pt);
    double dSi_num = (Si1-Si0)/dT;
    d = std::log(std::fabs(dSidT + 1e-15)) - std::log(std::fabs(dSi_num + 1e-15));
    if (verbose || std::fabs(d) > tol) { print("dSi/dT", std::vector<double>{dSidT, dSi_num, d}); error_output++; }

    // // Derivative of residual entropy Sr w.r.t. temperature
    // Sr1 = this->Sr(p0, T0+dT, n0, 0, pt);
    // double dSr_num = (Sr1-Sr0)/dT;
    // d = std::log(std::fabs(dSrdT + 1e-15)) - std::log(std::fabs(dSr_num + 1e-15));
    // if (verbose || std::fabs(d) > tol) { print("dSr/dT", std::vector<double>{dSrdT, dSr_num, d}); error_output++; }

    // // Derivative of enthalpy Hr w.r.t. temperature
    // S1 = this->S(p0, T0+dT, n0, 0, pt);
    // double dS_num = (S1-S0)/dT;
    // d = std::log(std::fabs(dSdT + 1e-15)) - std::log(std::fabs(dS_num + 1e-15));
    // if (verbose || std::fabs(d) > tol) { print("dS/dT", std::vector<double>{dSdT, dS_num, d}); error_output++; }

    // Second derivative with respect to temperature and composition d2/dTdni
    std::vector<double> dGidn1 = this->dGi_dni(p0, T0+dT, n0, 0, pt);
    std::vector<double> dGrdn1 = this->dGr_dni(p0, T0+dT, n0, 0, pt);
    std::vector<double> dGdn1 = this->dG_dni(p0, T0+dT, n0, 0, pt);
    std::vector<double> dHidn1 = this->dHi_dni(T0+dT, n0, 0, pt);
    std::vector<double> dHrdn1 = this->dHr_dni(p0, T0+dT, n0, 0, pt);
    std::vector<double> dHdn1 = this->dH_dni(p0, T0+dT, n0, 0, pt);
    std::vector<double> dSidn1 = this->dSi_dni(p0, T0+dT, n0, 0, pt);
    std::vector<double> dSrdn1 = this->dSr_dni(p0, T0+dT, n0, 0, pt);
    std::vector<double> dSdn1 = this->dS_dni(p0, T0+dT, n0, 0, pt);

    for (int j = 0; j < nc; j++)
    {
        // Derivative of ideal molar entropy si(PT) w.r.t. temperature
        this->solve_PT(p0, T0, n0, 0, true);
        double dsidT = this->dsi_dT(p0, T0, j, true);
        double si0 = this->si(p0, T0, j, true);
        this->solve_PT(p0, T0+dT, n0, 0, true);
        double si1 = this->si(p0, T0+dT, j, true);
        double dsi_num = (si1-si0)/dT;
        d = std::log(std::fabs(dsidT + 1e-15)) - std::log(std::fabs(dsi_num + 1e-15));
        if (verbose || std::fabs(d) > tol) { print("dsi/dT", std::vector<double>{static_cast<double>(j), dsidT, dsi_num, d}); error_output++; }

        // // Second derivative of ideal Gibbs free energy Gi w.r.t. temperature and composition
        // double d2Gi_num = (dGidn1[j]-dGidn[j])/dT;
        // d = std::log(std::fabs(d2GidTdn[j] + 1e-15)) - std::log(std::fabs(d2Gi_num + 1e-15));
        // if (verbose || std::fabs(d) > tol) { print("d2Gi/dTdnj", std::vector<double>{static_cast<double>(j), d2GidTdn[j], d2Gi_num, d}); error_output++; }

        // // Second derivative of residual Gibbs free energy Gr w.r.t. temperature and composition
        // double d2Gr_num = (dGrdn1[j]-dGrdn[j])/dT;
        // d = std::log(std::fabs(d2GrdTdn[j] + 1e-15)) - std::log(std::fabs(d2Gi_num + 1e-15));
        // if (verbose || std::fabs(d) > tol) { print("d2Gr/dTdnj", std::vector<double>{static_cast<double>(j), d2GrdTdn[j], d2Gr_num, d}); error_output++; }

        // // Second derivative of total Gibbs free energy G w.r.t. temperature and composition
        // double d2G_num = (dGdn1[j]-dGdn[j])/dT;
        // d = std::log(std::fabs(d2GdTdn[j] + 1e-15)) - std::log(std::fabs(d2Gi_num + 1e-15));
        // if (verbose || std::fabs(d) > tol) { print("d2G/dTdnj", std::vector<double>{static_cast<double>(j), d2GdTdn[j], d2G_num, d}); error_output++; }

        // Second derivative of ideal enthalpy Hi w.r.t. temperature and composition
        double d2Hi_num = (dHidn1[j]-dHidn[j])/dT;
        d = std::log(std::fabs(d2HidTdn[j] + 1e-15)) - std::log(std::fabs(d2Hi_num + 1e-15));
        if (verbose || std::fabs(d) > tol) { print("d2Hi/dTdnj", std::vector<double>{static_cast<double>(j), d2HidTdn[j], d2Hi_num, d}); error_output++; }

        // // Second derivative of residual enthalpy Hr w.r.t. temperature and composition
        // double d2Hr_num = (dHrdn1[j]-dGrdn[j])/dT;
        // d = std::log(std::fabs(d2HrdTdn[j] + 1e-15)) - std::log(std::fabs(d2Hr_num + 1e-15));
        // if (verbose || std::fabs(d) > tol) { print("d2Hr/dTdnj", std::vector<double>{static_cast<double>(j), d2HrdTdn[j], d2Hr_num, d}); error_output++; }

        // // Second derivative of total enthalpy H w.r.t. temperature and composition
        // double d2H_num = (dGdn1[j]-dGdn[j])/dT;
        // d = std::log(std::fabs(d2HdTdn[j] + 1e-15)) - std::log(std::fabs(d2Gi_num + 1e-15));
        // if (verbose || std::fabs(d) > tol) { print("d2H/dTdnj", std::vector<double>{static_cast<double>(j), d2HdTdn[j], d2H_num, d}); error_output++; }

        // Second derivative of ideal entropy Si w.r.t. temperature and composition
        double d2Si_num = (dSidn1[j]-dSidn[j])/dT;
        d = std::log(std::fabs(d2SidTdn[j] + 1e-15)) - std::log(std::fabs(d2Si_num + 1e-15));
        if (verbose || std::fabs(d) > tol) { print("d2Si/dTdnj", std::vector<double>{static_cast<double>(j), d2SidTdn[j], d2Si_num, d}); error_output++; }

        // // Second derivative of residual entropy Sr w.r.t. temperature and composition
        // double d2Sr_num = (dSrdn1[j]-dSrdn[j])/dT;
        // d = std::log(std::fabs(d2SrdTdn[j] + 1e-15)) - std::log(std::fabs(d2Sr_num + 1e-15));
        // if (verbose || std::fabs(d) > tol) { print("d2Sr/dTdnj", std::vector<double>{static_cast<double>(j), d2SrdTdn[j], d2Sr_num, d}); error_output++; }

        // // Second derivative of total entropy S w.r.t. temperature and composition
        // double d2S_num = (dSdn1[j]-dSdn[j])/dT;
        // d = std::log(std::fabs(d2SdTdn[j] * T0 + 1e-15)) - std::log(std::fabs(d2S_num + 1e-15));
        // if (verbose || std::fabs(d) > tol) { print("d2S/dTdnj", std::vector<double>{static_cast<double>(j), d2SdTdn[j], d2S_num, d}); error_output++; }
    }

    // Composition derivatives
    for (int j = 0; j < nc; j++)
    {
        if (n0[j] != 0.)
        {
            std::vector<double> nn = n0;
            double dnj = 1e-5 * n0[j];
            nn[j] += dnj;
            this->solve_PT(p0, T0, nn, 0, true);

            // Derivative of ideal Gibbs free energy Gi w.r.t. composition
            Gi1 = this->Gi(p0, T0, nn, 0, pt);
            dGi_num = (Gi1-Gi0)/dnj;
            d = std::log(std::fabs(dGidn[j] + 1e-15)) - std::log(std::fabs(dGi_num + 1e-15));
            if (verbose || std::fabs(d) > tol) { print("dGi/dnj", std::vector<double>{static_cast<double>(j), dGidn[j], dGi_num, d}); error_output++; }

            // Derivative of residual Gibbs free energy Gr w.r.t. composition
            Gr1 = this->Gr(p0, T0, nn, 0, pt);
            double dGr_num = (Gr1-Gr0)/dnj;
            d = std::log(std::fabs(dGrdn[j] + 1e-15)) - std::log(std::fabs(dGr_num + 1e-15));
            if (verbose || std::fabs(d) > tol) { print("dGr/dnj", std::vector<double>{static_cast<double>(j), dGrdn[j], dGr_num, d}); error_output++; }

            // Derivative of Gibbs free energy G w.r.t. composition
            G1 = this->G(p0, T0, nn, 0, pt);
            double dG_num = (G1-G0)/dnj;
            d = std::log(std::fabs(dGdn[j] + 1e-15)) - std::log(std::fabs(dG_num + 1e-15));
            if (verbose || std::fabs(d) > tol) { print("dG/dnj", std::vector<double>{static_cast<double>(j), dGdn[j], dG_num, d}); error_output++; }
            
            // Derivative of ideal enthalpy Hi w.r.t. composition
            Hi1 = this->Hi(T0, nn, 0, pt);
            dHi_num = (Hi1-Hi0)/dnj;
            d = std::log(std::fabs(dHidn[j] + 1e-15)) - std::log(std::fabs(dHi_num + 1e-15));
            if (verbose || std::fabs(d) > tol) { print("dHi/dnj", std::vector<double>{static_cast<double>(j), dHidn[j], dHi_num, d}); error_output++; }

            // Derivative of residual enthalpy Hr w.r.t. composition
            Hr1 = this->Hr(p0, T0, nn, 0, pt);
            double dHr_num = (Hr1-Hr0)/dnj;
            d = std::log(std::fabs(dHrdn[j] + 1e-15)) - std::log(std::fabs(dHr_num + 1e-15));
            if (verbose || std::fabs(d) > tol) { print("dHr/dnj", std::vector<double>{static_cast<double>(j), dHrdn[j], dHr_num, d}); error_output++; }

            // Derivative of enthalpy Hr w.r.t. composition
            H1 = this->H(p0, T0, nn, 0, pt);
            double dH_num = (H1-H0)/dnj;
            d = std::log(std::fabs(dHdn[j] + 1e-15)) - std::log(std::fabs(dH_num + 1e-15));
            if (verbose || std::fabs(d) > tol) { print("dH/dnj", std::vector<double>{static_cast<double>(j), dHdn[j], dH_num, d}); error_output++; }
            
            // Derivative of residual entropy Sr w.r.t. composition
            Si1 = this->Si(p0, T0, nn, 0, pt);
            dSi_num = (Si1-Si0)/dnj;
            d = std::log(std::fabs(dSidn[j] + 1e-15)) - std::log(std::fabs(dSi_num + 1e-15));
            if (verbose || std::fabs(d) > tol) { print("dSi/dnj", std::vector<double>{static_cast<double>(j), dSidn[j], dSi_num, d}); error_output++; }

            // Derivative of residual entropy Sr w.r.t. composition
            Sr1 = this->Sr(p0, T0, nn, 0, pt);
            double dSr_num = (Sr1-Sr0)/dnj;
            d = std::log(std::fabs(dSrdn[j] + 1e-15)) - std::log(std::fabs(dSr_num + 1e-15));
            if (verbose || std::fabs(d) > tol) { print("dSr/dnj", std::vector<double>{static_cast<double>(j), dSrdn[j], dSr_num, d}); error_output++; }

            // Derivative of enthalpy Hr w.r.t. composition
            S1 = this->S(p0, T0, nn, 0, pt);
            double dS_num = (S1-S0)/dnj;
            d = std::log(std::fabs(dSdn[j] + 1e-15)) - std::log(std::fabs(dS_num + 1e-15));
            if (verbose || std::fabs(d) > tol) { print("dS/dnj", std::vector<double>{static_cast<double>(j), dSdn[j], dS_num, d}); error_output++; }
        }
    }

    return error_output;
}

int EoS::derivatives_test(double p_, double T_, std::vector<double>& n_, double tol, bool verbose)
{
    // Default implementation of derivatives_test() tests nothing and returns 0
    (void) p_;
    (void) T_;
    (void) n_;
    (void) tol;
    (void) verbose;
    return 0;
}
