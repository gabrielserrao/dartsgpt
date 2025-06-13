#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <map>

#include "dartsflash/maths/maths.hpp"
#include "dartsflash/maths/geometry.hpp"
#include "dartsflash/flash/flash_params.hpp"
#include "dartsflash/flash/flash_results.hpp"
#include "dartsflash/flash/flash.hpp"

FlashResults::FlashResults(FlashParams& flashparams, double p_, double T_, std::vector<double>& z_, std::vector<TrialPhase>& comps, bool derivs)
{
    // Return p and T values
    this->flash_params = flashparams;
    this->pressure = p_;
    this->temperature = T_;
    this->zi = z_;

    // Phase fractions, phase compositions, eos and roots
    np = 0;
    int j = 0;
    for (std::string eos_name: flashparams.eos_order)
    {
        EoSParams* eos_params = &flashparams.eos_params[eos_name];
        for (EoS::RootFlag root: eos_params->root_order)
        {
            // Find compositions of particular EoS and root type
            bool phase_present = false;
            std::vector<int> comp_order = {};
            for (int jj = 0; jj < static_cast<int>(comps.size()); jj++)
            {
                // If eos_name and root type correspond, add jj'th comp to idxs
                if (comps[jj].eos_name == eos_name && comps[jj].nu > 0. && (root == EoS::RootFlag::STABLE || comps[jj].root == root ||
                                                                           (comps[jj].root == EoS::RootFlag::STABLE && eos_params->eos->is_root_type(comps[jj].X.begin()) == root)))
                {
                    phase_present = true;
                    comp_order.push_back(jj);
                }
            }
            // If phase is not present, add -1 idx of length of comp order or of length 1 if comp order has not been specified
            if (!phase_present)
            {
                comp_order = std::vector<int>{-1};
            }

            // If root type is liquid and rich phase order has been specified, put rich-phases in the right order
            if (root == EoS::RootFlag::MIN && !eos_params->rich_phase_order.empty())
            {
                // If no liquid phase is present, comp_order is set to {-1}; clear it
                if (comp_order.size() == 1 && comp_order[0] == -1)
                {
                    comp_order.clear();
                }

                // Loop over rich phases to find the right order
                size_t ii = 0;  // ii'th rich phase
                for (int i: eos_params->rich_phase_order)
                {
                    size_t idx = (i == -1) ? ii : -1;
                    double Xi = NAN;

                    // If i'th rich phase is specified
                    if (i > -1)
                    {
                        // Find if there is a i'th-component rich phase
                        for (size_t jj = ii; jj < comp_order.size(); jj++)
                        {
                            double Xj = comps[comp_order[jj]].X[i];
                            // If Xj is rich phase, set idx and continue loop to find max
                            if (Xj >= eos_params->rich_phase_composition && (std::isnan(Xi) || Xj > Xi))
                            {
                                idx = jj;
                                Xi = Xj;
                            }
                        }
                    }
                    // Else, if non-rich phase is specified, check if there are phases left in the set of compositions
                    else
                    {
                        if (ii < comp_order.size())
                        {
                            idx = ii;
                            Xi = 1.;
                        }
                    }

                    // If the (rich) phase exists, put it to the right location in comp_order
                    if (!std::isnan(Xi))
                    {
                        std::rotate(comp_order.begin() + ii, comp_order.begin() + idx, comp_order.end());
                    }
                    // Else, rich phase is not present, add -1 index in the right place
                    else
                    {
                        comp_order.insert(comp_order.begin() + ii, -1);
                    }

                    ii++;
                }
            }

            // Add jj compositions to the results vectors
            for (int jj: comp_order)
            {
                // If (rich) phase is present, it has been assigned a non-negative index
                if (jj >= 0)
                {
                    double nu_j = comps[jj].nu;
                    nuj.push_back(nu_j);
                    Xij.insert(Xij.end(), comps[jj].X.begin(), comps[jj].X.end());
                    nij.resize(Xij.size());
                    std::transform(Xij.end() - flashparams.ns, Xij.end(), nij.end() - flashparams.ns, 
                                        [&nu_j](double element) { return element *= nu_j; });
                    eos_idx.push_back(j);
                    root_type.push_back(root);
                    phase_idxs.push_back(static_cast<int>(nuj.size()-1));  // global index of jj'th phase type);
                    np++;
                }
                // If not present, they have been assigned a -1 index
                else
                {
                    nuj.push_back(0.);
                    std::vector<double> Xnan(flashparams.ns, NAN);
                    Xij.insert(Xij.end(), Xnan.begin(), Xnan.end());
                    nij.insert(nij.end(), Xnan.begin(), Xnan.end());
                    eos_idx.push_back(j);
                    root_type.push_back(root);
                }
            }
        }
        j++;
    }
    this->np_tot = static_cast<int>(this->nuj.size());

    // If needed, get derivatives
    if (derivs)
    {
        // Solve linear system to obtain derivatives of flash for simulation
        // Simulation requires derivatives of phase fractions and phase compositions w.r.t. primary variables X
        // Need to apply chain rule on derivatives w.r.t. pressure, temperature and composition
        this->calc_matrix_and_inverse();
    }
}

double FlashResults::phase_prop(EoS::Property prop, int phase_idx)
{
    // Calculate phase molar property
    EoS* eosj = flash_params.eos_params[flash_params.eos_order[eos_idx[phase_idx]]].eos;
    if (prop == EoS::Property::ENTROPY)
    {
        return eosj->S(this->pressure, this->temperature, this->Xij, phase_idx * flash_params.ns, true) * M_R;
    }
    else if (prop == EoS::Property::GIBBS)
    {
        return eosj->G(this->pressure, this->temperature, this->Xij, phase_idx * flash_params.ns, true) * M_R;
    }
    else // if (prop == EoS::Property::ENTHALPY)
    {
        return eosj->H(this->pressure, this->temperature, this->Xij, phase_idx * flash_params.ns, true) * M_R;
    }
}
double FlashResults::total_prop(EoS::Property prop)
{
    // Calculate total molar property
    double result = 0.;
    for (int j = 0; j < np; j++)
    {
        result += this->nuj[phase_idxs[j]] * this->phase_prop(prop, phase_idxs[j]);
    }
    return result;
}

void FlashResults::calc_matrix_and_inverse()
{
    // Calculate generic matrix and inverse to calculate derivatives of flash equations w.r.t. state variables and composition

    // Mass balance equations: 1 - Σj nu_j = 0 [1]; lnfi0 - lnfij = 0 [(NP-1)*NC]; zi - Σj xij nu_j = 0 [NC], Σi (xi0-xij) = 0 [NP-1] -> total [(NC+1)*NP]
    // Unknowns: dnuj/dX [NP], dxij/dX [NP*NC], dnuj/dzk [NP * NC], dxij/dzk [(NP*NC) * NC]
    int nc = flash_params.nc;
    int n_eq = np * (nc+1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n_eq, n_eq);

    // Solve reference phase fugacity coefficients
    std::vector<std::vector<double>> dlnphiijdn(this->np);
    for (int j = 0; j < this->np; j++)
    {
        // Solve phase j fugacity coefficients
        std::string eosname = flash_params.eos_order[eos_idx[phase_idxs[j]]];
        EoS* eosj = flash_params.eos_params[eosname].eos;
        eosj->set_root_flag(root_type[phase_idxs[j]]);
        eosj->solve_PT(pressure, temperature, nij, phase_idxs[j]*nc, true);
        dlnphiijdn[j] = eosj->dlnphi_dn();
    }

    // Construct A and B

    // [1 eq] [0] 1 - Σj nu_j = 0
    // d/dX: dnu0/dX + Σj dnuj/dX = 0
    for (int j = 0; j < this->np; j++)
    {
        A(0, j) += 1.;
    }

    // [NP-1 eq] [1, NP] Σi (xi0-xij) = 0
    // Σi (dx_i0/dX - dx_ij/dX) = 0
    for (int j = 1; j < this->np; j++)
    {
        for (int i = 0; i < nc; i++)
        {
            A(j, np + i) += 1.;
            A(j, np + j*nc + i) -= 1.;
        }
    }

    // [NC eq] [NP, NP+NC] zi - Σj nu_j xij = 0
    // Σj x_ij dnu_j/dzk + nu_j dx_ij/dzk = dzi/dX {1-zi if i == k else -zi}
    for (int i = 0; i < nc; i++)
    {
        int idxi = np + i;
        for (int j = 0; j < this->np; j++)
        {
            // nu_j contribution: Σj x_ij
            A(idxi, j) += Xij[phase_idxs[j]*nc + i];

            // xij contribution: Σj nu_j
            A(idxi, np + j*nc + i) += nuj[phase_idxs[j]];
        }
    }

    // [(NP-1)*NC eq] [NP+NC, (NC+1)*NP] lnfi0 - lnfij = 0
    // dlnphi_i0/dX + Σk dlnphi_i0/dn_k0 dn_k0/dX + 1/xi0 dx_i0/dX - dlnphi_ij/dX - Σk dlnphi_ij/dn_kj dn_kj/dX - 1/xij dx_ij/dX  = 0
    // dn_ij/dX = dnu_j/dX x_ij + nu_j dx_ij/dX
    for (int j = 1; j < this->np; j++)
    {
        for (int i = 0; i < nc; i++)
        {
            int idxi = np + j*nc + i;

            // dnuj/dX contribution: Σk dlnphi_i0/dn_k0 x_k0 dnu_0/dX - Σk dlnphi_ij/dn_kj x_kj dnu_j/dX
            for (int k = 0; k < nc; k++)
            {
                A(idxi, 0) += dlnphiijdn[0][i*nc + k] * Xij[phase_idxs[0]*nc + k];
                A(idxi, j) -= dlnphiijdn[j][i*nc + k] * Xij[phase_idxs[j]*nc + k];
            }

            // dxi0/dX contribution: Σi dlnphi_i0/dn_k0 nu_0 dx_k0/dX + 1/xi0 dx_i0/dX
            A(idxi, this->np + i) += 1./Xij[phase_idxs[0]*nc + i];
            for (int k = 0; k < nc; k++)
            {
                A(idxi, this->np + k) += dlnphiijdn[0][i*nc + k] * nuj[phase_idxs[0]];
            }

            // dxij/dX contribution: - Σk dlnphi_ij/dn_kj nu_j dx_kj/dX - 1/xij dx_ij/dX
            A(idxi, this->np + j*nc + i) -= 1./Xij[phase_idxs[j]*nc + i];
            for (int k = 0; k < nc; k++)
            {
                A(idxi, this->np + j*nc + k) -= dlnphiijdn[j][i*nc + k] * nuj[phase_idxs[j]];
            }
        }
    }

    // Solve Ay = b for dP and dT and AY = B for dzk
    this->LUofA.compute(A);
}

void FlashResults::dP_derivs(std::vector<double>& dnudP, std::vector<double>& dxdP)
{
    // Derivatives of flash w.r.t. pressure, temperature and composition zk (mole fractions)

    // Mass balance equations: 1 - Σj nu_j = 0 [1]; lnfi0 - lnfij = 0 [(NP-1)*NC]; zi - Σj xij nu_j = 0 [NC], Σi (xi0-xij) = 0 [NP-1] -> total [(NC+1)*NP]
    // Unknowns: dnuj/dX [NP], dxij/dX [NP*NC], dnuj/dzk [NP * NC], dxij/dzk [(NP*NC) * NC]
    int nc = flash_params.nc;
    int n_eq = np * (nc+1);
    dnudP = std::vector<double>(np_tot, 0.);
    dxdP = std::vector<double>(np_tot * nc, 0.);
    Eigen::VectorXd bP = Eigen::VectorXd::Zero(n_eq);

    // Solve reference phase fugacity coefficients
    std::vector<std::vector<double>> dlnphiijdP(this->np);
    for (int j = 0; j < this->np; j++)
    {
        // Solve phase j fugacity coefficients
        std::string eosname = flash_params.eos_order[eos_idx[phase_idxs[j]]];
        EoS* eosj = flash_params.eos_params[eosname].eos;
        eosj->set_root_flag(root_type[phase_idxs[j]]);
        eosj->solve_PT(pressure, temperature, nij, phase_idxs[j]*nc, true);
        dlnphiijdP[j] = eosj->dlnphi_dP();
    }

    // Construct B

    // [1 eq] [0] 1 - Σj nu_j = 0
    // d/dX: dnu0/dX + Σj dnuj/dX = 0

    // [NP-1 eq] [1, NP] Σi (xi0-xij) = 0
    // Σi (dx_i0/dX - dx_ij/dX) = 0

    // [NC eq] [NP, NP+NC] zi - Σj nu_j xij = 0
    // Σj x_ij dnu_j/dzk + nu_j dx_ij/dzk = dzi/dX {1 if i == k else -1}

    // [(NP-1)*NC eq] [NP+NC, (NC+1)*NP] lnfi0 - lnfij = 0
    // dlnphi_i0/dX + Σk dlnphi_i0/dn_k0 dn_k0/dX + 1/xi0 dx_i0/dX - dlnphi_ij/dX - Σk dlnphi_ij/dn_kj dn_kj/dX - 1/xij dx_ij/dX  = 0
    // dn_ij/dX = dnu_j/dX x_ij + nu_j dx_ij/dX
    for (int j = 1; j < this->np; j++)
    {
        for (int i = 0; i < nc; i++)
        {
            int idxi = np + j*nc + i;

            // bi contribution: -(dlnphi_i0/dX - dlnphi_ij/dX)
            bP(idxi) -= dlnphiijdP[0][i] - dlnphiijdP[j][i];
        }
    }

    // Solve Ay = b for dP and dT and AY = B for dzk
    Eigen::VectorXd yP = this->LUofA.solve(bP);

    for (int j = 0; j < np; j++)
    {
        // P and T derivatives
        dnudP[phase_idxs[j]] = yP(j);
        for (int i = 0; i < nc; i++)
        {
            dxdP[phase_idxs[j]*nc + i] = yP(np + j*nc + i);
        }
    }
    return;
}

void FlashResults::dT_derivs(std::vector<double>& dnudT, std::vector<double>& dxdT)
{
    // Derivatives of flash w.r.t. pressure, temperature and composition zk (mole fractions)

    // Mass balance equations: 1 - Σj nu_j = 0 [1]; lnfi0 - lnfij = 0 [(NP-1)*NC]; zi - Σj xij nu_j = 0 [NC], Σi (xi0-xij) = 0 [NP-1] -> total [(NC+1)*NP]
    // Unknowns: dnuj/dX [NP], dxij/dX [NP*NC], dnuj/dzk [NP * NC], dxij/dzk [(NP*NC) * NC]
    int nc = flash_params.nc;
    int n_eq = np * (nc+1);
    dnudT = std::vector<double>(np_tot, 0.);
    dxdT = std::vector<double>(np_tot * nc, 0.);
    Eigen::VectorXd bT = Eigen::VectorXd::Zero(n_eq);

    // Solve reference phase fugacity coefficients
    std::vector<std::vector<double>> dlnphiijdT(this->np);
    for (int j = 0; j < this->np; j++)
    {
        // Solve phase j fugacity coefficients
        std::string eosname = flash_params.eos_order[eos_idx[phase_idxs[j]]];
        EoS* eosj = flash_params.eos_params[eosname].eos;
        eosj->set_root_flag(root_type[phase_idxs[j]]);
        eosj->solve_PT(pressure, temperature, nij, phase_idxs[j]*nc, true);
        dlnphiijdT[j] = eosj->dlnphi_dT();
    }

    // Construct A and B

    // [1 eq] [0] 1 - Σj nu_j = 0
    // d/dX: dnu0/dX + Σj dnuj/dX = 0

    // [NP-1 eq] [1, NP] Σi (xi0-xij) = 0
    // Σi (dx_i0/dX - dx_ij/dX) = 0

    // [NC eq] [NP, NP+NC] zi - Σj nu_j xij = 0
    // Σj x_ij dnu_j/dzk + nu_j dx_ij/dzk = dzi/dX {1 if i == k else -1}

    // [(NP-1)*NC eq] [NP+NC, (NC+1)*NP] lnfi0 - lnfij = 0
    // dlnphi_i0/dX + Σk dlnphi_i0/dn_k0 dn_k0/dX + 1/xi0 dx_i0/dX - dlnphi_ij/dX - Σk dlnphi_ij/dn_kj dn_kj/dX - 1/xij dx_ij/dX  = 0
    // dn_ij/dX = dnu_j/dX x_ij + nu_j dx_ij/dX
    for (int j = 1; j < this->np; j++)
    {
        for (int i = 0; i < nc; i++)
        {
            int idxi = np + j*nc + i;

            // bi contribution: -(dlnphi_i0/dX - dlnphi_ij/dX)
            bT(idxi) -= dlnphiijdT[0][i] - dlnphiijdT[j][i];
        }
    }

    // Solve Ay = b for dP and dT and AY = B for dzk
    Eigen::VectorXd yT = this->LUofA.solve(bT);
    
    for (int j = 0; j < np; j++)
    {
        // P and T derivatives
        dnudT[phase_idxs[j]] = yT(j);
        for (int i = 0; i < nc; i++)
        {
            dxdT[phase_idxs[j]*nc + i] = yT(np + j*nc + i);
        }
    }
    return;
}

void FlashResults::dz_derivs(std::vector<double>& dnudzk, std::vector<double>& dxdzk)
{
    // Derivatives of flash w.r.t. pressure, temperature and composition zk (mole fractions)

    // Mass balance equations: 1 - Σj nu_j = 0 [1]; lnfi0 - lnfij = 0 [(NP-1)*NC]; zi - Σj xij nu_j = 0 [NC], Σi (xi0-xij) = 0 [NP-1] -> total [(NC+1)*NP]
    // Unknowns: dnuj/dX [NP], dxij/dX [NP*NC], dnuj/dzk [NP * NC], dxij/dzk [(NP*NC) * NC]
    int nc = flash_params.nc;
    int n_eq = np * (nc+1);
    dnudzk = std::vector<double>(np_tot * nc, 0.);
    dxdzk = std::vector<double>(np_tot * nc * nc, 0.);
    Eigen::MatrixXd Bz = Eigen::MatrixXd::Zero(n_eq, nc);

    // Solve reference phase fugacity coefficients
    std::vector<std::vector<double>> dlnphiijdP(this->np), dlnphiijdT(this->np), dlnphiijdn(this->np);
    for (int j = 0; j < this->np; j++)
    {
        // Solve phase j fugacity coefficients
        std::string eosname = flash_params.eos_order[eos_idx[phase_idxs[j]]];
        EoS* eosj = flash_params.eos_params[eosname].eos;
        eosj->set_root_flag(root_type[phase_idxs[j]]);
        eosj->solve_PT(pressure, temperature, nij, phase_idxs[j]*nc, true);
        dlnphiijdn[j] = eosj->dlnphi_dn();
    }

    // Construct A and B

    // [1 eq] [0] 1 - Σj nu_j = 0
    // d/dX: dnu0/dX + Σj dnuj/dX = 0

    // [NP-1 eq] [1, NP] Σi (xi0-xij) = 0
    // Σi (dx_i0/dX - dx_ij/dX) = 0

    // [NC eq] [NP, NP+NC] zi - Σj nu_j xij = 0
    // Σj x_ij dnu_j/dzk + nu_j dx_ij/dzk = dzi/dX {1-zi if i == k else -zk}
    for (int i = 0; i < nc; i++)
    {
        int idxi = np + i;

        // Bik contribution: {1-zi if i == k else -zk}
        for (int k = 0; k < nc; k++)
        {
            Bz(idxi, k) = -zi[i];
        }
        Bz(idxi, i) = 1.-zi[i];
    }

    // [(NP-1)*NC eq] [NP+NC, (NC+1)*NP] lnfi0 - lnfij = 0
    // dlnphi_i0/dX + Σk dlnphi_i0/dn_k0 dn_k0/dX + 1/xi0 dx_i0/dX - dlnphi_ij/dX - Σk dlnphi_ij/dn_kj dn_kj/dX - 1/xij dx_ij/dX  = 0
    // dn_ij/dX = dnu_j/dX x_ij + nu_j dx_ij/dX

    // Solve Ay = b for dP and dT and AY = B for dzk
    Eigen::MatrixXd Yz = LUofA.solve(Bz);

    for (int j = 0; j < np; j++)
    {
        // Composition derivatives
        for (int k = 0; k < nc; k++)
        {
            dnudzk[k * np_tot + phase_idxs[j]] = Yz(j, k);
            for (int i = 0; i < nc; i++)
            {
                dxdzk[k * np_tot * nc + phase_idxs[j]*nc + i] = Yz(np + j*nc + i, k);
            }
        }
    }
    return;
}

double FlashResults::dXdT(EoS::Property prop, std::vector<double>& dnudT, std::vector<double>& dxdT)
{
    // Calculate derivative of total enthalpy or entropy with respect to temperature

    // Enthalpy: dH/dT = d/dT [Σj H_j] = Σj dHj/dT
    //                 = Σj [(dHj/dT)_P,n + Σk (dHj/dnk)_P,T dnk/dT] = Σj [Cpj + Σk dHj/dnk (nuj dxk/dT + xk dnuj/dT]
    if (prop == EoS::Property::ENTHALPY)
    {
        double dHdT = 0.;
        for (int j = 0; j < np; j++)
        {
            EoS* eosj = flash_params.eos_params[flash_params.eos_order[eos_idx[phase_idxs[j]]]].eos;
            int start_idx = phase_idxs[j]*flash_params.ns;

            // Derivative of enthalpy with respect to temperature (dHj/dT)_P,n: Heat capacity Cpr
            dHdT += eosj->Cp(pressure, temperature, nij, start_idx, true) * M_R;

            // Derivative of enthalpy with respect to composition (dHj/dnk)_P,T
            std::vector<double> dHj_dnk = eosj->dH_dni(pressure, temperature, nij, start_idx, true);
            for (int i = 0; i < flash_params.nc; i++)
            {
                dHdT += dHj_dnk[i] * (dnudT[phase_idxs[j]] * Xij[start_idx + i] 
                                    + nuj[phase_idxs[j]] * dxdT[start_idx + i]) * M_R;
            }
        }
        return dHdT;
    }
    // Entropy: dS/dT = d/dT [Σj Sj] = Σj dSj/dT
    //                = Σj [(dSj/dT)_P,n + Σk (dSj/dnk)_P,T dnk/dT] = Σj [Cpj/T + Σk dSj/dnk (nuj dxk/dT + xk dnuj/dT]
    else
    {
        double dSdT = 0.;
        for (int j = 0; j < np; j++)
        {
            EoS* eosj = flash_params.eos_params[flash_params.eos_order[eos_idx[phase_idxs[j]]]].eos;
            int start_idx = phase_idxs[j]*flash_params.ns;

            // Derivative of entropy with respect to temperature (dSj/dT)_P,n: Heat capacity Cpr/T
            dSdT += eosj->Cp(pressure, temperature, nij, start_idx, true) * M_R / temperature;

            // Derivative of enthalpy with respect to composition (dHj/dnk)_P,T
            std::vector<double> dSj_dnk = eosj->dS_dni(pressure, temperature, nij, start_idx, true);
            for (int i = 0; i < flash_params.nc; i++)
            {
                dSdT += dSj_dnk[i] * (dnudT[phase_idxs[j]] * Xij[start_idx + i] 
                                    + nuj[phase_idxs[j]] * dxdT[start_idx + i]) * M_R;
            }
        }
        return dSdT;
    }
}

void FlashResults::dX_derivs(EoS::Property prop, std::vector<double>& dnudT, std::vector<double>& dxdT, 
                                                 std::vector<double>& dnudX, std::vector<double>& dxdX)
{
    // Derivatives of flash with respect to other state specifications: H, S
    dnudX = dnudT;
    dxdX = dxdT;
    
    // Enthalpy or entropy: dnu/dX = dnu/dT dT/dX
    //                      dx/dX = dx/dT dT/dX
    double dTdX = 1./this->dXdT(prop, dnudT, dxdT);
    
    std::transform(dnudT.begin(), dnudT.end(), dnudX.begin(),
                   [&dTdX](double element) { return element *= dTdX; });
    std::transform(dxdT.begin(), dxdT.end(), dxdX.begin(),
                   [&dTdX](double element) { return element *= dTdX; });
    return;
}
