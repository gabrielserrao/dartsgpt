//--------------------------------------------------------------------------
#ifndef OPENDARTS_FLASH_FLASH_FLASHRESULTS_H
#define OPENDARTS_FLASH_FLASH_FLASHRESULTS_H
//--------------------------------------------------------------------------

#include <vector>
#include "dartsflash/global/global.hpp"
#include "dartsflash/flash/flash_params.hpp"
#include "dartsflash/flash/flash_results.hpp"
#include "dartsflash/flash/trial_phase.hpp"

struct FlashResults
{
	int np, np_tot;
	double pressure, temperature;
	std::vector<double> zi, nuj, Xij, nij;
	std::vector<int> eos_idx, phase_idxs;
	std::vector<EoS::RootFlag> root_type;
	Eigen::PartialPivLU<Eigen::MatrixXd> LUofA;
	FlashParams flash_params;

	FlashResults(FlashParams& flashparams,
				 double p_, double T_, std::vector<double>& z_, std::vector<TrialPhase>& comps, bool derivs);

	// Phase and total properties
	double phase_prop(EoS::Property prop, int phase_idx);  // phase property
	double total_prop(EoS::Property prop);  // total property

	// Derivatives
	void calc_matrix_and_inverse();

	// Primary flash variables: P, T, z
	void dP_derivs(std::vector<double>& dnudP, std::vector<double>& dxdP);
	void dT_derivs(std::vector<double>& dnudT, std::vector<double>& dxdT);
	void dz_derivs(std::vector<double>& dnudzk, std::vector<double>& dxdzk);

	// Secondary flash variables and derivatives of flash with respect to other specifications: H, S
	double dXdT(EoS::Property prop, std::vector<double>& dnudT, std::vector<double>& dxdT);
	void dX_derivs(EoS::Property prop, std::vector<double>& dnudT, std::vector<double>& dxdT, 
									   std::vector<double>& dnudX, std::vector<double>& dxdX);

	void print_results()
	{
		print("p, T", std::vector<double>{pressure, temperature});
		print("nu", nuj);
		print("X", Xij, np);
		print("eos_idx", eos_idx);
		print("roots", root_type);
	}
};

//--------------------------------------------------------------------------
#endif // OPENDARTS_FLASH_FLASH_FLASHRESULTS_H
//--------------------------------------------------------------------------
